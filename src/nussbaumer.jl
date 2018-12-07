"""
Polynomial multiplication algorithm by
H.J. Nussbaumer, "Fast Polynomial Transform Algorithms for Digital Convolution"
IEEE Transactions on Acoustics, Speech, and Signal Processing, 28(2), 205–215 (1980)
doi:10.1109/TASSP.1980.1163372

The algorithm in a more clear form is provided in Knuth's TAOCP Vol.2, Exercise 4.6.4-59.
"""

struct NussbaumerStage
    log_len :: Int
    r :: Int
    m :: Int
    batch :: Int
    log_r :: Int
    log_m :: Int
end


struct NussbaumerPlan{T <: Integer}
    stages :: Array{NussbaumerStage, 1}
    buffers :: Array{Array{T, 1}, 1}
    output1 :: Array{T, 1}
    output2 :: Array{T, 1}
    shift_buffer :: Array{T, 1}
    scale :: T
    inverse_scale :: Bool
    final_batch :: Int

    function NussbaumerPlan(tp::Type{<: Integer}, len::Int, negacyclic::Bool)
        log_len = trailing_zeros(len)
        @assert 2^log_len == len
        @assert len >= 4 # The last stage multiplies length-4 polynomials.

        if log_len > 2
            log_lens = [log_len]
            while true
                new_log_len = cld(log_lens[end], 2)
                if new_log_len == 2
                    break
                end
                push!(log_lens, new_log_len)
            end
        else
            log_lens = []
        end

        stages = NussbaumerStage[]
        buffers = Array{tp, 1}[]
        for (i, log_len) in enumerate(log_lens)
            m = 2^fld(log_len, 2)
            r = 2^cld(log_len, 2)
            push!(stages, NussbaumerStage(log_len, r, m, len * 2^i ÷ (r * 2m),
                cld(log_len, 2), fld(log_len, 2)))

            # output1 and output2 are used as the buffers for the last stage
            if i != length(log_lens)
                push!(buffers, Array{tp}(undef, r * 2m * (len * 2^i) ÷ (r * 2m)))
            end
        end

        println(len, ": ", [stage.r for stage in stages])

        final_batch = len * 2^(length(stages)) ÷ 4

        if length(stages) > 0
            output1 = Array{tp}(undef, len * 2^length(stages))
            output2 = Array{tp}(undef, len * 2^length(stages))
            shift_buffer = Array{tp}(undef, stages[1].r)
        else
            output1 = Array{tp}(undef, len)
            output2 = Array{tp}(undef, len)
            shift_buffer = Array{tp}(undef, 0)
        end

        # Since we did not divide by 2 in the internal multiplication functions
        # (because it may be slow for some residue ring element representations),
        # we need to rescale here.
        log_scale = get_log_scale(log_len, negacyclic)

        scale = 1 << log_scale
        if tp <: AbstractRRElem
            # Assuming here that gcd(scale, modulus) == 1
            # Since scale is a power of 2, it is enough for the modulus to be odd.
            scale = invmod(scale, rr_modulus_simple(tp))
            inv_scale = true
        else
            inv_scale = false
        end

        new{tp}(stages, buffers, output1, output2, shift_buffer, scale, inv_scale, final_batch)
    end
end


const _nussbaumer_plans = Dict{Tuple{Type, Int, Bool}, NussbaumerPlan}()


function get_nussbaumer_plan(tp::Type{<: Integer}, len::Int, negacyclic::Bool)
    key = (tp, len, negacyclic)
    if !haskey(_nussbaumer_plans, key)
        plan = NussbaumerPlan(tp, len, negacyclic)
        _nussbaumer_plans[key] = plan
        plan
    else
        _nussbaumer_plans[key]
    end
end


@inline function nussbaumer_mul_negacyclic(x::Array{T, 1}, y::Array{T, 1}, rescale::Bool) where T
    plan = get_nussbaumer_plan(T, length(x), true)
    res = similar(x)
    nussbaumer_negacyclic_forward!(plan, plan.output1, x)
    nussbaumer_negacyclic_forward!(plan, plan.output2, y)
    nussbaumer_negacyclic_mul_4!(plan.output1, plan.output1, plan.output2, plan.final_batch)
    nussbaumer_negacyclic_inverse!(plan, res, plan.output1)

    if rescale
        if plan.inverse_scale
            res .*= plan.scale
        else
            res .÷= plan.scale
        end
    end

    res
end


@inline function nussbaumer_mul_cyclic(x::Array{T, 1}, y::Array{T, 1}, rescale::Bool) where T
    # TODO: since cyclic convolution is not of interest at the moment,
    # a slower recursive version is implemented.
    # Most probably can be sped up.

    # assuming that the polynomial length is power of 2
    n = trailing_zeros(length(x))

    if n == 2
        z = similar(x)
        nussbaumer_cyclic_mul_4!(z, x, y, 1)
        return z
    end

    m = length(x) >> 1

    xx = copy(x)
    yy = copy(y)

    for k in 1:m
        t = x[m + k]
        xx[m + k] = x[k] - t
        xx[k] = x[k] + t

        t = y[m + k]
        yy[m + k] = y[k] - t
        yy[k] = y[k] + t
    end

    z = similar(x)
    z[1:m] = nussbaumer_mul_cyclic(xx[1:m], yy[1:m], false)
    z[m+1:2m] = nussbaumer_mul_negacyclic(xx[m+1:2m], yy[m+1:2m], false)

    for k in 1:m
        t = z[m + k]
        z[m + k] = (z[k] - t)
        z[k] = (z[k] + t)
    end

    # To avoid final rescaling `z` must be divided by 2 here.
    # Instead, in order for the accumulated scaling from cyclic multiplication for the length `2^n`
    # to be the same as for the negacyclic multiplication, we add a multiplication by 2
    # for certain values of `n`.
    if 2^trailing_zeros(n - 1) == n - 1
        z .*= 2
    end

    if rescale
        plan = get_nussbaumer_plan(T, length(x), true)
        if plan.inverse_scale
            z .*= plan.scale
        else
            z .÷= plan.scale
        end
    end

    z
end


@inline function bitreverse32(x::UInt32)
    x = ((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1)
    x = ((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2)
    x = ((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4)
    x = ((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8)
    (x >> 16) | (x << 16)
end


@inline function nussbaumer_negacyclic_forward!(
        plan::NussbaumerPlan, res::Array{T, 1}, x::Array{T, 1}) where T

    if length(plan.stages) == 0
        res .= x
        return
    end

    @inbounds @simd for i in 1:length(plan.stages)

        stage = plan.stages[i]

        n = stage.log_len
        m = stage.m
        r = stage.r

        log_m = stage.log_m
        log_r = stage.log_r
        log_mr = stage.log_m + stage.log_r
        log_r_mask = 1<<log_r - 1

        batch = stage.batch

        if i == 1
            input = x
        else
            input = plan.buffers[i-1]
        end
        inp = input

        if i == length(plan.stages)
            X = res
        else
            X = plan.buffers[i]
        end

        Xt = plan.shift_buffer

        # Transpose and duplicate (the latter allows us to skip one FFT stage)
        @simd for b in 0:batch-1
            @simd for j in 0:m-1
                @simd for i in 0:r-1
                    t = inp[j + i<<log_m + b<<log_mr + 1]
                    X[i + j<<log_r + b<<(log_mr+1) + 1] = t
                    X[i + (j+m)<<log_r + b<<(log_mr+1) + 1] = t
                end
            end
        end

        @simd for j in log_m-1:-1:0
            @simd for st in 0:m-1
                s = (st >> j) << (j + 1)
                t = st & ((1 << j) - 1)

                sp = bitreverse32(UInt32(s)) >> (32 - log_m - 1 - j) # Remove hardcoding

                k = sp << (log_r - log_m)

                cycle = isodd(k >> log_r)
                k = k & log_r_mask

                i1 = s + t
                i2 = i1 + 1 << j

                #=
                Shift and linearly combine:
                t = mul_by_root_of_one(X[:,i2,:], k)
                X[:,i2,:] = X[:,i1,:] - t
                X[:,i1,:] = X[:,i1,:] + t
                Where i1=st_, i2=st_+2^j, and the multipliction by the root of one
                is just the negacyclic shift.
                =#

                @simd for b in 0:batch-1
                    i1_offset = i1<<log_r + b<<(log_mr+1) + 1
                    i2_offset = i2<<log_r + b<<(log_mr+1) + 1

                    first_offset = r - k + 1
                    last_offset = -k + 1

                    @simd for q in 0:r-1
                        Xt[q + 1] = X[q + i2_offset]
                    end

                    if cycle
                        @simd for q in 0:k-1
                            t = X[q + i1_offset]
                            t1 = Xt[q + first_offset]
                            X[q + i2_offset] = t - t1
                            X[q + i1_offset] = t + t1
                        end
                        @simd for q in k:r-1
                            t = X[q + i1_offset]
                            t1 = Xt[q + last_offset]
                            X[q + i2_offset] = t + t1
                            X[q + i1_offset] = t - t1
                        end
                    else
                        @simd for q in 0:k-1
                            t = X[q + i1_offset]
                            t1 = Xt[q + first_offset]
                            X[q + i2_offset] = t + t1
                            X[q + i1_offset] = t - t1
                        end
                        @simd for q in k:r-1
                            t = X[q + i1_offset]
                            t1 = Xt[q + last_offset]
                            X[q + i2_offset] = t - t1
                            X[q + i1_offset] = t + t1
                        end
                    end
                end
            end
        end
    end
end


@inline function nussbaumer_negacyclic_inverse!(
        plan::NussbaumerPlan, res::Array{T, 1}, x::Array{T, 1}) where T

    if length(plan.stages) == 0
        res .= x
        return
    end

    @inbounds for i in length(plan.stages):-1:1
        stage = plan.stages[i]

        if i == length(plan.stages)
            input = x
        else
            input = plan.buffers[i]
        end

        if i == 1
            output = res
        else
            output = plan.buffers[i-1]
        end

        Zt = plan.shift_buffer

        data_size = length(input)
        n = stage.log_len
        m = stage.m
        r = stage.r

        log_m = stage.log_m
        log_r = stage.log_r
        log_mr = stage.log_m + stage.log_r
        log_r_mask = 1<<log_r - 1

        batch = stage.batch

        Z = input
        Z_new = output

        @simd for j = 0:log_m
            j_mask = (1 << j) - 1
            @simd for st in 0:m-1
                s = (st >> j) << (j + 1)
                t = st & j_mask

                sp = bitreverse32(UInt32(s)) >> (32 - log_m - 1 - j) # Remove hardcoding

                k = -(sp << (log_r - log_m))

                cycle = isodd(k >> log_r)
                k = k & log_r_mask

                i1 = s + t
                i2 = i1 + 1 << j

                #=
                Shift and linearly combine:
                Z[:,i2,:] = mul_by_root_of_one(Z[:,i1,:] - Z[:,i2,:], k)
                Z[:,i1,:] = Z[:,i1,:] + Z[:,i2,:]
                Where i1=st_, i2=st_+2^j, and the multipliction by the root of one
                is just the negacyclic shift.
                =#
                @simd for b in 0:batch-1
                    @simd for q in 0:r-1
                        t = Z[q + i2<<log_r + b<<(log_mr+1) + 1]
                        Zt[q + 1] = Z[q + i1<<log_r + b<<(log_mr+1) + 1] - t
                        Z[q + i1<<log_r + b<<(log_mr+1) + 1] += t
                    end

                    if cycle
                        @simd for q in 0:k-1
                            Z[q + i2<<log_r + b<<(log_mr+1) + 1] = Zt[r-k+q + 1]
                        end
                        @simd for q in k:r-1
                            Z[q + i2<<log_r + b<<(log_mr+1) + 1] = -Zt[q-k + 1]
                        end
                    else
                        @simd for q in 0:k-1
                            Z[q + i2<<log_r + b<<(log_mr+1) + 1] = -Zt[r-k+q + 1]
                        end
                        @simd for q in k:r-1
                            Z[q + i2<<log_r + b<<(log_mr+1) + 1] = Zt[q-k + 1]
                        end
                    end
                end

                # To avoid final rescaling `Z` must be divided by 2 here.
            end
        end

        r_mask = (1 << log_r) - 1
        @simd for b in 0:batch-1
            @simd for i in 0:m-1
                @simd for j in 0:r-1
                    t1 = Z[j + i<<log_r + b<<(log_mr+1) + 1]
                    t2 = Z[(j-1) & r_mask + (i+m)<<log_r + b<<(log_mr+1) + 1]
                    if j == 0
                        t2 = -t2
                    end

                    Z_new[i + j<<log_m + b<<log_mr + 1] = t1 + t2
                end
            end
        end
    end
end


@inline function nussbaumer_cyclic_mul_4!(
        output::Array{T, 1}, x::Array{T, 1}, y::Array{T, 1}, batch::Int) where T
    # A straightforward version, by definition.
    # There are optimized low-mul versions available, but they don't seem to be working faster.
    @inbounds @simd for b in 0:batch-1
        i1 = 1 + b<<2
        i2 = i1 + 1
        i3 = i2 + 1
        i4 = i3 + 1

        a0 = x[i1]
        a1 = x[i2]
        a2 = x[i3]
        a3 = x[i4]
        b0 = y[i1]
        b1 = y[i2]
        b2 = y[i3]
        b3 = y[i4]

        z0 = a0*b0 + a1*b3 + a2*b2 + a3*b1
        z1 = a0*b1 + a1*b0 + a2*b3 + a3*b2
        z2 = a0*b2 + a1*b1 + a2*b0 + a3*b3
        z3 = a0*b3 + a1*b2 + a2*b1 + a3*b0

        output[i1] = z0
        output[i2] = z1
        output[i3] = z2
        output[i4] = z3
    end
end


@inline function nussbaumer_negacyclic_mul_4!(
        output::Array{T, 1}, x::Array{T, 1}, y::Array{T, 1}, batch::Int) where T
    # A straightforward version, by definition.
    # There are optimized low-mul versions available, but they don't seem to be working faster.
    @inbounds @simd for b in 0:batch-1
        i1 = 1 + b<<2
        i2 = i1 + 1
        i3 = i2 + 1
        i4 = i3 + 1

        a0 = x[i1]
        a1 = x[i2]
        a2 = x[i3]
        a3 = x[i4]
        b0 = y[i1]
        b1 = y[i2]
        b2 = y[i3]
        b3 = y[i4]

        z0 = a0*b0 - a1*b3 - a2*b2 - a3*b1
        z1 = a0*b1 + a1*b0 - a2*b3 - a3*b2
        z2 = a0*b2 + a1*b1 + a2*b0 - a3*b3
        z3 = a0*b3 + a1*b2 + a2*b1 + a3*b0

        output[i1] = z0
        output[i2] = z1
        output[i3] = z2
        output[i4] = z3
    end
end


function get_log_scale(log_len, negacyclic)
    if log_len <= 2
        0
    else
        if negacyclic
            fld(log_len, 2) + 1 + get_log_scale(cld(log_len, 2), true)
        else
            1 + get_log_scale(log_len - 1, true) + (2^trailing_zeros(log_len - 1) == log_len - 1)
        end
    end
end
