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


struct NussbaumerScale{T <: Integer}
    scale :: T
    inverse :: Bool

    function NussbaumerScale(tp::Type{<:Integer}, log_scale::Int)
        scale = 1 << log_scale
        if tp <: AbstractRRElem
            # Assuming here that gcd(scale, modulus) == 1
            # Since scale is a power of 2, it is enough for the modulus to be odd.
            scale = invmod(scale, rr_modulus_simple(tp))
            inverse = true
        else
            inverse = false
        end
        new{tp}(scale, inverse)
    end
end


function get_log_scale(log_len, negacyclic, log_kernel_size)
    if log_len <= log_kernel_size
        0
    else
        log_r = max(cld(log_len, 2), log_kernel_size)
        log_m = log_len - log_r

        if negacyclic
            # (log_m + 1) is the accumulated scale from inverse FFT stages
            # get_log_scale() is the scale from the recursive call
            log_m + 1 + get_log_scale(log_r, true, log_kernel_size)
        else
            # 1 is the accumulated scale from the joining stage after the recursion
            # get_log_scale() is the scale from the recursive call
            # An additional factor is needed so that
            # `get_log_scale(x, false, y) == get_log_scale(x, true, y)` for all `x` and `y`.
            # (that is, the accumulated scale is the same for cyclic and negacyclic convolutions)
            q, k = divrem(log_len - 1, log_kernel_size)
            adjust = k == 0 && 2^trailing_zeros(q) == q

            (1 + get_log_scale(log_len - 1, true, log_kernel_size) + adjust)
        end
    end
end


struct NussbaumerPlan{T <: Integer}
    stages :: Array{NussbaumerStage, 1}
    buffers :: Array{Array{T, 1}, 1}
    output1 :: Array{T, 1}
    output2 :: Array{T, 1}
    shift_buffer :: Array{T, 1}
    cyclic_scale :: NussbaumerScale{T}
    negacyclic_scale :: NussbaumerScale{T}
    final_batch :: Int
    log_kernel_size :: Int
    kernel_size_val :: Val

    function NussbaumerPlan(tp::Type{<: Integer}, len::Int)

        log_kernel_size = 3

        log_len = trailing_zeros(len)
        @assert 2^log_len == len
        @assert len >= 2^log_kernel_size

        stages = NussbaumerStage[]
        buffers = Array{tp, 1}[]
        c_log_len = log_len
        if c_log_len > log_kernel_size
            while true
                i = length(stages) + 1
                log_r = max(cld(c_log_len, 2), log_kernel_size)
                log_m = c_log_len - log_r

                m = 2^log_m
                r = 2^log_r

                push!(
                    stages,
                    NussbaumerStage(c_log_len, r, m, len * 2^i ÷ (r * 2m), log_r, log_m))

                push!(buffers, Array{tp}(undef, r * 2m * (len * 2^i) ÷ (r * 2m)))

                if log_r == log_kernel_size
                    break
                end
                c_log_len = log_r
            end

            # output1 and output2 are used as the buffers for the last stage
            pop!(buffers)
        end

        final_batch = len * 2^(length(stages)) ÷ (2^log_kernel_size)

        if length(stages) > 0
            output1 = Array{tp}(undef, len * 2^length(stages))
            output2 = Array{tp}(undef, len * 2^length(stages))
            shift_buffer = Array{tp}(undef, stages[1].r)
        else
            output1 = Array{tp}(undef, len)
            output2 = Array{tp}(undef, len)
            shift_buffer = Array{tp}(undef, 0)
        end

        #log_scale = get_log_scale(log_len, negacyclic)

        new{tp}(
            stages, buffers, output1, output2, shift_buffer,
            NussbaumerScale(tp, get_log_scale(log_len, false, log_kernel_size)),
            NussbaumerScale(tp, get_log_scale(log_len, true, log_kernel_size)),
            final_batch, log_kernel_size, Val(2^log_kernel_size))
    end
end


const _nussbaumer_plans = Dict{Tuple{Type, Int}, NussbaumerPlan}()


function get_nussbaumer_plan(tp::Type{<: Integer}, len::Int)
    key = (tp, len)
    if !haskey(_nussbaumer_plans, key)
        plan = NussbaumerPlan(tp, len)
        _nussbaumer_plans[key] = plan
        plan
    else
        _nussbaumer_plans[key]
    end
end


@inline function nussbaumer_mul_negacyclic(x::Array{T, 1}, y::Array{T, 1}, rescale::Bool) where T
    plan = get_nussbaumer_plan(T, length(x))
    res = similar(x)
    nussbaumer_negacyclic_forward!(plan, plan.output1, x)
    nussbaumer_negacyclic_forward!(plan, plan.output2, y)
    nussbaumer_negacyclic_kernel!(
        plan.output1, plan.output1, plan.output2, plan.final_batch, plan.kernel_size_val)
    nussbaumer_negacyclic_inverse!(plan, res, plan.output1)

    # Since we did not divide by 2 in the internal multiplication functions
    # (because it may be slow for some residue ring element representations),
    # we need to rescale here.
    if rescale
        if plan.negacyclic_scale.inverse
            res .*= plan.negacyclic_scale.scale
        else
            res .÷= plan.negacyclic_scale.scale
        end
    end

    res
end


@inline function nussbaumer_mul_cyclic(x::Array{T, 1}, y::Array{T, 1}, rescale::Bool) where T
    # TODO: since cyclic convolution is not of interest at the moment,
    # a slower recursive version is implemented.
    # Most probably can be sped up.
    plan = get_nussbaumer_plan(T, length(x))

    # assuming that the polynomial length is power of 2
    n = trailing_zeros(length(x))

    if length(x) == 1 << plan.log_kernel_size
        z = similar(x)
        nussbaumer_cyclic_kernel!(z, x, y, 1, plan.kernel_size_val)
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
    q, k = divrem(n - 1, plan.log_kernel_size)
    adjust = k == 0 && 2^trailing_zeros(q) == q
    if k == 0 && 2^trailing_zeros(q) == q
        z .*= 2
    end

    if rescale
        if plan.cyclic_scale.inverse
            z .*= plan.cyclic_scale.scale
        else
            z .÷= plan.cyclic_scale.scale
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



@inline function nussbaumer_negacyclic_kernel!(
        output::Array{T, 1}, x::Array{T, 1}, y::Array{T, 1},
        batch::Int, kernel_size::Val{2}) where T
    @inbounds @simd for b in 0:batch-1
        i1 = 1 + b<<1
        i2 = i1 + 1
        z1, z2 = mul2_negacyclic(
            (x[i1], x[i2]),
            (y[i1], y[i2]))
        output[i1] = z1
        output[i2] = z2
    end
end


@inline function nussbaumer_cyclic_kernel!(
        output::Array{T, 1}, x::Array{T, 1}, y::Array{T, 1},
        batch::Int, kernel_size::Val{2}) where T
    @inbounds @simd for b in 0:batch-1
        i1 = 1 + b<<1
        i2 = i1 + 1
        z1, z2 = mul2_cyclic(
            (x[i1], x[i2]),
            (y[i1], y[i2]))
        output[i1] = z1
        output[i2] = z2
    end
end


@inline function nussbaumer_negacyclic_kernel!(
        output::Array{T, 1}, x::Array{T, 1}, y::Array{T, 1},
        batch::Int, kernel_size::Val{4}) where T
    @inbounds @simd for b in 0:batch-1
        i1 = 1 + b<<2
        i2 = i1 + 1
        i3 = i2 + 1
        i4 = i3 + 1
        z1, z2, z3, z4 = mul4_negacyclic_karatsuba(
            (x[i1], x[i2], x[i3], x[i4]),
            (y[i1], y[i2], y[i3], y[i4]))
        output[i1] = z1
        output[i2] = z2
        output[i3] = z3
        output[i4] = z4
    end
end


@inline function nussbaumer_cyclic_kernel!(
        output::Array{T, 1}, x::Array{T, 1}, y::Array{T, 1},
        batch::Int, kernel_size::Val{4}) where T
    @inbounds @simd for b in 0:batch-1
        i1 = 1 + b<<2
        i2 = i1 + 1
        i3 = i2 + 1
        i4 = i3 + 1
        z1, z2, z3, z4 = mul4_cyclic_karatsuba(
            (x[i1], x[i2], x[i3], x[i4]),
            (y[i1], y[i2], y[i3], y[i4]))
        output[i1] = z1
        output[i2] = z2
        output[i3] = z3
        output[i4] = z4
    end
end


@inline function nussbaumer_negacyclic_kernel!(
        output::Array{T, 1}, x::Array{T, 1}, y::Array{T, 1},
        batch::Int, kernel_size::Val{8}) where T
    @inbounds @simd for b in 0:batch-1
        i1 = 1 + b<<3
        i2 = i1 + 1
        i3 = i2 + 1
        i4 = i3 + 1
        i5 = i4 + 1
        i6 = i5 + 1
        i7 = i6 + 1
        i8 = i7 + 1
        z1, z2, z3, z4, z5, z6, z7, z8 = mul8_negacyclic_karatsuba(
            (x[i1], x[i2], x[i3], x[i4], x[i5], x[i6], x[i7], x[i8]),
            (y[i1], y[i2], y[i3], y[i4], y[i5], y[i6], y[i7], y[i8]))
        output[i1] = z1
        output[i2] = z2
        output[i3] = z3
        output[i4] = z4
        output[i5] = z5
        output[i6] = z6
        output[i7] = z7
        output[i8] = z8
    end
end


@inline function nussbaumer_cyclic_kernel!(
        output::Array{T, 1}, x::Array{T, 1}, y::Array{T, 1},
        batch::Int, kernel_size::Val{8}) where T
    @inbounds @simd for b in 0:batch-1
        i1 = 1 + b<<3
        i2 = i1 + 1
        i3 = i2 + 1
        i4 = i3 + 1
        i5 = i4 + 1
        i6 = i5 + 1
        i7 = i6 + 1
        i8 = i7 + 1
        z1, z2, z3, z4, z5, z6, z7, z8 = mul8_cyclic_karatsuba(
            (x[i1], x[i2], x[i3], x[i4], x[i5], x[i6], x[i7], x[i8]),
            (y[i1], y[i2], y[i3], y[i4], y[i5], y[i6], y[i7], y[i8]))
        output[i1] = z1
        output[i2] = z2
        output[i3] = z3
        output[i4] = z4
        output[i5] = z5
        output[i6] = z6
        output[i7] = z7
        output[i8] = z8
    end
end
