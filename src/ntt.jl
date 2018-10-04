using Primes: factor


rr_modulus(::Type{RRElem{T, M}}) where {T, M} = M
rr_modulus(::Type{RRElemMontgomery{T, M}}) where {T, M} = M


function rr_modulus_simple(tp::Type{<:AbstractRRElem})
    m = rr_modulus(tp)
    convert(encompassing_type(typeof(m)), m)
end


function ff_inverse(x::T) where T <: AbstractRRElem
    m = rr_modulus_simple(T)
    x^(m - 2)
end


function get_root(::Type{T}) where T <: AbstractRRElem
    modulus = rr_modulus_simple(T)
    factors = keys(factor(modulus - 1))
    for w in 2:modulus-1
        found = true
        gw = T(w)
        for q in factors
            if gw^(div(modulus - 1, q)) == one(T)
                found = false
                break
            end
        end
        if found
            return gw
        end
    end
    zero(T)
end


function get_twiddle_base(::Type{T}, N, inverse) where T <: AbstractRRElem
    m = rr_modulus_simple(T)
    w = get_root(T) # w^(modulus(T) - 1) = 1
    b = w^div(m - 1, N)
    if inverse
        ff_inverse(b)
    else
        b
    end
end


function get_inverse_coeff(::Type{T}, N) where T <: AbstractRRElem
    # To compare with the DFT
    m = rr_modulus_simple(T)
    T(m - div(m - 1, N))
end


function bitreverse(x, l)
    # Slow, but simple function. Needs some optimized algorithm in practice.
    parse(Int, "0b" * reverse(bitstring(x)[end-l+1:end]))
end


function prepare_swap_indices(len)
    log_len = round(Int, log2(len))
    indices = Tuple{Int, Int}[]
    for i in 1:len
        j = bitreverse(i-1, log_len) + 1
        if j > i
            push!(indices, (i, j))
        end
    end
    indices
end


function prepare_twiddle_factors(tp, len, inverse)
    log_len = round(Int, log2(len))
    w = get_twiddle_base(tp, len, inverse)
    twiddles = Array{tp, 1}[]
    for stage in 1:log_len
        mmax = 2^(stage-1)
        push!(twiddles, w.^((0:mmax-1) * 2^(log_len - stage)))
    end
    twiddles
end


struct NTTPlan{T <: AbstractRRElem}
    forward_coeffs :: Array{T, 1}
    use_forward_coeffs :: Bool
    inverse_coeffs :: Array{T, 1}
    use_inverse_coeffs :: Bool
    forward_twiddle_coeffs :: Array{Array{T, 1}}
    inverse_twiddle_coeffs :: Array{Array{T, 1}}
    swap_indices :: Array{Tuple{Int, Int}, 1}

    function NTTPlan(tp::Type{<: AbstractRRElem}, len::Int, tangent::Bool)

        if len < 2 || len & (len - 1) != 0
            error("len must be >=2 and a power of 2")
        end

        if tangent
            w = get_twiddle_base(tp, 2 * len, false)
            idx = collect(0:len-1)

            forward_coeffs = w.^idx
            use_forward_coeffs = true
            inverse_coeffs = get_inverse_coeff(tp, len) .* (w.^(mod.(2 * len .- idx, 2 * len)))
        else
            forward_coeffs = ones(tp, len)
            use_forward_coeffs = false
            inverse_coeffs = ones(tp, len) * get_inverse_coeff(tp, len)
        end

        new{tp}(
            forward_coeffs,
            use_forward_coeffs,
            inverse_coeffs,
            true,
            prepare_twiddle_factors(tp, len, false),
            prepare_twiddle_factors(tp, len, true),
            prepare_swap_indices(len))
    end

end


const _ntt_plans = Dict{Tuple{Type, Int, Bool}, NTTPlan}()


function get_ntt_plan(tp::Type{<: AbstractRRElem}, len::Int, tangent::Bool)
    key = (tp, len, tangent)
    if !haskey(_ntt_plans, key)
        plan = NTTPlan(tp, len, tangent)
        _ntt_plans[key] = plan
        plan
    else
        _ntt_plans[key]
    end
end


@Base.propagate_inbounds function ntt!(plan::NTTPlan{T}, data::Array{T, 1}, inverse::Bool) where T
    len = length(data)

    if inverse
        twiddle_coeffs = plan.inverse_twiddle_coeffs
    else
        twiddle_coeffs = plan.forward_twiddle_coeffs
    end

    if !inverse && plan.use_forward_coeffs
        data .*= plan.forward_coeffs
    end

    for (i, j) in plan.swap_indices
        data[j], data[i] = data[i], data[j]
    end

    logn = round(Int, log2(len))
    for stage in 1:logn
        mmax = 2^(stage-1)
        istep = mmax * 2
        ws = twiddle_coeffs[stage]
        for m = 1:mmax
            w = ws[m]
            for i = m:istep:len
                j = i + mmax
                temp = w * data[j]
                data[j] = data[i] - temp
                data[i] += temp
            end
        end
    end

    if inverse && plan.use_inverse_coeffs
        data .*= plan.inverse_coeffs
    end
end


function ntt(data::Array{T, 1}, inverse::Bool) where T <: AbstractRRElem
    plan = get_ntt_plan(T, length(data), false)
    data = copy(data)
    ntt!(plan, data, inverse)
    data
end
