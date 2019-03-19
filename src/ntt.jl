using Primes: factor


"""
Returns the inverse of `x` modulo `rr_modulus(T)`.
The modulus must be a prime number.
"""
function ff_inverse(x::T) where T <: AbstractRRElem
    m = rr_modulus_simple(T)
    x^(m-2)
    # TODO: can be implemented with `invmod()`, which is generally faster
    # for non-Montgomery numbers, and works for any moduli, not just for prime ones.
    # Currently blocked by Julia issue #29971.
    #val = convert(encompassing_type(T), x)
    #T(invmod(val, m))
end


"""
Finds a generator element for a finite field defined by type `T`
(that is, we assume that the modulus in `T` is prime).
This means that every power of the returned `g` from `1` to `M-1` produces
all the elements of the field (integers from `1` to `M-1`), and `g^(M-1) = 1`.
"""
function get_generator(::Type{T}) where T <: AbstractRRElem
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


"""
Returns the root of one for NTT
(an analogue of the root of one in FFT, `e^(-2pi * im / N)`;
the returned value also has the property `w^N = 1`).
`(M-1)` must be divisible by `N`, where `M` is the modulus of the type `T`.
"""
function get_root_of_one(::Type{T}, N::Integer, inverse::Bool) where T <: AbstractRRElem
    m = rr_modulus_simple(T)
    if mod(m - 1, N) != 0
        error("(modulus - 1) must be divisible by the NTT length")
    end
    g = get_generator(T) # g^(modulus(T) - 1) = 1
    w = g^div(m - 1, N)
    if inverse
        ff_inverse(w)
    else
        w
    end
end


"""
Returns the scaling coefficient for the inverse NTT.
Similarly to FFT, it's `1/N`, but in our case we need to take the finite field inverse.
"""
function get_inverse_coeff(::Type{T}, N::Integer) where T <: AbstractRRElem
    # Can also be calculated as `N^(-1) mod M == (M - (M-1) ÷ N)`
    ff_inverse(T(N))
end


"""
Reverses the order of the lowest `l` bits in `x` and fills the rest with 0s.
"""
function bitreverse(x::T, l::Integer) where T <: Integer
    # Slow, but simple function.
    # TODO: can be optimized if its performance becomes critical.
    parse(T, "0b" * reverse(bitstring(x)[end-l+1:end]))
end


"""
Returns the list of paris of indices for the elements that must be swapped
at the start of the FFT algorithm.
"""
function prepare_swap_indices(len::Integer)
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


"""
Retruns the array of arrays of twiddle factors for each stage of the FFT algorithm.
"""
function prepare_twiddle_factors(tp::Type, len::Integer, inverse::Bool)
    log_len = round(Int, log2(len))
    w = get_root_of_one(tp, len, inverse)
    twiddles = Array{tp, 1}[]
    for stage in 1:log_len
        mmax = 2^(stage-1)
        push!(twiddles, w.^((0:mmax-1) * 2^(log_len - stage)))
    end
    twiddles
end


"""
Prepared NTT plan.
"""
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
            w = get_root_of_one(tp, 2 * len, false)
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
