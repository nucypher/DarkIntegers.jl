using Primes: factor


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
        inv(w)
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
    inv(T(N))
end


"""
Reverses the order of the lowest `l` bits in `x` and fills the rest with 0s.
"""
function bitreverse(x::T, l::Integer) where T <: Integer
    # Slow, but simple function.
    # TODO: can be optimized if its performance becomes critical.
    if l == 0
        zero(T)
    else
        parse(T, "0b" * reverse(bitstring(x)[end-l+1:end]))
    end
end


"""
Prepared NTT plan.

If `tangent` is `false`, the NTT can be used to multiply polynomials modulo `X^len-1` (cyclic),
if it is `true`, it can be used to multiply polynomials modulo `X^len+1` (negacyclic).
"""
struct NTTPlan{T <: Unsigned, M}
    inverse_coeff :: RRElemMontgomery{T, M}
    forward_twiddle_coeffs :: Array{RRElemMontgomery{T, M}, 1}
    inverse_twiddle_coeffs :: Array{RRElemMontgomery{T, M}, 1}

    function NTTPlan(::Type{T}, modulus::T, len::Int, tangent::Bool) where T

        #=
        We store coefficients in Montgomery representation regardless of the actual array type.
        Since Montgomery multiplication for `x` and `y` gives `x * y * R^(-1) mod M`,
        it works both when `x` is an `RRElem` and `y` is an `RRElemMontgomery`,
        and when `x` and `y` are `RRElemMontgomery`.
        This makes multiplication by coefficients faster.
        =#
        coeffs_tp = RRElemMontgomery{T, modulus}

        if len < 2 || len & (len - 1) != 0
            error("len must be >=2 and a power of 2")
        end

        log_len = trailing_zeros(len)

        forward_twiddle_coeffs = Array{coeffs_tp}(undef, len)
        inverse_twiddle_coeffs = Array{coeffs_tp}(undef, len)

        if tangent
            w = get_root_of_one(coeffs_tp, 2 * len, false)
        else
            w = get_root_of_one(coeffs_tp, len, false)
        end
        w_inv = inv(w)

        forward_twiddle_coeffs[1] = 0 # unused
        inverse_twiddle_coeffs[1] = 0 # unused
        for stage in 0:log_len-1
            m = 1 << stage
            for i in 0:m-1
                if tangent
                    pwr = bitreverse(m + i, log_len)
                else
                    pwr = bitreverse(i, stage) << (log_len - stage - 1)
                end

                forward_twiddle_coeffs[m + i + 1] = w^pwr
                inverse_twiddle_coeffs[m + i + 1] = w_inv^pwr
            end
        end

        # Similarly to FFT, the scaling coefficient for the inverse transform is also `1/N`,
        # but in our case we need to take the finite field inverse.
        # Can also be calculated as `N^(-1) mod M == (M - (M-1) ÷ N)`.
        inverse_coeff = get_inverse_coeff(coeffs_tp, len)

        new{T, modulus}(
            inverse_coeff,
            forward_twiddle_coeffs,
            inverse_twiddle_coeffs)
    end

end


const _ntt_plans = Dict{Tuple{Type, Val, Int, Bool}, NTTPlan}()


function _get_ntt_plan(::Type{T}, modulus::T, len::Int, tangent::Bool) where T
    key = (T, Val(modulus), len, tangent)
    if !haskey(_ntt_plans, key)
        plan = NTTPlan(T, modulus, len, tangent)
        _ntt_plans[key] = plan
        plan
    else
        _ntt_plans[key]
    end
end


get_ntt_plan(::Type{RRElem{T, M}}, len::Int, tangent::Bool) where {T, M} =
    _get_ntt_plan(T, M, len, tangent)
get_ntt_plan(::Type{RRElemMontgomery{T, M}}, len::Int, tangent::Bool) where {T, M} =
    _get_ntt_plan(T, M, len, tangent)


#=
We are using a modification of the "decimation-in-frequency" FFT algorithm.
This layout seems to be optimal for automatic vectorization.

Since the primary purpose of NTT is polynomial multiplication,
the order of the elements in the result is not important, so we skip the rearrangement step.
=#


@Base.propagate_inbounds function ntt!(
        plan::NTTPlan{T, M}, output::Array{V, 1}, res::Array{V, 1}) where {T, M, V <: AbstractRRElem{T, M}}

    len = length(res)
    log_len = trailing_zeros(len)

    output .= res

    @inbounds for stage in 0:log_len-1
        m = 1 << stage
        logt1 = log_len - stage
        t = 1 << (logt1 - 1)
        @simd for i in 0:m-1
            j1 = i << logt1
            j2 = j1 + t - 1
            w = plan.forward_twiddle_coeffs[m + i + 1]
            for j in j1:j2
                a = output[j+1]
                temp = output[j+t+1] * w
                output[j+1] = a + temp
                output[j+t+1] = a - temp
            end
        end
    end
end


@Base.propagate_inbounds function intt!(
        plan::NTTPlan{T, M}, output::Array{V, 1}, res::Array{V, 1}) where {T, M, V <: AbstractRRElem{T, M}}

    len = length(res)
    log_len = trailing_zeros(len)

    output .= res

    @inbounds for stage in log_len-1:-1:0
        m = 1 << stage
        logt1 = log_len - stage
        t = 1 << (logt1 - 1)
        @simd for i in 0:m-1
            j1 = i << logt1
            j2 = j1 + t - 1
            w = plan.inverse_twiddle_coeffs[m + i + 1]
            for j in j1:j2
                a = output[j+1]
                b = output[j+t+1]
                output[j+1] = a + b
                output[j+t+1] = (a - b) * w
            end
        end
    end

    output .*= plan.inverse_coeff
end


function ntt(data::Array{T, 1}; inverse::Bool=false, negacyclic::Bool=false) where T <: AbstractRRElem
    plan = get_ntt_plan(T, length(data), negacyclic)
    output = similar(data)
    if inverse
        intt!(plan, output, data)
    else
        ntt!(plan, output, data)
    end
    output
end
