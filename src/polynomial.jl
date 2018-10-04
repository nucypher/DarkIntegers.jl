"""
Polynomials modulo `x^n+1`.
The element type can be a modulo number (as long as it has arithmetic operators defined for it).
"""
struct Polynomial{T}
    coeffs :: Array{T, 1}
    negacyclic :: Bool

    function Polynomial(
            ::Type{T}, coeffs::AbstractArray{V, 1}, negacyclic) where T where V <: Integer
        coeffs_rm = T.(coeffs)
        new{T}(coeffs_rm, negacyclic)
    end

    function Polynomial(coeffs::Array{T, 1}, negacyclic) where T
        new{T}(coeffs, negacyclic)
    end
end


# It has to match a polynomial with any size and modulus,
# So it can't be a `Polynomial` object.
struct ZeroPolynomial
end


Base.zero(::Polynomial{T}) where T = ZeroPolynomial()

Base.length(p::Polynomial{T}) where T = length(p.coeffs)


function Base.:(==)(p1::Polynomial{T}, p2::Polynomial{T}) where T
    p1.negacyclic == p2.negacyclic && p1.coeffs == p2.coeffs
end


function Base.:*(p1::Polynomial{T}, p2::Polynomial{T}) where T
    karatsuba_mul(p1, p2)
    #fast_reference_mul(p1.coeffs, p2.coeffs)
end

Base.:*(p1::Polynomial, p2::ZeroPolynomial) = ZeroPolynomial()
Base.:*(p1::ZeroPolynomial, p2::Polynomial) = ZeroPolynomial()


function Base.:*(p1::Polynomial{T}, p2::Integer) where T
    Polynomial(p1.coeffs .* convert(T, p2), p1.negacyclic)
end

function Base.:*(p1::Polynomial{T}, p2::V) where T where V
    Polynomial(p1.coeffs .* convert(T, p2), p1.negacyclic)
end


function Base.:*(p1::Integer, p2::Polynomial)
    p2 * p1
end

function Base.:*(p1::Polynomial{T}, p2::T) where T
    Polynomial(p1.coeffs .* p2, p1.negacyclic)
end


function Base.:+(p1::Polynomial{T}, p2::Polynomial{T}) where T
    Polynomial(p1.coeffs .+ p2.coeffs, p1.negacyclic)
end

function Base.:+(p1::Polynomial{T}, p2::T) where T
    coeffs = copy(p1.coeffs)
    coeffs[1] += p2
    Polynomial(coeffs, p1.negacyclic)
end


function Base.:-(p1::Polynomial{T}, p2::Polynomial{T}) where T
    Polynomial(p1.coeffs .- p2.coeffs, p1.negacyclic)
end


function Base.:-(p1::Polynomial{T}, p2::Unsigned) where T
    Polynomial(p1.coeffs .- T(p2), p1.negacyclic)
end


Base.:-(p1::Polynomial, p2::ZeroPolynomial) = p1


with_modulus(p::Polynomial{T}, new_modulus::V) where T where V =
    # TODO: technically, we need to only convert the modulus from Integer once
    Polynomial(with_modulus.(p.coeffs, new_modulus), p.negacyclic)


function with_length(p::Polynomial{T}, new_length::Integer) where T
    @assert new_length >= length(p)
    Polynomial([p.coeffs; zeros(eltype(p.coeffs), new_length - length(p))], p.negacyclic)
end



function modulus_reduction(p::Polynomial{T}, new_modulus::Unsigned) where T
    # TODO: technically, we need to only convert the modulus from Integer once
    Polynomial(modulus_reduction.(p.coeffs, new_modulus), p.negacyclic)
end


@Base.propagate_inbounds function shift_polynomial(p::Polynomial{T}, shift::Integer) where T
    if shift == 0
        p
    else
        n = length(p)
        cycle = mod(fld(shift, n), 2)
        shift = mod(shift, n)

        global_neg = (cycle == 1 && p.negacyclic)
        shift_neg = p.negacyclic

        new_coeffs = similar(p.coeffs)
        coeffs = p.coeffs
        for j in 1:shift
            new_coeffs[j] = !(global_neg && shift_neg) ? -coeffs[n-shift+j] : coeffs[n-shift+j]
        end
        for j in shift+1:n
            new_coeffs[j] = global_neg ? -coeffs[j-shift] : coeffs[j-shift]
        end

        Polynomial(new_coeffs, p.negacyclic)
    end
end



function reference_mul(p1::Polynomial{T}, p2::Polynomial{T}) where T
    res = Polynomial(zeros(T, length(p1)), p1.negacyclic)
    for (j, c) in enumerate(p1.coeffs)
        res = res + shift_polynomial(p2, j - 1) * c
    end
    res
end


@Base.propagate_inbounds @inline function fast_reference_mul(
        p1::Polynomial{T}, p2::Polynomial{T}) where T

    res = zeros(T, length(p1))

    for j in 1:length(p1)
        c = p1.coeffs[j]

        cc = p1.negacyclic ? -c : c
        for k in 1:j-1
            res[k] += p2.coeffs[end-j+1+k] * cc
        end
        for k in j:length(p1)
            res[k] += p2.coeffs[k-j+1] * c
        end

    end
    Polynomial(res, p1.negacyclic)
end


@inline function mul_with_overflow(l, res, res_s, p1, p1_s, p2, p2_s)
    @simd for j in 1:l
        for k in 1:l
            res[res_s+j+k-2] += p2[p2_s+k-1] * p1[p1_s+j-1]
        end
    end
end


@inline function _karatsuba_mul(full_len, res, res_s, p1c, p1_s, p2c, p2_s, buf, buf_s,
        buf2, buf2_s)

    if full_len <= 8
        mul_with_overflow(full_len, res, res_s, p1c, p1_s, p2c, p2_s)
        return
    end

    half_len = div(full_len, 2)

    _karatsuba_mul(
        half_len, res, res_s, p1c, p1_s, p2c, p2_s, buf, buf_s+full_len, buf2, buf2_s)
    _karatsuba_mul(
        half_len, res, res_s+full_len, p1c, p1_s+half_len, p2c, p2_s+half_len, buf, buf_s+full_len,
        buf2, buf2_s)

    z = zero(eltype(p1c))
    @simd for i in buf_s:buf_s+full_len-1
        buf[i] = z
    end

    @simd for i in 1:half_len
        buf2[buf2_s+i-1] = p1c[p1_s+i-1] + p1c[p1_s+half_len+i-1]
        buf2[buf2_s+half_len+i-1] = p2c[p2_s+i-1] + p2c[p2_s+half_len+i-1]
    end

    _karatsuba_mul(half_len, buf, buf_s, buf2, buf2_s, buf2, buf2_s+half_len, buf, buf_s+full_len,
        buf2, buf2_s+full_len)

    @simd for i in 1:full_len
        buf2[buf2_s+i-1] = buf[buf_s+i-1] - res[res_s+i-1] - res[res_s+i+full_len-1]
    end

    @simd for i in 1:full_len
        res[res_s+i+half_len-1] += buf2[buf2_s+i-1]
    end
end


@Base.propagate_inbounds @inline function karatsuba_mul(
        p1::Polynomial{T}, p2::Polynomial{T}) where T

    full_len = length(p1)
    half_len = div(length(p1), 2)

    z = zero(T)
    r0 = similar(p1.coeffs)
    r1 = similar(p1.coeffs)
    r2 = similar(p1.coeffs)
    buf = similar(p1.coeffs)
    buf2 = similar(p1.coeffs)
    r3 = similar(p1.coeffs)

    @simd for i in 1:full_len
        r0[i] = z
        r1[i] = z
        r2[i] = z
    end

    _karatsuba_mul(half_len, r0, 1, p1.coeffs, 1, p2.coeffs, 1, buf, 1, buf2, 1)
    _karatsuba_mul(half_len, r2, 1, p1.coeffs, half_len+1, p2.coeffs, half_len+1, buf, 1, buf2, 1)

    @simd for i in 1:half_len
        r3[i] = p1.coeffs[i] + p1.coeffs[i+half_len]
        r3[i+half_len] = p2.coeffs[i] + p2.coeffs[i+half_len]
    end
    _karatsuba_mul(half_len, r1, 1, r3, 1, r3, half_len+1, buf, 1, buf2, 1)
    r1 .-= r2 .+ r0

    # TODO: assuming p.negacyclic here, adding a variant for -1 should be simple
    @simd for i in 1:half_len
        r0[i+half_len] += r1[i]
        r0[i] -= r1[i+half_len]
    end
    @simd for i in 1:full_len
        r0[i] -= r2[i]
    end

    Polynomial(r0, p1.negacyclic)
end


# Multiplication of negacyclic polynomials based on tangent NTT
# TODO: add support for posicyclic polynomials, using regular NTT
function ntt_mul(p1::Polynomial{T}, p2::Polynomial{T}) where T
    @assert p1.negacyclic && p2.negacyclic
    plan = get_ntt_plan(T, length(p1), true)
    c1 = copy(p1.coeffs)
    ntt!(plan, c1, false)
    c2 = copy(p2.coeffs)
    ntt!(plan, c2, false)
    c1 .*= c2
    ntt!(plan, c1, true)
    Polynomial(c1, true)
end
