using Primes: isprime


"""
Select the best available polynomial multiplication function
based on the element type an polynomial length.
"""
@inline @generated function _get_polynomial_mul_function(
        ::Type{T}, ::Val{N}, ::Val{NC}) where {T, N, NC}
    if T <: AbstractRRElem
        m = rr_modulus_simple(T)
        # Regular NTT needs (m - 1) to be a multiple of N,
        # tangent NTT (for negacyclic polynomials) needs it to be a multiple of 2N.
        factor = NC ? 2 * N : N
        if rem(m - 1, factor) == 0 && isprime(m)
            :( ntt_mul )
        else
            :( karatsuba_mul )
        end
    else
        :( karatsuba_mul )
    end
end


"""
Polynomials modulo `x^n-1` (cyclic) or `x^n+1` (negacyclic).
Supports any type that has arithmetic operators defined for it,
including [`RRElem`](@ref) and [`RRElemMontgomery`](@ref).

    Polynomial(coeffs::AbstractArray{T, 1}, negacyclic::Bool) where T

Create a polynomial given the array of coefficients
(the `i`-th coefficient corresponds to the `(i-1)`-th power).
"""
struct Polynomial{T}
    coeffs :: Array{T, 1}
    negacyclic :: Bool
    mul_function :: Function

    @inline function Polynomial(coeffs::Array{T, 1}, negacyclic::Bool) where T
        len = length(coeffs)
        mul_function = _get_polynomial_mul_function(T, Val(len), Val(negacyclic))
        new{T}(coeffs, negacyclic, mul_function)
    end

    @inline function Polynomial{T}(coeffs::Array{V, 1}, negacyclic::Bool) where {T, V}
        Polynomial(convert.(T, coeffs), negacyclic)
    end

    @inline function Polynomial(coeffs::Array{T, 1}, negacyclic::Bool, mul_function) where T
        new{T}(coeffs, negacyclic, mul_function)
    end
end


# Required for broadcasting


Base.length(x::Polynomial{T}) where T = 1


Base.iterate(x::Polynomial{T}) where T = (x, nothing)
Base.iterate(x::Polynomial{T}, state) where T = nothing


# It has to match a polynomial with any size and modulus,
# So it can't be a `Polynomial` object.
struct ZeroPolynomial
end


@inline Base.zero(::Polynomial{T}) where T = ZeroPolynomial()


@inline function Base.:(==)(p1::Polynomial{T}, p2::Polynomial{T}) where T
    p1.negacyclic == p2.negacyclic && p1.coeffs == p2.coeffs
end


@inline function Base.:*(p1::Polynomial{T}, p2::Polynomial{T}) where T
    @assert p1.negacyclic == p2.negacyclic && length(p1.coeffs) == length(p2.coeffs)
    p1.mul_function(p1, p2)
end

@inline Base.:*(p1::Polynomial, p2::ZeroPolynomial) = ZeroPolynomial()
@inline Base.:*(p1::ZeroPolynomial, p2::Polynomial) = ZeroPolynomial()


@inline Base.:+(p1::Polynomial, p2::ZeroPolynomial) = p1
@inline Base.:+(p1::ZeroPolynomial, p2::Polynomial) = p2


@inline function Base.:*(p1::Polynomial{T}, p2::Integer) where T
    Polynomial(p1.coeffs .* convert(T, p2), p1.negacyclic, p1.mul_function)
end

@inline function Base.:*(p1::Polynomial{T}, p2::V) where T where V
    Polynomial(p1.coeffs .* convert(T, p2), p1.negacyclic, p1.mul_function)
end


@inline function Base.:*(p1::Integer, p2::Polynomial)
    p2 * p1
end

@inline function Base.:*(p1::Polynomial{T}, p2::T) where T
    Polynomial(p1.coeffs .* p2, p1.negacyclic, p1.mul_function)
end


@inline function Base.div(p1::Polynomial{T}, p2::Integer) where T
   Polynomial(div.(p1.coeffs, p2), p1.negacyclic, p1.mul_function)
end


@inline function Base.:+(p1::Polynomial{T}, p2::Polynomial{T}) where T
    Polynomial(p1.coeffs .+ p2.coeffs, p1.negacyclic, p1.mul_function)
end

@inline function Base.:+(p1::Polynomial{T}, p2::T) where T
    coeffs = copy(p1.coeffs)
    coeffs[1] += p2
    Polynomial(coeffs, p1.negacyclic, p1.mul_function)
end


@inline function Base.:-(p1::Polynomial{T}, p2::Polynomial{T}) where T
    Polynomial(p1.coeffs .- p2.coeffs, p1.negacyclic, p1.mul_function)
end

@inline function Base.:-(p1::Polynomial{T}, p2::Unsigned) where T
    Polynomial(p1.coeffs .- T(p2), p1.negacyclic, p1.mul_function)
end

@inline function Base.:-(p1::Polynomial{T}) where T
    Polynomial(.-p1.coeffs, p1.negacyclic, p1.mul_function)
end

@inline Base.:-(p1::Polynomial, p2::ZeroPolynomial) = p1


"""
    shift_polynomial(p::Polynomial{T}, shift::Integer) where T

Multiply the polynomial by `x^shift`. `shift` can be negative.
"""
@Base.propagate_inbounds @inline function shift_polynomial(
        p::Polynomial{T}, shift::Integer) where T

    if shift == 0
        p
    else
        n = length(p.coeffs)
        cycle = isodd(fld(shift, n))
        shift = mod(shift, n)

        shift_first = p.negacyclic && (!cycle)
        shift_last = p.negacyclic && cycle

        new_coeffs = similar(p.coeffs)
        coeffs = p.coeffs
        for j in 1:shift
            new_coeffs[j] = shift_first ? -coeffs[n-shift+j] : coeffs[n-shift+j]
        end
        for j in shift+1:n
            new_coeffs[j] = shift_last ? -coeffs[j-shift] : coeffs[j-shift]
        end

        Polynomial(new_coeffs, p.negacyclic, p.mul_function)
    end
end


"""
Multiplies two polynomials of length `l` whose coefficients are stored in arrays
`p1c` and `p2c` starting from locations `p1_s` and `p2_s`, respectively.
The result is placed in the array `res` starting from the location `res_s`.
`res` is assumed to have enough space to store the full result of length `2l`.
"""
@inline function mul_naive(
        l::Int, res::Array{T, 1}, res_s::Int,
        p1c::Array{T, 1}, p1_s::Int, p2c::Array{T, 1}, p2_s::Int) where T
    @simd for j in 1:l
        for k in 1:l
            res[res_s+j+k-2] += p2c[p2_s+k-1] * p1c[p1_s+j-1]
        end
    end
end


"""
Recursive Karatsuba multiplication function.
Multiplies two polynomials of length `full_len` (must be a power of 2)
whose coefficients are stored in arrays
`p1c` and `p2c` starting from locations `p1_s` and `p2_s`, respectively.
The result is placed in the array `res` starting from the location `res_s`.
`res` is assumed to have enough space to store the full result of length `2 * full_len`.
`buf` and `buf2` are temporary buffers with enough space for `2 * full_len` elements
starting from locations `buf_s` and `buf2_s`, respectively.
"""
@inline function _karatsuba_mul(
        full_len::Int, res::Array{T, 1}, res_s::Int,
        p1c::Array{T, 1}, p1_s::Int,
        p2c::Array{T, 1}, p2_s::Int,
        buf::Array{T, 1}, buf_s::Int,
        buf2::Array{T, 1}, buf2_s::Int) where T

    if full_len <= 8
        mul_naive(full_len, res, res_s, p1c, p1_s, p2c, p2_s)
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


"""
Multiplies two polynomials using Karatsuba algorithm.
Assumes the polynomials have the same length and the same value of the `negacyclic` field.
"""
@Base.propagate_inbounds @inline function karatsuba_mul(
        p1::Polynomial{T}, p2::Polynomial{T}) where T

    full_len = length(p1.coeffs)
    half_len = div(length(p1.coeffs), 2)

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

    if p1.negacyclic
        @simd for i in 1:half_len
            r0[i+half_len] += r1[i]
            r0[i] -= r1[i+half_len]
        end
        @simd for i in 1:full_len
            r0[i] -= r2[i]
        end
    else
        @simd for i in 1:half_len
            r0[i+half_len] += r1[i]
            r0[i] += r1[i+half_len]
        end
        @simd for i in 1:full_len
            r0[i] += r2[i]
        end
    end

    Polynomial(r0, p1.negacyclic, p1.mul_function)
end


"""
Multiplies two polynomials using NTT.
Assumes the polynomials have the same length and the same value of the `negacyclic` field.
"""
@inline function ntt_mul(p1::Polynomial{T}, p2::Polynomial{T}) where T
    plan = get_ntt_plan(T, length(p1.coeffs), p1.negacyclic)
    c1 = similar(p1.coeffs)
    ntt!(plan, c1, p1.coeffs)
    c2 = similar(p2.coeffs)
    ntt!(plan, c2, p2.coeffs)
    c1 .*= c2
    intt!(plan, c2, c1)
    Polynomial(c2, p1.negacyclic, p1.mul_function)
end


"""
Multiplies two polynomials using Nussbaumer convolution algorithm.
Assumes the polynomials have the same length and the same value of the `negacyclic` field.
"""
@Base.propagate_inbounds @inline function nussbaumer_mul(
        p1::Polynomial{T}, p2::Polynomial{T}) where T
    if p1.negacyclic
        res = nussbaumer_mul_negacyclic(p1.coeffs, p2.coeffs, true)
    else
        res = nussbaumer_mul_cyclic(p1.coeffs, p2.coeffs, true)
    end
    Polynomial(res, p1.negacyclic, p1.mul_function)
end
