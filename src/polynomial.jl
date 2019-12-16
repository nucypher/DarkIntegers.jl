abstract type AbstractPolynomialModulus end


abstract type AbstractCyclicModulus <: AbstractPolynomialModulus end


# Polynomials modulo x^N+1, where N-1 is the polynomial degree
struct NegacyclicModulus <: AbstractCyclicModulus end


const negacyclic_modulus = NegacyclicModulus()


# Polynomials modulo x^N-1, where N-1 is the polynomial degree
struct CyclicModulus <: AbstractCyclicModulus end


const cyclic_modulus = CyclicModulus()


"""
Select the best available polynomial multiplication function
based on the element type an polynomial length.
"""
@inline function _get_polynomial_mul_function(
        ::Type{T}, ::Val{N}, ::AbstractPolynomialModulus) where {T, N}
    karatsuba_mul
end

@inline @generated function _get_polynomial_mul_function(
        ::Type{T}, ::Val{N}, pm::AbstractCyclicModulus) where {T <: AbstractModUInt, N}
    m = modulus_as_builtin(T)
    # Regular NTT needs (m - 1) to be a multiple of N,
    # tangent NTT (for negacyclic polynomials) needs it to be a multiple of 2N.
    factor = pm == NegacyclicModulus ? 2 * N : N
    if rem(m - 1, factor) == 0 && isprime(m)
        :( ntt_mul )
    else
        :( karatsuba_mul )
    end
end


"""
Fixed-size polynomials (with the degree limited by `N-1`).
Currently supports moduli `x^n-1` (`cyclic_modulus`) or `x^n+1` (`negacyclic_modulus`).
Supports any type that has arithmetic operators defined for it,
including [`ModUInt`](@ref) and [`MgModUInt`](@ref).

    Polynomial(coeffs::AbstractArray{T, 1}, modulus::AbstractCyclicModulus) where T

Create a polynomial given the array of coefficients
(the `i`-th coefficient corresponds to the `(i-1)`-th power).
"""
struct Polynomial{T, N} <: AbstractArray{T, 1}
    coeffs :: Array{T, 1}
    modulus :: AbstractPolynomialModulus
    mul_function :: Function

    @inline function Polynomial(
            coeffs::Array{T, 1}, modulus::AbstractCyclicModulus, mul_function) where T
        new{T, length(coeffs)}(coeffs, modulus, mul_function)
    end

    @inline function Polynomial(coeffs::Array{T, 1}, modulus::AbstractCyclicModulus) where T
        len = length(coeffs)
        mul_function = _get_polynomial_mul_function(T, Val(len), modulus)
        Polynomial(coeffs, modulus, mul_function)
    end
end


# It has to match a polynomial with any modulus,
# So it can't be a `Polynomial` object.
struct ZeroPolynomial{T, N}
end


@inline Base.zero(::Polynomial{T, N}) where {T, N} = ZeroPolynomial{T, N}()


@inline function Base.:(==)(p1::Polynomial{T, N}, p2::Polynomial{T, N}) where {T, N}
    p1.modulus == p2.modulus && p1.coeffs == p2.coeffs
end


@inline function Base.:*(p1::Polynomial{T, N}, p2::Polynomial{T, N}) where {T, N}
    @assert p1.modulus == p2.modulus && length(p1.coeffs) == length(p2.coeffs)
    p1.mul_function(p1, p2)
end

@inline Base.:*(p1::Polynomial{T, N}, p2::ZeroPolynomial{T, N}) where {T, N} =
    ZeroPolynomial{T, N}()

@inline Base.:*(p1::ZeroPolynomial{T, N}, p2::Polynomial{T, N}) where {T, N} =
    ZeroPolynomial{T, N}()

@inline function Base.:*(p1::Polynomial{T, N}, p2::Integer) where {T, N}
    Polynomial(p1.coeffs .* convert(T, p2), p1.modulus, p1.mul_function)
end

@inline function Base.:*(p1::Polynomial{T, N}, p2::V) where {T, N, V}
    Polynomial(p1.coeffs .* convert(T, p2), p1.modulus, p1.mul_function)
end

@inline function Base.:*(p1::Integer, p2::Polynomial)
    p2 * p1
end

@inline function Base.:*(p1::Polynomial{T, N}, p2::T) where {T, N}
    Polynomial(p1.coeffs .* p2, p1.modulus, p1.mul_function)
end


@inline function Base.div(p1::Polynomial{T, N}, p2::Integer) where {T, N}
   Polynomial(div.(p1.coeffs, p2), p1.modulus, p1.mul_function)
end


@inline function Base.:+(p1::Polynomial{T, N}, p2::Polynomial{T}) where {T, N}
    Polynomial(p1.coeffs .+ p2.coeffs, p1.modulus, p1.mul_function)
end

@inline function Base.:+(p1::Polynomial{T, N}, p2::T) where {T, N}
    coeffs = copy(p1.coeffs)
    coeffs[1] += p2
    Polynomial(coeffs, p1.modulus, p1.mul_function)
end

@inline Base.:+(p1::Polynomial{T, N}, p2::ZeroPolynomial{T, N}) where {T, N} = p1

@inline Base.:+(p1::ZeroPolynomial{T, N}, p2::Polynomial{T, N}) where {T, N} = p2


@inline function Base.:-(p1::Polynomial{T, N}, p2::Polynomial{T, N}) where {T, N}
    Polynomial(p1.coeffs .- p2.coeffs, p1.modulus, p1.mul_function)
end

@inline function Base.:-(p1::Polynomial{T, N}, p2::Unsigned) where {T, N}
    Polynomial(p1.coeffs .- convert(T, p2), p1.modulus, p1.mul_function)
end

@inline function Base.:-(p1::Polynomial{T, N}) where {T, N}
    Polynomial(.-p1.coeffs, p1.modulus, p1.mul_function)
end

@inline Base.:-(p1::Polynomial{T, N}, p2::ZeroPolynomial{T, N}) where {T, N} = p1

@inline Base.:-(p1::ZeroPolynomial{T, N}, p2::Polynomial{T, N}) where {T, N} = -p2


"""
    mul_by_monomial(p::Polynomial, power::Integer)

Multiply the polynomial by `x^power`.
If `power` lies outside `[0, 2 * N)` where `N-1` is the maximum degree of the polynomial,
a modulo `2 * N` will be taken.
"""
@Base.propagate_inbounds @inline function mul_by_monomial(p::Polynomial, power::Integer)

    @assert isa(p.modulus, AbstractCyclicModulus)

    if power == 0
        p
    else
        n = length(p.coeffs)
        cycle = isodd(fld(power, n))
        power = mod(power, n)

        shift_first = p.modulus == negacyclic_modulus && (!cycle)
        shift_last = p.modulus == negacyclic_modulus && cycle

        new_coeffs = similar(p.coeffs)
        coeffs = p.coeffs
        for j in 1:power
            new_coeffs[j] = shift_first ? -coeffs[n-power+j] : coeffs[n-power+j]
        end
        for j in power+1:n
            new_coeffs[j] = shift_last ? -coeffs[j-power] : coeffs[j-power]
        end

        Polynomial(new_coeffs, p.modulus, p.mul_function)
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
Assumes the polynomials have the same length and the same value of the `modulus` field.
"""
@Base.propagate_inbounds @inline function karatsuba_mul(
        p1::Polynomial{T, N}, p2::Polynomial{T, N}) where {T, N}

    @assert p1.modulus == p2.modulus
    @assert isa(p1.modulus, AbstractCyclicModulus)

    full_len = N
    half_len = N ÷ 2

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

    if p1.modulus == negacyclic_modulus
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

    Polynomial(r0, p1.modulus, p1.mul_function)
end


"""
Multiplies two polynomials using NTT.
Assumes the polynomials have the same length and the same value of the `modulus` field.
"""
@inline function ntt_mul(p1::Polynomial{T, N}, p2::Polynomial{T, N}) where {T, N}

    @assert p1.modulus == p2.modulus
    @assert isa(p1.modulus, AbstractCyclicModulus)

    plan = get_ntt_plan(T, N, p1.modulus == negacyclic_modulus)
    c1 = similar(p1.coeffs)
    ntt!(plan, c1, p1.coeffs)
    c2 = similar(p2.coeffs)
    ntt!(plan, c2, p2.coeffs)
    c1 .*= c2
    intt!(plan, c2, c1)
    Polynomial(c2, p1.modulus, p1.mul_function)
end


"""
Multiplies two polynomials using Nussbaumer convolution algorithm.
Assumes the polynomials have the same length and the same value of the `modulus` field.
"""
@Base.propagate_inbounds @inline function nussbaumer_mul(
        p1::Polynomial{T, N}, p2::Polynomial{T, N}) where {T, N}

    @assert p1.modulus == p2.modulus
    @assert isa(p1.modulus, AbstractCyclicModulus)

    if p1.modulus == negacyclic_modulus
        res = nussbaumer_mul_negacyclic(p1.coeffs, p2.coeffs, true)
    else
        res = nussbaumer_mul_cyclic(p1.coeffs, p2.coeffs, true)
    end
    Polynomial(res, p1.modulus, p1.mul_function)
end


# Broadcasting machinery
# Is this really the simplest way to do it?

Base.length(x::Polynomial) = length(x.coeffs)


Base.size(x::Polynomial) = size(x.coeffs)


Base.getindex(x::Polynomial, inds::Vararg{Int, N}) where N = x.coeffs[inds...]


Base.setindex!(x::Polynomial, val, inds::Vararg{Int, N}) where N = x.coeffs[inds...] = val


Base.BroadcastStyle(::Type{<:Polynomial}) = Broadcast.ArrayStyle{Polynomial}()


function Base.similar(
        bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{Polynomial}},
        ::Type{ElType}) where ElType
    x = _find_polynomial(bc)
    @assert length(axes(bc)) == 1
    Polynomial(similar(Array{ElType}, axes(bc)), x.modulus)
end


_find_polynomial(bc::Base.Broadcast.Broadcasted) = _find_polynomial(bc.args)
_find_polynomial(args::Tuple) = _find_polynomial(_find_polynomial(args[1]), Base.tail(args))
_find_polynomial(x) = x
_find_polynomial(x::Base.Broadcast.Extruded) = x.x
_find_polynomial(::Tuple{}) = nothing
_find_polynomial(a::Polynomial, rest) = a
_find_polynomial(::Any, rest) = _find_polynomial(rest)


Base.iterate(x::Polynomial) = iterate(x.coeffs)

Base.iterate(x::Polynomial, state) = iterate(x.coeffs, state)
