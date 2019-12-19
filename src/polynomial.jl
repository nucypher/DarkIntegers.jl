abstract type AbstractPolynomialModulus end


struct _UndefinedModulus <: AbstractPolynomialModulus end


const _undefined_modulus = _UndefinedModulus()


abstract type AbstractCyclicModulus <: AbstractPolynomialModulus end


# Polynomials modulo x^N+1, where N-1 is the polynomial degree
struct NegacyclicModulus <: AbstractCyclicModulus end


"""
A constant denoting negacyclic polynomial modulus (`x^N + 1`),
to be supplied to the [`Polynomial`](@ref) constructor.
"""
const negacyclic_modulus = NegacyclicModulus()


# Polynomials modulo x^N-1, where N-1 is the polynomial degree
struct CyclicModulus <: AbstractCyclicModulus end


"""
A constant denoting cyclic polynomial modulus (`x^N - 1`),
to be supplied to the [`Polynomial`](@ref) constructor.
"""
const cyclic_modulus = CyclicModulus()


"""
    known_isprime(::Val{X})

A method of this function can be defined by the user for a certain `X` to avoid
`isprime()` being called on it when determining the multiplication function for a polynomial
(which can reduce start-up time if `X` is very large).
"""
@generated function known_isprime(::Val{X}) where X
    res = isprime(as_builtin(X))
    :( $res )
end


"""
Select the best available polynomial multiplication function
based on the element type an polynomial length.
"""
@inline function _get_polynomial_mul_function(
        ::Type{T}, ::Val{N}, ::AbstractPolynomialModulus) where {T, N}
    nothing
end

@inline function _get_polynomial_mul_function(
        ::Type{T}, ::Val{N}, ::AbstractCyclicModulus) where {T, N}
    karatsuba_mul
end

# Ideally it should have been a @generated function,
# but in this case it does not pick up user-defined `known_isprime()` methods.
@inline function _get_polynomial_mul_function(
        ::Type{T}, ::Val{N}, pm::AbstractCyclicModulus) where {T <: AbstractModUInt, N}
    m = modulus(T)
    tp = eltype(T)
    # Regular NTT needs (m - 1) to be a multiple of N,
    # tangent NTT (for negacyclic polynomials) needs it to be a multiple of 2N.
    factor = pm == NegacyclicModulus() ? double(N) : N

    # Often `factor` will be a power of 2, in which case we can make the remainder check faster
    log2_factor = trailing_zeros(factor)
    m_minus_one = m - one(tp)
    if factor == 1 << log2_factor
        no_rem = iszero(m_minus_one & ((one(tp) << log2_factor) - one(tp)))
    else
        no_rem = iszero(rem(m_minus_one, factor))
    end

    if no_rem && known_isprime(Val(m))
        ntt_mul
    else
        karatsuba_mul
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
struct Polynomial{T, N} <: Number
    coeffs :: Array{T, 1}
    modulus :: AbstractPolynomialModulus
    mul_function :: Union{Function, Nothing}
    is_zero :: Bool

    @inline function Polynomial{T, N}(
            coeffs::Array{T, 1}, modulus::AbstractPolynomialModulus,
            mul_function, is_zero::Bool=false) where {T, N}
        @assert is_zero || modulus != _undefined_modulus
        @assert length(coeffs) == N
        new{T, N}(coeffs, modulus, mul_function, is_zero)
    end

    @inline function Polynomial{T, N}(
            coeffs::Array{T, 1}, modulus::AbstractPolynomialModulus,
            is_zero::Bool=false) where {T, N}
        len = length(coeffs)
        mul_function = _get_polynomial_mul_function(T, Val(len), modulus)
        Polynomial{T, N}(coeffs, modulus, mul_function, is_zero)
    end

    @inline function Polynomial(
            coeffs::Array{T, 1}, modulus::AbstractPolynomialModulus,
            is_zero::Bool=false) where T
        len = length(coeffs)
        Polynomial{T, len}(coeffs, modulus, is_zero)
    end
end


function Base.convert(::Type{Polynomial{T, N}}, x::Polynomial{V, M}) where {T, N, V, M}
    if x.is_zero
        zero(Polynomial{T, N})
    else
        @assert N >= M
        Polynomial{T, N}(convert.(T, x.coeffs), x.modulus)
    end
end


Base.promote_type(::Type{Polynomial{T, N}}, ::Type{Polynomial{V, M}}) where {T, N, V, M} =
    Polynomial{promote_type(T, V), max(N, M)}


@inline Base.zero(::Type{Polynomial{T, N}}) where {T, N} =
    Polynomial{T, N}(Array{T}(undef, N), _undefined_modulus, true)

Base.zero(::Polynomial{T, N}) where {T, N} = zero(Polynomial{T, N})


"""
    with_modulus(p::Polynomial{T, N}, new_modulus::AbstractPolynomialModulus)

Returns a new polynomial object with a changed modulus.
"""
function with_modulus(p::Polynomial{T, N}, m::Union{CyclicModulus, NegacyclicModulus}) where {T, N}
    if p.is_zero
        coeffs = zeros(T, N)
    else
        coeffs = p.coeffs
    end
    Polynomial{T, N}(coeffs, m)
end


@inline function Base.:(==)(p1::Polynomial{T, N}, p2::Polynomial{T, N}) where {T, N}
    if p1.is_zero
        iszero(p2.coeffs)
    elseif p2.is_zero
        iszero(p1.coeffs)
    else
        p1.modulus == p2.modulus && p1.coeffs == p2.coeffs
    end
end


@inline function Base.:*(p1::Polynomial{T, N}, p2::Polynomial{T, N}) where {T, N}
    if p1.is_zero || p2.is_zero
        p1
    else
        @assert p1.modulus == p2.modulus
        if p1.modulus == _undefined_modulus
            throw(Exception("Cannot multiply polynomials with undefined modulus"))
        end
        # Making sure `Base.promote_op()` can figure out what the return type is.
        res = p1.mul_function(p1, p2)
        Polynomial{T, N}(res.coeffs, p1.modulus, p1.mul_function)
    end
end

@inline function Base.:*(p1::Polynomial{T, N}, p2::Integer) where {T, N}
    if p1.is_zero
        p1
    else
        Polynomial{T, N}(p1.coeffs .* convert(T, p2), p1.modulus, p1.mul_function)
    end
end

@inline function Base.:*(p1::Integer, p2::Polynomial)
    p2 * p1
end


@inline function Base.div(p1::Polynomial{T, N}, p2::Integer) where {T, N}
    if p1.is_zero
        p1
    else
        Polynomial{T, N}(div.(p1.coeffs, p2), p1.modulus, p1.mul_function)
    end
end


@inline function Base.:+(p1::Polynomial{T, N}, p2::Polynomial{T, N}) where {T, N}
    if p1.is_zero && p2.is_zero
        p1
    elseif p1.is_zero
        p2
    elseif p2.is_zero
        p1
    else
        @assert p1.modulus == p2.modulus
        Polynomial{T, N}(p1.coeffs .+ p2.coeffs, p1.modulus, p1.mul_function)
    end
end

@inline function Base.:+(p1::Polynomial{T, N}, p2::Integer) where {T, N}
    if p1.is_zero
        coeffs = zeros(T, N)
    else
        coeffs = copy(p1.coeffs)
    end
    coeffs[1] += convert(T, p2)
    Polynomial{T, N}(coeffs, p1.modulus, p1.mul_function)
end


@inline function Base.:-(p1::Polynomial{T, N}, p2::Polynomial{T, N}) where {T, N}
    if p1.is_zero && p2.is_zero
        p1
    elseif p1.is_zero
        -p2
    elseif p2.is_zero
        p1
    else
        @assert p1.modulus == p2.modulus
        Polynomial{T, N}(p1.coeffs .- p2.coeffs, p1.modulus, p1.mul_function)
    end
end

@inline function Base.:-(p1::Polynomial{T, N}, p2::Integer) where {T, N}
    if p1.is_zero
        coeffs = zeros(T, N)
    else
        coeffs = copy(p1.coeffs)
    end
    coeffs[1] -= convert(T, p2)
    Polynomial{T, N}(coeffs, p1.modulus, p1.mul_function)
end

@inline function Base.:-(p1::Polynomial{T, N}) where {T, N}
    if p1.is_zero
        p1
    else
        Polynomial{T, N}(.-p1.coeffs, p1.modulus, p1.mul_function)
    end
end


"""
    mul_by_monomial(p::Polynomial, power::Integer)

Multiply the polynomial by `x^power`.
If `power` lies outside `[0, 2 * N)` where `N-1` is the maximum degree of the polynomial,
a modulo `2 * N` will be taken.
"""
@Base.propagate_inbounds @inline function mul_by_monomial(
        p::Polynomial{T, N}, power::Integer) where {T, N}

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

        Polynomial{T, N}(new_coeffs, p.modulus, p.mul_function)
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

    Polynomial{T, N}(r0, p1.modulus, p1.mul_function)
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
    Polynomial{T, N}(c2, p1.modulus, p1.mul_function)
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
    Polynomial{T, N}(res, p1.modulus, p1.mul_function)
end


"""
    resize(p::Polynomial, new_length::Integer)

Returns a polynomial with changed length.
If `new_length` is greater than the current length, coefficients for higher powers will be set to 0.
If `new_length` is smaller, coefficients for higher powers will be discarded.
"""
@inline function resize(p::Polynomial{T, N}, new_length::Integer) where {T, N}
    if new_length > N
        Polynomial{T, new_length}([p.coeffs; zeros(T, new_length - N)], p.modulus)
    elseif new_length < N
        Polynomial{T, new_length}(p.coeffs[1:new_length], p.modulus)
    else
        p
    end
end


#=
Broadcasting machinery.

We want to allow bot expressions like
    p1 .+ [p2, p3, p4, p5] # add a polynomial to an array of polynomials
and
    p1 .+ [1, 2, 3, 4] # add an array to polynomial coefficients
    p2 = func.(p1, some_arg)
(that is both polynomial objects acting as scalars, and as broadcastables).

It seems that the current state of broadcasting in Julia does not allow to resolve this ambiguity.
Therefore, we're marking `Polynomial` as a scalar so that the standard API is invoked
for the former, and adding some custom broadcasting functions
that can be used for the latter.
=#


"""
    broadcast_into_polynomial(func, args...)

Treat any polynomials in `args` as 1-dimensional arrays and apply `func.(args...)`,
saving the result into a new polynomial.

The moduli of the polynomials in `args` must be the same.
"""
function broadcast_into_polynomial(func, args...)
    processed_args = _extract_coeffs.(args)
    new_modulus = _process_polynomials(args)
    if new_modulus === nothing
        new_modulus = cyclic_modulus
    end
    Polynomial(func.(processed_args...), new_modulus)
end


"""
    broadcast_into_polynomial!(func, p::Polynomial, args...)

Treat `p` and any polynomials in `args` as 1-dimensional arrays and apply `p .= func.(args...)`.

The moduli of the polynomials in `args` and the modulus of `p` must be the same.
"""
function broadcast_into_polynomial!(func, p::Polynomial, args...)
    processed_args = _extract_coeffs.(args)
    new_modulus = _process_polynomials(args)
    if !(new_modulus === nothing) && new_modulus != p.modulus
        throw(Exception(
            "Modulus mismatch in polynomial broadcasting: " *
            "target has $(p.modulus), derived $new_modulus"))
    end
    @assert new_modulus == p.modulus
    p.coeffs .= func.(processed_args...)
end


_extract_coeffs(x) = x
_extract_coeffs(x::Polynomial) = x.coeffs


_process_polynomials(args::Tuple) = _join_modulus(
    _extract_modulus(args[1]), _process_polynomials(Base.tail(args)))
_process_polynomials(::Tuple{}) = nothing


_extract_modulus(x::Polynomial) = x.modulus
_extract_modulus(x) = nothing


_join_modulus(x::AbstractPolynomialModulus, y::AbstractPolynomialModulus) =
    if x == y
        x
    else
        throw(Exception("Modulus mismatch in polynomial broadcasting: $x and $y"))
    end
_join_modulus(x::AbstractPolynomialModulus, y) = x
_join_modulus(x, y::AbstractPolynomialModulus) = y
_join_modulus(x, y) = nothing


Base.Broadcast.broadcastable(x::Polynomial) = (x,)


Base.eltype(x::Polynomial{T, N}) where {T, N} = T
