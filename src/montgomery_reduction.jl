"""
Underlying algorithms for Montgomery representation of integers.
"""

using Base: setindex


"""
Calculate the coefficient used for multiplication and conversion from M. representation.
Namely, `m' = -m^(-1) mod b`, where `b = typemax(T)+1`.
"""
function get_montgomery_coeff(modulus::T) where T <: Unsigned
    # have to widen the type, since we need to use `typemax(T)+1`.
    T2 = widen(T)
    -T(invmod(convert(T2, modulus), convert(T2, typemax(T)) + 1))
end


"""
An implementation of Montgomery coefficient for a multi-precision number.
Using the fact that `-m^(-1) mod R == -(mod(m, b+1))^(-1) mod b`,
where `b = typemax(T)+1` is the radix, and `R = b^N = typemax(MPNumber{N, T})+1`
is the total size of the type.
"""
function get_montgomery_coeff(modulus::MPNumber{N, T}) where {N, T}
    get_montgomery_coeff(modulus[1])
end


"""
Calculates `x + c * v`, where `x` is a multi-precision number of length `N+1`
(passed as a number of length `N` + the highmost digit),
`v` is a multi-precision number of length `N`, and `c` is a single radix digit.
Returns a tuple of the result and the carry bit (if there's an overflow in the `N+1`-th digit).
"""
@inline @generated function _mul_addto_g(
        x::MPNumber{N, T}, x_hi::T, c::T, v::MPNumber{N, T}) where {N, T}

    # A generated version is a bit uglier, but works noticeably faster.

    loop_body = []
    for j in 1:N
        push!(loop_body, quote
            new_hi, lo = mulhilo(c, v[$j])

            r, new_hi_carry1 = _addc(x[$j], lo)
            r, new_hi_carry2 = _addc(r, hi)
            r, new_hi_carry3 = _addc(r, T(hi_carry1 + hi_carry2 + hi_carry3))

            out = setindex(out, r, $j)
            hi = new_hi
            hi_carry1 = new_hi_carry1
            hi_carry2 = new_hi_carry2
            hi_carry3 = new_hi_carry3
        end)
    end

    quote
        hi_carry1 = false
        hi_carry2 = false
        hi_carry3 = false
        hi = zero(T)
        out = zero(MPNumber{N, T})

        $(loop_body...)

        out_hi, c1 = _addc(x_hi, hi)
        out_hi, c2 = _addc(out_hi, T(hi_carry1 + hi_carry2 + hi_carry3))

        # TODO: it seems that we can have at most 1 carried over to the n+2-th digit
        out, out_hi, c1 || c2
    end
end


"""
Montgomery multiplication (or Montgomery reduction algorithm).
For `x = x' * R mod m` and `y = y' * R mod m`
calculates `x' * y' * R mod m`, where `R = typemax(MPNumber{N, T}) + 1`.
`m_prime` is the Montgomery coefficient (see [`get_montgomery_coeff`](@ref)).
"""
@Base.propagate_inbounds @inline function mulmod_montgomery(
        x::MPNumber{N, T}, y::MPNumber{N, T}, m::MPNumber{N, T}, m_prime::T) where {N, T}

    a = zero(MPNumber{N, T})
    a_hi = zero(T)
    for i in 1:N

        u = (a[1] + x[i] * y[1]) * m_prime

        # A = A + x[i] * y + u * m
        a, a_hi, c1 = _mul_addto_g(a, a_hi, x[i], y)
        a, a_hi, c2 = _mul_addto_g(a, a_hi, u, m)

        # A = A / b
        for j in 1:N-1
            a = setindex(a, a[j+1], j)
        end
        a = setindex(a, a_hi, N)

        # TODO: is it actually ever greater than 0? It seems to always be in tests.
        # Setting it to just zero improves the performance a lot.
        a_hi = T(c1 + c2)
    end

    # It can be proven that `a` is at most `2m - 1`.
    if a_hi > 0 || a >= m
        a - m
    else
        a
    end
end


"""
An implementation of Montgomery reduction for simple types.
Treats them as `MPNumber`s of length 1.
"""
@inline function mulmod_montgomery(x::T, y::T, m::T, m_prime::T) where T <: Unsigned
    mulmod_montgomery(MPNumber((x,)), MPNumber((y,)), MPNumber((m,)), m_prime)[1]
end


"""
Calculates the coefficient for conversion of an integer to Montgomery representation.
The conversion formula is `x' = x * R mod m`, but `R = typemax(T) + 1`
does not fit in the type of `x`.
So we have to calculate `c = R mod m` and use it as `x' = x * c mod m`.
"""
function get_to_montgomery_coeff(m::T) where T <: Unsigned
    # have to widen the type, since we need to use `typemax(T)+1`.
    T2 = widen(T)
    T(mod(convert(T2, typemax(T)) + 1, convert(T2, m)))
end


"""
Converts an integer to Montgomery representation
(that is, calculates `x * R mod m == x * coeff mod m` where `coeff = R mod m`,
where `R = typemax(T) + 1`).
"""
@inline function to_montgomery(x::T, m::T, coeff::T) where T <: Unsigned
    mulmod(x, coeff, m)
end


"""
Recovers an integer from Montgomery representation
(that is, calculates `x` given `x * R mod m`, where `R = typemax(MPNumber{N, T}) + 1`).
"""
@Base.propagate_inbounds @inline function from_montgomery(
        x::MPNumber{N, T}, m::MPNumber{N, T}, m_prime::T) where {N, T}
    # Montgomery multiplication of `1` and `x` in M. representation (`x * R`)
    # results in `1 * (x * R) / R = x`.
    mulmod_montgomery(one(MPNumber{N, T}), x, m, m_prime)
end


@Base.propagate_inbounds @inline function from_montgomery(x::T, m::T, m_prime::T) where T <: Unsigned
    from_montgomery(MPNumber((x,)), MPNumber((m,)), m_prime)[1]
end
