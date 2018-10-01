using Base: setindex


function get_montgomery_coeff(modulus::MPNumber{N, T}) where {N, T}
    # calculate -m^(-1) mod b, where b = typemax(T)+1
    get_montgomery_coeff(modulus[1])
end


function get_montgomery_coeff(modulus::T) where T <: Unsigned
    # calculate -m^(-1) mod b, where b = typemax(T)+1
    -T(invmod(BigInt(modulus), BigInt(typemax(T)) + 1))
end


# `res += c * v`, where `res` is a multi-precision number of length `n+1`,
# `v` is a multi-precision number of length `n`, and `c` is a single radix digit.
# Modifies `res` and returns the carry bit (if there's an overflow in the `n+1`-th digit).
@inline function _mul_addto(
        res::MPNumber{N, T}, res_hi::T, c::T, v::MPNumber{N, T}) where {N, T}

    hi_carry1 = false
    hi_carry2 = false
    hi_carry3 = false
    hi = zero(T)
    out = zero(MPNumber{N, T})
    for j in 1:N
        new_hi, lo = mulhilo(c, v[j])

        r, new_hi_carry1 = _addc(res[j], lo)
        r, new_hi_carry2 = _addc(r, hi)
        r, new_hi_carry3 = _addc(r, T(hi_carry1 + hi_carry2 + hi_carry3))

        out = setindex(out, r, j)
        hi = new_hi
        hi_carry1 = new_hi_carry1
        hi_carry2 = new_hi_carry2
        hi_carry3 = new_hi_carry3
    end

    out_hi, c1 = _addc(res_hi, hi)
    out_hi, c2 = _addc(out_hi, T(hi_carry1 + hi_carry2 + hi_carry3))
    out, out_hi, c1 || c2 # TODO: it seems that we can have at most 1 carried over to the n+2-th digit
end


@inline @generated function _mul_addto_g(
        res::MPNumber{N, T}, res_hi::T, c::T, v::MPNumber{N, T}) where {N, T}

    loop_body = []
    for j in 1:N
        push!(loop_body, quote
            new_hi, lo = mulhilo(c, v[$j])

            r, new_hi_carry1 = _addc(res[$j], lo)
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

        out_hi, c1 = _addc(res_hi, hi)
        out_hi, c2 = _addc(out_hi, T(hi_carry1 + hi_carry2 + hi_carry3))

        # TODO: it seems that we can have at most 1 carried over to the n+2-th digit
        out, out_hi, c1 || c2
    end
end



# Montgomery multiplication algorithm for multi-precision numbers.
@Base.propagate_inbounds @inline function mulmod_montgomery(
        x::MPNumber{N, T}, y::MPNumber{N, T}, m::MPNumber{N, T}, m_prime::T) where {N, T}
    """
    Montgomery multiplication: for `x = x' * R mod m` and `y = y' * R mod m`
    calculate `x' * y' * R mod m`, where `R = typemax(MPNumber{N, T}) + 1`.
    `m_prime` is the Montgomery coefficient (see `get_montgomery_coeff()`).
    """

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


@inline function mulmod_montgomery(x::T, y::T, m::T, m_prime::T) where T <: Unsigned
    mulmod_montgomery(MPNumber((x,)), MPNumber((y,)), MPNumber((m,)), m_prime)[1]
end


# Calculate `R mod m` for use in `to_montgomery()`
function get_to_montgomery_coeff(modulus::T) where T <: Unsigned
    T(mod(convert(BigInt, typemax(T)) + 1, convert(BigInt, modulus)))
end


# Given `x`, return `x * R mod m`,
# where `R = b^n`, `b` is the digit size, and `n` is the number length.
@inline function to_montgomery(x::T, m::T, coeff::T) where T <: Unsigned
    mulmod(x, coeff, m)
end


# Given `x`, return `x / R mod m` (Montgomery reduction),
# where `R = b^n`, `b` is the digit size, and `n` is the number length.
# Obtained from `mulmod_montgomery()` by setting `x=1` and renaming `y` to `x`.
@Base.propagate_inbounds @inline function from_montgomery(
        x::MPNumber{N, T}, m::MPNumber{N, T}, m_prime::T) where {N, T}
    mulmod_montgomery(one(MPNumber{N, T}), x, m, m_prime)
end


@Base.propagate_inbounds @inline function from_montgomery(x::T, m::T, m_prime::T) where T <: Unsigned
    from_montgomery(MPNumber((x,)), MPNumber((m,)), m_prime)[1]
end
