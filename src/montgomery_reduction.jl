"""
Underlying algorithms for Montgomery representation of integers.

The algorithm is slightly different to the classical one (found e.g. in Knuth):
instead of `-m^(-1) mod b` as the coefficient, we have just `m^(-1) mod b`.
Consequently, this leads to a plus replaced with minus at one stage during multiplication,
and one less condition and addition for single-limb numbers.

The idea is that during Montgomery reduction or multiplication,
in the sum where the coefficient is used, the lowest limb gets turned to 0 during an addition
(due to the choice of the coefficient).
If `x+y mod b == 0`, `x+y` can be either `0` or `b`,
so we have to check for carry and add it to the result.
With our new coefficient, we have `x-y mod b == 0`, so there's no carry needed.

This leads to 4x performance increase for arrays of UInt64, for example,
because the new version can be vectorized much more efficiently.
"""

using Base: setindex


"""
Calculate the coefficient used for multiplication and conversion from M. representation.
Namely, `m' = m^(-1) mod b`, where `b = typemax(T)+1`.
"""
function get_montgomery_coeff(modulus::T) where T <: Unsigned
    # have to widen the type, since we need to use `typemax(T)+1`.
    T2 = widen(T)
    T(invmod(convert(T2, modulus), convert(T2, typemax(T)) + 1))
end


"""
An implementation of Montgomery coefficient for a multi-precision number.
Using the fact that `m^(-1) mod R == (mod(m, b+1))^(-1) mod b`,
where `b = typemax(T)+1` is the radix, and `R = b^N = typemax(MPNumber{N, T})+1`
is the total size of the type.
"""
function get_montgomery_coeff(modulus::MPNumber{N, T}) where {N, T}
    get_montgomery_coeff(modulus[1])
end


"""
Calculates `x * y`, where `x` is a multi-precision number of length `N`,
and `y` is a single radix digit.
Returns the `N+1`-long result as tuple of the lowest digit
and a multi-precision number with the `N` highest digits.
"""
@inline @generated function _mul_by_single(x::MPNumber{N, T}, y::T) where {N, T}

    body = []

    for j in 1:N
        push!(body, quote
            hi, lo = mulhilo(x[$j], y)
        end)

        if j > 1
            push!(body, quote
                hi, lo = addhilo(hi, lo, c)
                hi, lo = addhilo(hi, lo, w[$(j-1)])
            end)
        end

        push!(body, quote
            $(j == 1 ? :(w_lo = lo) : :(w = setindex(w, lo, $(j-1))))
            c = hi
        end)
    end

    quote
        w = zero(MPNumber{N, T})
        c = zero(T)
        w_lo = zero(T)

        $(body...)

        w = setindex(w, c, N)
        w_lo, w
    end
end


"""
Montgomery multiplication (or Montgomery reduction algorithm).
For `x = x' * R mod m` and `y = y' * R mod m`
calculates `x' * y' * R mod m`, where `R = typemax(MPNumber{N, T}) + 1`.
`m_prime` is the Montgomery coefficient (see [`get_montgomery_coeff`](@ref)).
"""
@Base.propagate_inbounds @inline @generated function mulmod_montgomery(
        x::MPNumber{N, T}, y::MPNumber{N, T}, m::MPNumber{N, T}, m_prime::T) where {N, T}

    body = []

    for i in 1:N
        push!(body, quote
            p1_lo, p1 = _mul_by_single(y, x[$i])
            t = p1_lo
        end)

        if i > 1
            push!(body, quote
                t += a_lo
                a += p1
            end)
        else
            push!(body, quote
                a = p1
            end)
        end

        push!(body, quote
            u = t * m_prime
            p2_lo, p2 = _mul_by_single(m, u)
        end)

        # a_lo + p1_lo - p2_lo is always zero mod radix (sizeof(T)+1)
        # We just want to know if we need to carry 1 or not to the higher limbs.
        # If i == 1, a_lo = 0, so the carry is also 0.
        if i > 1
            push!(body, quote
                if (a_lo + p1_lo) < a_lo
                    a += one(MPNumber{N, T})
                end
            end)
        end

        push!(body, quote
            # TODO: can `a` be greater than `m` at this point?
            a = submod(a, p2, m)
        end)

        if i == N
            push!(body, quote
                a
            end)
        else
            new_a = [:(a[$j]) for j in 2:N]
            push!(body, quote
                a_lo = a[1]
                a = MPNumber{N, T}(($(new_a...), zero(T)))
            end)
        end
    end

    quote
        $(body...)
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
"""
function get_to_montgomery_coeff(m::T) where T <: Unsigned
    #=
    The conversion formula is `x' = x * R mod m`.
    We want to use the quick Montgomery multiplication, so we precalculate `R^2 mod m`,
    and later use it as `x' = montgomery_mul(x, R^2 mod m) = x * R^2 * R^(-1) mod m = x * R mod m`.
    (here R is the container type size)
    =#
    T2 = widen(T)
    R = convert(T2, typemax(T)) + 1
    m2 = convert(T2, m)
    R_mod_m = mod(R, m2)
    T(mod(R_mod_m^2, m2))
end


"""
Converts an integer to Montgomery representation
(that is, calculates `x * R mod m == x * coeff mod m` where `coeff = R mod m`,
where `R = typemax(T) + 1`).
"""
@inline function to_montgomery(
        x::MPNumber{N, T}, m::MPNumber{N, T}, m_prime::T, coeff::MPNumber{N, T}) where {N, T}
    mulmod_montgomery(x, coeff, m, m_prime)
end

@inline function to_montgomery(x::T, m::T, m_prime::T, coeff::T) where T <: Unsigned
    mulmod_montgomery(x, coeff, m, m_prime)
end


"""
Recovers an integer from Montgomery representation
(that is, calculates `x` given `x * R mod m`, where `R = typemax(MPNumber{N, T}) + 1`).
"""
@inline function from_montgomery(x::MPNumber{N, T}, m::MPNumber{N, T}, m_prime::T) where {N, T}
    # Montgomery multiplication of `1` and `x` in M. representation (`x * R`)
    # results in `1 * (x * R) / R = x`.
    mulmod_montgomery(one(MPNumber{N, T}), x, m, m_prime)
end

@inline function from_montgomery(x::T, m::T, m_prime::T) where T <: Unsigned
    mulmod_montgomery(one(T), x, m, m_prime)
end
