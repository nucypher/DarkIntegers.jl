"""
Montgomery multiplication algorithms in constant time.

Non-constant time is faster for 1-2 limbs, but by 3-4 limbs the difference evens out.

TODO: at the moment these functions are not joined with high-level types;
      this will be done as a part of ConstantTime integration.
"""


"""
Calculate the coefficient used for multiplication and conversion from M. representation.
Namely, `m' = -m^(-1) mod b`, where `b = typemax(T)+1`.

Note that this coefficient is different from what the non-constant-time version uses.
"""
function get_montgomery_coeff_ct(modulus::T) where T <: Unsigned
    # have to widen the type, since we need to use `typemax(T)+1`.
    T2 = widen(T)
    #invmod(modulus % T2, typemax(T) % T2 + 1) % T
    -invmod(modulus % T2, typemax(T) % T2 + one(T2)) % T
end


"""
An implementation of Montgomery coefficient for a multi-precision number.
Using the fact that `m^(-1) mod R == (mod(m, b+1))^(-1) mod b`,
where `b = typemax(T)+1` is the radix, and `R = b^N = typemax(MLUInt{N, T})+1`
is the total size of the type.
"""
function get_montgomery_coeff_ct(modulus::MLUInt{N, T}) where {N, T}
    get_montgomery_coeff_ct(modulus[1])
end


"""
Returns `x + y * z`, where `x` and `y` are multi-limb integers.
"""
@inline function muladd_ml(x::MLUInt{N, T}, y::MLUInt{N, T}, z::T) where {N, T}
    res = zero(MLUInt{N, T})

    hi, lo = muladdcarry(x[1], y[1], z)
    res = setindex(res, lo, 1)
    for i in 2:N
        hi, lo = muladdcarry(x[i], y[i], z, hi)
        res = setindex(res, lo, i)
    end
    hi, res
end


# A specialized version of `muladd` for `x = 0`.
@inline function muladd_ml(y::MLUInt{N, T}, z::T) where {N, T}
    res = zero(MLUInt{N, T})

    hi, lo = mulhilo(y[1], z)
    res = setindex(res, lo, 1)
    for i in 2:N
        hi, lo = muladdcarry(hi, y[i], z)
        res = setindex(res, lo, i)
    end
    hi, res
end


"""
And extended version of `subborrow()` for two multi-limb integers `(x..., x_hi)`
and `(y..., y_hi)` (`x` and `y` represent lower limbs of the corresponding numbers).
Returned `borrow` is the same as for the single-limb `subborrow()`.
"""
@inline function subborrow_ml(x_hi::T, x::MLUInt{N, T}, y_hi::T, y::MLUInt{N, T}) where {N, T}
    res = zero(MLUInt{N, T})
    borrow, lo = subborrow(x[1], y[1])
    res = setindex(res, lo, 1)
    for i in 2:N
        borrow, lo = subborrow(x[i], y[i], borrow)
        res = setindex(res, lo, i)
    end
    borrow, _ = subborrow(x_hi, y_hi, borrow)
    borrow, res
end


"""
A multi-limb version of `addcarry`.
"""
@inline function addcarry_ml(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    res = zero(MLUInt{N, T})
    hi, lo = addcarry(x[1], y[1])
    res = setindex(res, lo, 1)
    for i in 2:N
        hi, lo = addcarry(x[i], y[i], hi)
        res = setindex(res, lo, i)
    end
    hi, res
end


"""
Returns `mod((x..., x_hi) - (y..., y_hi), m)`.
Assumes `-m <= (x..., x_hi) - (y..., y_hi) < m`.
"""
@inline function submod_inner(
        x_hi::T, x::MLUInt{N, T}, y_hi::T, y::MLUInt{N, T}, m::MLUInt{N, T}) where {N, T}

    borrow, w = subborrow_ml(x_hi, x, y_hi, y)

    masked_m = m
    for i in 1:N
        masked_m = setindex(masked_m, m[i] & borrow, i)
    end

    # If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
    # borrow = 0x000...000. Thus, we use it as a mask to conditionally add the
    # modulus.
    _, res = addcarry_ml(w, masked_m)

    res
end



@Base.propagate_inbounds function mulmod_montgomery_ct(
        x::MLUInt{N, T}, y::MLUInt{N, T}, m::MLUInt{N, T}, m_p::T) where {N, T}

    u = (x[1] * y[1]) * m_p
    a_hi1, a_lo = muladd_ml(m, u)
    a_hi2, a_lo = muladd_ml(a_lo, y, x[1])
    hi, lo = addcarry(a_hi1, a_hi2)

    for j in 1:N-1
        a_lo = setindex(a_lo, a_lo[j+1], j)
    end
    a_lo = setindex(a_lo, lo, N)
    a_hi = hi

    for i in 2:N
        u = (a_lo[1] + x[i] * y[1]) * m_p
        a_hi1, a_lo = muladd_ml(a_lo, m, u)
        a_hi2, a_lo = muladd_ml(a_lo, y, x[i])
        hi, lo = addcarry(a_hi, a_hi1, a_hi2)

        for j in 1:N-1
            a_lo = setindex(a_lo, a_lo[j+1], j)
        end
        a_lo = setindex(a_lo, lo, N)
        a_hi = hi
    end

    submod_inner(a_hi, a_lo, zero(T), m, m)

end


@inline function mulmod_montgomery_ct(x::T, y::T, m::T, m_prime::T) where T <: Unsigned
    mulmod_montgomery_ct(MLUInt((x,)), MLUInt((y,)), MLUInt((m,)), m_prime)[1]
end
