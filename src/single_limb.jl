"""
Arithmetic on opaque unsigned types.
"""


"""
Returns the pair `(hi, lo)` of `x + y + z`.
That is, `hi * R + lo = x + y + z`,
where `0 <= x, y, z, lo < R`, `hi <= 2`, and `R = typemax(T) + 1`.

If `z <= 1`, `hi <= 1` (so it can be used as `new_carry, res = addcarry(x, y, carry)`).
"""
@inline function addcarry(x::T, y::T, z::T) where T <: Unsigned
    T2 = widen(T)
    res = (x % T2) + (y % T2) + (z % T2)
    (res >> bitsizeof(T)) % T, res % T
end


"""
A specialized version for `z == 0`.
"""
@inline function addcarry(x::T, y::T) where T <: Unsigned
    T2 = widen(T)
    ret = (x % T2) + (y % T2)
    (ret >> bitsizeof(T)) % T, ret % T
end


"""
Since `widen(UInt128) == BigInt`, we need a specialized version to avoid allocations.
"""
@inline function addcarry(x::UInt128, y::UInt128, z::UInt128)
    x_lo = x % UInt64
    x_hi = (x >> 64) % UInt64
    y_lo = y % UInt64
    y_hi = (y >> 64) % UInt64
    z_lo = z % UInt64
    z_hi = (z >> 64) % UInt64

    t1, r1 = addcarry(x_lo, y_lo, z_lo) # t1 <= 1
    t2, r2 = addcarry(x_hi, y_hi, t1) # t2 <= 1
    t3, r2 = addcarry(r2, z_hi) # t3 <= 1

    lo = (r1 % UInt128) | ((r2 % UInt128) << 64)
    hi = (t2 + t3) % UInt128 # t2, t3 <=1 -> no overflow here

    hi, lo
end


@inline function addcarry(x::UInt128, y::UInt128)
    x_lo = x % UInt64
    x_hi = (x >> 64) % UInt64
    y_lo = y % UInt64
    y_hi = (y >> 64) % UInt64

    t1, r1 = addcarry(x_lo, y_lo)
    t2, r2 = addcarry(x_hi, y_hi, t1)

    lo = (r1 % UInt128) | ((r2 % UInt128) << 64)
    hi = t2 % UInt128

    hi, lo
end


"""
Returns the pair `(new_borrow, res)` of `x - (y + borrow >> (sizeof(T) - 1))`.
`borrow` and `new_borrow` can be either `0` or `typemax(T)`,
the latter meaning that the borrow occurred during subtraction.

Note that it is not an analogue of `addhilo` for subtraction.
"""
@inline function subborrow(x::T, y::T, borrow::T) where T <: Unsigned
    T2 = widen(T)
    ret = (x % T2) - ((y % T2) + ((borrow >> (bitsizeof(T) - 1)) % T2))
    (ret >> bitsizeof(T)) % T, ret % T
end


@inline function subborrow(x::UInt128, y::UInt128, borrow::UInt128)
    x_lo = x % UInt64
    x_hi = (x >> 64) % UInt64
    y_lo = y % UInt64
    y_hi = (y >> 64) % UInt64
    b_lo = borrow % UInt64 # borrow is always 0 or typemax()

    nb1, r1 = subborrow(x_lo, y_lo, b_lo)
    nb2, r2 = subborrow(x_hi, y_hi, nb1)

    res = (r1 % UInt128) | ((r2 % UInt128) << 64)
    new_borrow = (nb2 % UInt128) | ((nb2 % UInt128) << 64)

    new_borrow, res
end


@inline function subborrow(x::T, y::T) where T <: Unsigned
    T2 = widen(T)
    ret = (x % T2) - (y % T2)
    (ret >> bitsizeof(T)) % T, ret % T
end


@inline function subborrow(x::UInt128, y::UInt128)
    x_lo = x % UInt64
    x_hi = (x >> 64) % UInt64
    y_lo = y % UInt64
    y_hi = (y >> 64) % UInt64

    nb1, r1 = subborrow(x_lo, y_lo)
    nb2, r2 = subborrow(x_hi, y_hi, nb1)

    res = (r1 % UInt128) | ((r2 % UInt128) << 64)
    new_borrow = (nb2 % UInt128) | ((nb2 % UInt128) << 64)

    new_borrow, res
end


"""
Returns the pair `(hi, lo)` of `x + y * z + w`.
That is, `hi * R + lo = x + y * z + w`,
where `0 <= x, y, z, w, hi, lo < R`, and `R = typemax(T) + 1`.
"""
@inline function muladdcarry(x::T, y::T, z::T, w::T) where T
    hi, lo = mulhilo(y, z)
    t, r1 = addcarry(x, lo, w)
    _, r2 = addcarry(hi, t)
    r2, r1
end


@inline function muladdcarry(x::T, y::T, z::T) where T
    hi, lo = mulhilo(y, z)
    t, r1 = addcarry(x, lo)
    _, r2 = addcarry(hi, t)
    r2, r1
end


# Addition of unsigned numbers with carry
@inline function _addc(x::T, y::T) where T <: Unsigned
    # Base.Checked.add_with_overflow is slower
    r = x + y
    r, r < x
end


# Subtraction of unsigned numbers with carry
@inline function _subc(x::T, y::T) where T <: Unsigned
    # Base.Checked.sub_with_overflow is slower
    r = x - y
    r, x < y
end


# Multiplication of unsigned numbers with remembering overflow
@inline function _mulc(x::T, y::T) where T <: Unsigned
    Base.Checked.mul_with_overflow(x, y)
end

@inline function _mulc(x::UInt4, y::UInt4) where T <: Unsigned
    r = x.value * y.value
    UInt4(r & 0xf), r > 0xf
end


"""
    addhilo(x_hi::T, x_lo::T, y::T) where T <: Unsigned

Calculates `hi * B + lo = x_hi * B + x_lo + y`, where `B == typemax(T) + 1`.
Returns the result as a pair `(hi::T, lo::T)`.
An overflow in `hi` (if any) is ignored.
"""
@inline function addhilo(x_hi::T, x_lo::T, y::T) where T <: Unsigned
    r, c = _addc(x_lo, y)
    if c
        x_hi += one(T)
    end
    x_hi, r
end


# The mask for the lower half of the type's bits
_low_mask(::Type{UInt8}) = UInt8(0xf)
_low_mask(::Type{UInt16}) = UInt16(0xff)
_low_mask(::Type{UInt32}) = UInt32(0xffff)
_low_mask(::Type{UInt64}) = UInt64(0xffffffff)
_low_mask(::Type{UInt128}) = UInt128(0xffffffffffffffff)


_low_shift(tp::Type{<:Unsigned}) = sizeof(tp) * 4


# Adopted from one of the methods of `widemul()` from Julia base library
@inline function mulhilo_same_type(x::T, y::T) where T <: Unsigned

    m = _low_mask(T)
    shift = _low_shift(T)

    x0 = x & m
    x1 = x >> shift
    y0 = y & m
    y1 = y >> shift

    w0 = x0 * y0
    t = x1 * y0 + (w0 >> shift)
    w2 = t >> shift
    w1 = x0 * y1 + (t & m)

    hi = x1 * y1 + w2 + (w1 >> shift)
    lo = w0 & m + (w1 << shift)

    hi, lo
end


# Multiplication using conversion to a larger type
@inline function mulhilo_widemul(x::T, y::T) where T <: Unsigned
    # Works for the types for which `widemul()` is defined in the base library
    # (that is, for which a type with double the bitsize exists)
    r = widemul(x, y)
    (r >> bitsizeof(T)) % T, r % T
end


@doc """
    mulhilo(x::T, y::T) where T <: Unsigned

Calculates `hi * B + lo = x * y`, where `B == typemax(T) + 1`.
Returns the result as a pair `(hi::T, lo::T)`.
""" mulhilo()

# `mulhilo_widemul()` seems to work faster when it's available.
@inline mulhilo(x::T, y::T) where T <: Unsigned = mulhilo_widemul(x, y)

# `widen()` for `UInt128` is `BigInt`, so have to use a more complicated algorithm.
@inline mulhilo(x::UInt128, y::UInt128) = mulhilo_same_type(x, y)


@inline function divremhilo_widen(x_hi::T, x_lo::T, y::T) where T <: Unsigned
    T2 = widen(T)
    x_hi_w = x_hi % T2
    x_lo_w = x_lo % T2

    x_w = (x_hi_w << bitsizeof(T)) | x_lo_w
    y_w = y % T2

    q_w, r_w = divrem(x_w, y_w)
    q_lo = q_w & (typemax(T) % T2)
    overflow = q_lo != q_w
    q_lo % T, r_w % T, overflow
end


# Divide typemax(T)+1 over `x`
@inline function _div_max_by(x::T) where T
    d, r = divrem(typemax(T), x)
    r += one(T)
    if r == 0
        d += one(T)
    end
    d, r
end


@inline function divremhilo_same_type(x_hi::T, x_lo::T, y::T) where T <: Unsigned

    # b = y * a + alpha
    a, alpha = _div_max_by(y)

    # x_lo = y * c + gamma
    c, gamma = divrem(x_lo, y)

    # x = x_hi * b + x_lo
    #   = x_hi * (y * a + alpha) + y * c + gamma
    #   = y * (x_hi * a + c) + x_hi * alpha + gamma

    q, o1 = _mulc(x_hi, a)
    q, o2 = _addc(q, c)
    overflow = o1 | o2

    hi, lo = mulhilo(x_hi, alpha)
    hi, lo = addhilo(hi, lo, gamma)
    if hi == 0
        q2, r = divrem(lo, y)
        o3 = false
    else
        q2, r, o3 = divremhilo_same_type(hi, lo, y)
    end

    q, o4 = _addc(q, q2)

    q, r, overflow | o3 | o4
end


@doc """
    divremhilo(x_hi::T, x_lo::T, y::T) where T <: Unsigned

Calculates `divrem(x_hi * B + x_lo, y)`, where `B == typemax(T) + 1`.
Returns a tuple `(q::T, r::T, o::Bool)`, where `q` is the quotient
(the part of it fitting into the bits of `T`), `r` is the remainder is `o` is the overflow flag.
""" divremhilo()


# `mulhilo_widemul()` seems to work faster when it's available.
@inline divremhilo(x_hi::T, x_lo::T, y::T) where T <: Unsigned = divremhilo_widen(x_hi, x_lo, y)

# `widen()` for `UInt128` is `BigInt`, so have to use a more complicated algorithm.
@inline divremhilo(x_hi::UInt128, x_lo::UInt128, y::UInt128) = divremhilo_same_type(x_hi, x_lo, y)


"""
    divhilo(x_hi::T, x_lo::T, y::T) where T <: Unsigned

Calculates `div(x_hi * B + x_lo, y)`, where `B == typemax(T) + 1`.
Returns a tuple `(q::T, o::Bool)`, where `q` is the quotient
(the part of it fitting into the bits of `T`), and `o` is the overflow flag.
"""
@inline function divhilo(x_hi::T, x_lo::T, y::T) where T <: Unsigned
    q, r, o = divremhilo(x_hi, x_lo, y)
    q, o
end


"""
    remhilo(x_hi::T, x_lo::T, y::T) where T <: Unsigned

Calculates `rem(x_hi * B + x_lo, y)`, where `B == typemax(T) + 1`.
"""
@inline function remhilo(x_hi::T, x_lo::T, y::T) where T <: Unsigned
    q, r, o = divremhilo(x_hi, x_lo, y)
    r
end
