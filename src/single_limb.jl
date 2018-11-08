"""
Arithmetic on opaque unsigned types.
"""


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
    T(r >> bitsizeof(T)), T(r & typemax(T))
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
    x_hi_w = convert(T2, x_hi)
    x_lo_w = convert(T2, x_lo)

    x_w = (x_hi_w << bitsizeof(T)) | x_lo_w
    y_w = convert(T2, y)

    q_w, r_w = divrem(x_w, y_w)
    q_lo = q_w & convert(T2, typemax(T))
    overflow = q_lo != q_w
    T(q_lo), T(r_w), overflow
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
