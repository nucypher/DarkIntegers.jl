"""
Arithmetic on opaque unsigned types.
"""


@inline function addhilo(x_hi::T, x_lo::T, y::T) where T <: Unsigned
    r = x_lo + y
    if r < x_lo
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


"""
Multiplication of unsigned integers returning a pair of (high bits, low bits).
"""
@inline function mulhilo_widemul(x::T, y::T) where T <: Unsigned
    # Works for the types for which `widemul()` is defined in the base library
    # (that is, for which a type with double the bitsize exists)
    r = widemul(x, y)
    T(r >> bitsizeof(T)), T(r & typemax(T))
end


# `mulhilo_widemul()` seems to work faster when it's available.
@inline mulhilo(x::T, y::T) where T <: Unsigned = mulhilo_widemul(x, y)

# `widen()` for `UInt128` is `BigInt`, so have to use a more complicated algorithm.
@inline mulhilo(x::UInt128, y::UInt128) = mulhilo_same_type(x, y)


# Divide typemax(T)+1 over `x`
@inline function _div_max_by(x::T) where T
    d, r = divrem(typemax(T), x)
    r += one(T)
    if r == 0
        d += one(T)
    end
    d, r
end


@inline function divhilo(x_hi::T, x_lo::T, y::T) where T <: Unsigned
    # Assuming x_hi < y, otherwise the overflow will be ignored
    res = zero(T)

    a, alpha = _div_max_by(y)
    c, gamma = divrem(x_lo, y)
    res += x_hi * a + c

    hi, lo = mulhilo(x_hi, alpha)
    hi, lo = addhilo(hi, lo, gamma)
    if hi == 0
        res + div(lo, y)
    else
        res + divhilo(hi, lo, y)
    end
end


@inline function divremhilo(x_hi::T, x_lo::T, y::T) where T <: Unsigned
    # assuming x_hi < y
    # TODO: it may be possible to avoid `mulhilo()` and calculate the remainder in `divhilo()`
    q = divhilo(x_hi, x_lo, y)
    _, z = mulhilo(q, y)
    q, x_lo - z
end


@inline function modhilo(x_hi::T, x_lo::T, y::T) where T <: Unsigned
    # assuming x_hi < y
    d, r = divremhilo(x_hi, x_lo, y)
    r
end
