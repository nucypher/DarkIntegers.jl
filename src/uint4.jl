"""
A 4-bit unsigned type to simplify exhaustive testing of arithmetic operations.
"""

struct UInt4 <: Unsigned
    value :: UInt8

    UInt4(x::Integer) = new(x & 0xf)
end


Base.convert(tp::Type{<:Integer}, x::UInt4) = convert(tp, x.value)


Base.promote_type(tp::Type{<:Integer}, ::Type{UInt4}) = tp
Base.promote_type(::Type{UInt4}, tp::Type{<:Integer}) = tp


Base.string(x::UInt4) = uppercase(repr(x.value)[end])


Base.show(io::IO, x::UInt4) = print(io, string(x))


Base.zero(::Type{UInt4}) = UInt4(zero(UInt8))


Base.one(::Type{UInt4}) = UInt4(one(UInt8))


Base.typemin(::Type{UInt4}) = zero(UInt4)


Base.typemax(::Type{UInt4}) = UInt4(0xf)


Base.widen(::Type{UInt4}) = UInt8


Base.leading_zeros(x::UInt4) = leading_zeros(x.value) - 4


Base.:+(x::UInt4, y::UInt4) = UInt4(x.value + y.value)


Base.:-(x::UInt4, y::UInt4) = UInt4(x.value - y.value)
Base.:-(x::UInt4) = zero(UInt4) - x


Base.:*(x::UInt4, y::UInt4) = UInt4(x.value * y.value)


Base.:&(x::UInt4, y::UInt4) = UInt4(x.value & y.value)


Base.:|(x::UInt4, y::UInt4) = UInt4(x.value | y.value)


Base.:<<(x::UInt4, shift::Int) = UInt4(x.value << shift)


Base.:>>(x::UInt4, shift::Int) = UInt4(x.value >> shift)


Base.:<(x::UInt4, y::UInt4) = x.value < y.value


Base.:>(x::UInt4, y::UInt4) = x.value > y.value


Base.:<=(x::UInt4, y::UInt4) = x.value <= y.value


Base.:>=(x::UInt4, y::UInt4) = x.value >= y.value


Base.div(x::UInt4, y::UInt4) = UInt4(div(x.value, y.value))


Base.rem(x::UInt4, y::UInt4) = UInt4(rem(x.value, y.value))


function Base.divrem(x::UInt4, y::UInt4)
    d, r = divrem(x.value, y.value)
    UInt4(d), UInt4(r)
end


# `sizeof` gives the size in bytes, so it is too coarse for UInt4.
# We need the exact size in bits for various functions.
bitsizeof(::Type{UInt4}) = 4
bitsizeof(tp) = 8 * sizeof(tp)


# For mulhilo()
_low_mask(::Type{UInt4}) = UInt4(UInt8(0x3))
_low_shift(tp::Type{UInt4}) = 2
