"""
Modulo arithmetic on opaque unsigned types.
"""


@inline function addmod(x::T, y::T, modulus::T) where T <: Unsigned
    r = x + y
    if r < x || r >= modulus
        r - modulus
    else
        r
    end
end


@inline function submod(x::T, y::T, modulus::T) where T <: Unsigned
    r = x - y
    if x < y
        r + modulus
    else
        r
    end
end


@inline isone(x::T) where T <: Real = one(T) === x


@inline double(x::Unsigned) = x << 1


@inline halve(x::Unsigned) = x >> 1


@inline function mulmod_bitshift(x::T, y::T, modulus::T) where T <: Unsigned

    if iszero(x) || isone(y)
        return x
    end

    if iszero(y) || isone(x)
        return y
    end

    result  = zero(T)

    while !iszero(y)
        if isodd(y)
            result = addmod(result, x, modulus)
        end

        x = double(x)
        t = x
        if x < t || x >= modulus
            x = x - modulus
        end

        y = halve(y)
    end

    return result
end


@inline function mulmod_modhilo(x::T, y::T, modulus::T) where T <: Unsigned
    hi, lo = mulhilo(x, y)
    modhilo(hi, lo, modulus)
end


@inline function mulmod_widemul(x::T, y::T, modulus::T) where T <: Unsigned
    T(mod(widemul(x, y), widen(T)(modulus)))
end


@inline mulmod(x::T, y::T, modulus::T) where T <: Unsigned =
    mulmod_widemul(x, y, modulus)
@inline mulmod(x::T, y::T, modulus::T) where T <: Union{UInt64, UInt128} =
    mulmod_bitshift(x, y, modulus)
