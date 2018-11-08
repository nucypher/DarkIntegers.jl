"""
Modulo arithmetic on opaque unsigned types.
"""


"""
    addmod(x::T, y::T, modulus::T) where T <: Unsigned

Returns `mod(x + y, modulus)` (even if `x + y` overflows `T`).
"""
@inline function addmod(x::T, y::T, modulus::T) where T <: Unsigned
    r = x + y
    if r < x || r >= modulus
        r - modulus
    else
        r
    end
end


"""
    submod(x::T, y::T, modulus::T) where T <: Unsigned

Returns `mod(x - y, modulus)` (even if `x - y` underflows).
"""
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


@inline function doublemod(x::T, modulus::T) where T <: Unsigned
    # Could be done as `addmod(x, x)`, but this variant is faster.
    t = x
    x = double(x)
    if x < t || x >= modulus
        x = x - modulus
    end
end


@inline halve(x::Unsigned) = x >> 1


"""
Modulo multiplication using bitshift (by 1 bit) and modulo addition/subtraction.
"""
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
        x = doublemod(x, modulus)
        y = halve(y)
    end

    return result
end


"""
Modulo multiplication using wide multiplication into separate hi/lo variables.
"""
@inline function mulmod_remhilo(x::T, y::T, modulus::T) where T <: Unsigned
    hi, lo = mulhilo(x, y)
    remhilo(hi, lo, modulus)
end


"""
Modulo multiplication using wide multiplication into a single bigger type.
"""
@inline function mulmod_widemul(x::T, y::T, modulus::T) where T <: Unsigned
    T(mod(widemul(x, y), convert(widen(T), modulus)))
end


@doc """
    mulmod(x::T, y::T, modulus::T) where T <: Unsigned

Returns `mod(x * y, modulus)` (even if `x * y` overflows `T`).
""" mulmod()

@inline mulmod(x::T, y::T, modulus::T) where T <: Unsigned =
    mulmod_widemul(x, y, modulus)
@inline mulmod(x::T, y::T, modulus::T) where T <: Union{UInt64, UInt128} =
    mulmod_bitshift(x, y, modulus)
