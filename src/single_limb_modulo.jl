"""
Modulo arithmetic on opaque unsigned types.
"""


"""
    addmod(x::T, y::T, modulus::T) where T <: Unsigned

Returns `mod(x + y, modulus)` (even if `x + y` overflows `T`).
Assumes `x` and `y` are in range `0 ... modulus`.
"""
@inline function addmod(x::T, y::T, modulus::T) where T <: Unsigned
    r = x + y
    if r < x || !(r < modulus)
        r - modulus
    else
        r
    end
end


# Check that `x+y`, where `0 <= x, y < modulus` can't overflow the container type.
function addition_cant_overflow(modulus::T) where T
    modulus - 1 <= typemax(T) ÷ 2
end


"""
    addmod_no_overflow(x::T, y::T, modulus::T) where T <: Unsigned

Same as `addmod()`, but assumes that `x + y` doesn't overflow `T`
(that is, `modulus` is less than or equal to half of `typemax(T)+1`).
"""
@inline function addmod_no_overflow(x::T, y::T, modulus::T) where T <: Unsigned
    r = x + y
    if r < modulus
        r
    else
        r - modulus
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


@inline function doublemod(x::T, modulus::T) where T <: Unsigned
    # Could be done as `addmod(x, x)`, but this variant is faster.
    t = x
    x = double(x)
    if x < t || x >= modulus
        x = x - modulus
    end
    x
end


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

    result = zero(T)

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
