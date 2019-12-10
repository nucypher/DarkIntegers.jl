#=
Utility functions for integers.
=#

# Taken from https://github.com/JuliaLang/julia/pull/30515
# as a workaround for https://github.com/JuliaLang/julia/issues/29971
# Replace by `invmod` when it is merged.

_hastypemax(::Base.BitIntegerType) = true
_hastypemax(::Type{T}) where T = applicable(typemax, T)

function invmod_(n::T, m::T) where T <: Integer
    g, x, y = gcdx(n, m)
    g != 1 && throw(DomainError((n, m), "Greatest common divisor is $g."))
    m == 0 && throw(DomainError(m, "`m` must not be 0."))
    # Note that m might be negative here.
    r = (T <: Unsigned && _hastypemax(typeof(x)) && x > typemax(x)>>1) ? mod(x + m, m) : mod(x, m)
    # The postcondition is: mod(r * n, m) == mod(T(1), m) && div(r, m) == 0
    r
end


@inline halve(x::Integer) = x >> 1


@inline double(x::Integer) = x << 1


"""
    bitsizeof(::Type{T}) where T
    bitsizeof(x)

Size, in bits, of the canonical binary representation of the given type `T`.
Size, in bits, of object `x` if it is not a type.
"""
bitsizeof(::Type{T}) where T = sizeof(T) << 3
bitsizeof(::T) where T = bitsizeof(T)


"""
    num_bits(x)

Returns the number of bits in the representation of the absolute value of an integer.
"""
num_bits(x::T) where T <: Unsigned = bitsizeof(T) - leading_zeros(x)

function num_bits(x::T) where T <: Signed
    if x == typemin(T)
        bitsizeof(T)
    else
        num_bits(unsigned(signbit(x) ? -x : x))
    end
end

function num_bits(x::BigInt)
    if iszero(x)
        return 0
    end

    if signbit(x)
        x = -x
    end

    for i in abs(x.size):-1:1
        limb = unsafe_load(x.d, i)

        # BigInts seem to resize automatically, but we're venturing
        # into the undocumented territory here, so just in case
        # handle possible empty limbs.
        if !iszero(limb)
            limb_bits = sizeof(limb) << 3
            return (i - 1) * limb_bits + (limb_bits - leading_zeros(limb))
        end
    end
end


"""
    encompassing_type(tp::Type{<:Unsigned})

Returns the built-in type that covers all the range of `tp`
(not necessarily unsigned).
"""
encompassing_type(tp::Type{<:Unsigned}) = tp
