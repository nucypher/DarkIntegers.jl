"""
Residue ring element in Montgomery representation.
Multiplication is much faster than for [`RRElem`](@ref), addition and subtraction are the same,
but division and conversion to regular integers is slower.

Supports `+`, `-`, `*`, `divrem`, `div`, `rem`, `^`, `<`, `<=`, `>`, `>=`,
`zero`, `one` and `isodd`.
Note that the division is a regular division, not multiplication by the inverse
(which is not guaranteed to exist for any `M`).

    RRElemMontgomery{T, M}(x::Integer) where {T <: Unsigned, M}

Creates an `RRElemMontgomery` object. `M` must have the type `T` and be an odd number.
"""
struct RRElemMontgomery{T, M} <: AbstractRRElem{T, M}

    value :: T

    # This is the only method using `new` to ensure `M` has the type `T`
    # (since we cannot enforce it with Julia syntax)
    @inline function RRElemMontgomery(x::T, m::T, ::_NoConversion) where T <: Unsigned
        new{T, m}(x)
    end

    @inline function RRElemMontgomery{T, M}(x::T, ::_NoConversion) where {T <: Unsigned, M}
        RRElemMontgomery(x, M, _no_conversion)
    end

    @inline function RRElemMontgomery{T, M}(x::Integer) where {T <: Unsigned, M}
        # Need to take the modulus before converting `x` to `T`,
        # in case `x` does not fit in `T`.
        RRElemMontgomery{T, M}(
            to_montgomery(RRElemMontgomery{T, M}, convert(T, mod(x, M))), _no_conversion)
    end
end


"""
Caches the coefficient used for multiplication and conversion from Montgomery representation
based on the type.
"""
@inline @generated function montgomery_coeff(::Type{RRElemMontgomery{T, M}}) where {T, M}
    res = get_montgomery_coeff(M)
    :( $res )
end


@inline function from_montgomery(x::RRElemMontgomery{T, M}) where {T, M}
    from_montgomery(x.value, M, montgomery_coeff(RRElemMontgomery{T, M}))
end


"""
Caches the coefficient used for conversion to Montgomery representation based on the type.
"""
@inline @generated function to_montgomery_coeff(::Type{RRElemMontgomery{T, M}}) where {T, M}
    res = get_to_montgomery_coeff(M)
    :( $res )
end


@inline function to_montgomery(::Type{RRElemMontgomery{T, M}}, x::T) where {T, M}
    to_montgomery(x, M, to_montgomery_coeff(RRElemMontgomery{T, M}))
end


# Needed in cases when convert(RRElemMontgomery, MPNumber) is requested
# (including implicitly, e.g. in array element assignment)
# to prevent convert(::Type{<:Integer}, x::MPNumber) from firing.
@inline Base.convert(::Type{RRElemMontgomery{T, M}}, x::V) where {T, M, V <: MPNumber} =
    RRElemMontgomery{T, M}(convert(encompassing_type(V), x))
@inline Base.convert(::Type{RRElemMontgomery{T, M}}, x::T) where {T <: MPNumber, M} =
    RRElemMontgomery{T, M}(x)

@inline function Base.convert(::Type{RRElem{T, M}}, x::RRElemMontgomery{T, M}) where {T, M}
    RRElem(from_montgomery(x), M, _no_conversion)
end

@inline function Base.convert(::Type{RRElemMontgomery{T, M}}, x::RRElem{T, M}) where {T, M}
    RRElemMontgomery{T, M}(x.value)
end

@inline Base.convert(::Type{RRElemMontgomery{T, M}}, x::RRElemMontgomery{T, M}) where {T, M} = x

@inline function Base.convert(::Type{V}, x::RRElemMontgomery{T, M}) where V <: Integer where {T, M}
    convert(V, from_montgomery(x))
end


@inline Base.promote_type(
    ::Type{RRElemMontgomery{T, M}}, ::Type{RRElemMontgomery{T, M}}) where {T, M} =
    RRElemMontgomery{T, M}
@inline Base.promote_type(::Type{RRElemMontgomery{T, M}}, ::Type{<:Integer}) where {T, M} =
    RRElemMontgomery{T, M}
@inline Base.promote_type(::Type{<:Integer}, ::Type{RRElemMontgomery{T, M}}) where {T, M} =
    RRElemMontgomery{T, M}


# We need this to correctly process arithmetic operations on RRElemMontgomery and Int
# (which is signed and the default in Julia for number literals)
# without defining specific methods for each operator.
@inline Base.signed(x::RRElemMontgomery{T, M}) where {T, M} = x
@inline Base.unsigned(x::RRElemMontgomery{T, M}) where {T, M} = x


# Unlike `one(x)`, `zero(x)` does not have a fallback `= zero(typeof(x))` in the standard library
# and uses conversion instead. So we are defining our own.
@inline Base.zero(::Type{RRElemMontgomery{T, M}}) where {T, M} =
    RRElemMontgomery(zero(T), M, _no_conversion)
@inline Base.zero(::RRElemMontgomery{T, M}) where {T, M} =
    zero(RRElemMontgomery{T, M})


@inline Base.one(::Type{RRElemMontgomery{T, M}}) where {T, M} =
    RRElemMontgomery{T, M}(one(T))


@inline function Base.:+(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where {T, M}
    RRElemMontgomery(addmod(x.value, y.value, M), M, _no_conversion)
end


@inline function Base.:-(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where {T, M}
    RRElemMontgomery(submod(x.value, y.value, M), M, _no_conversion)
end


@inline function Base.:-(x::RRElemMontgomery{T, M}) where {T, M}
    # TODO: can be optimized
    zero(RRElemMontgomery{T, M}) - x
end


@inline function Base.:*(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where {T, M}
    res = mulmod_montgomery(x.value, y.value, M, montgomery_coeff(RRElemMontgomery{T, M}))
    RRElemMontgomery(res, M, _no_conversion)
end


function Base.isodd(x::RRElemMontgomery{T, M}) where {T, M}
    # TODO: optimize? Although currently it is not critical to the performance
    isodd(from_montgomery(x))
end


function Base.div(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where {T, M}
    x_T = from_montgomery(x)
    y_T = from_montgomery(y)
    RRElemMontgomery{T, M}(div(x_T, y_T))
end


function Base.div(x::RRElemMontgomery{T, M}, y::Unsigned) where {T, M}
    if y >= M
        zero(RRElemMontgomery{T, M})
    else
        x_T = from_montgomery(x)
        y_T = convert(T, y)
        RRElemMontgomery{T, M}(div(x_T, y_T))
    end
end


# Apparently we cannot just define a method for `y::Integer`, since there is a
# `div(Unsigned, Union{...})` in Base, resulting in ambiguity.
function Base.div(x::RRElemMontgomery{T, M}, y::Union{Int128, Int16, Int32, Int64, Int8}) where {T, M}
    y < 0 ? div(-x, unsigned(-y)) : div(x, unsigned(y))
end


function Base.divrem(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where {T, M}
    x_T = from_montgomery(x)
    y_T = from_montgomery(y)
    d, r = divrem(x_T, y_T)
    RRElemMontgomery{T, M}(d), RRElemMontgomery{T, M}(r)
end


Base.string(x::RRElemMontgomery{T, M}) where {T, M} =
    string(rr_value(change_representation(RRElem, x))) * "RRM"


Base.show(io::IO, x::RRElemMontgomery{T, M}) where {T, M} = print(io, string(x))


# Required for broadcasting


Base.length(x::RRElemMontgomery{T, M}) where {T, M} = 1


Base.iterate(x::RRElemMontgomery{T, M}) where {T, M} = (x, nothing)
Base.iterate(x::RRElemMontgomery{T, M}, state) where {T, M} = nothing
