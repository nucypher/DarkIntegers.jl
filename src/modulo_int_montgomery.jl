struct _NoModulo
end

const _no_modulo = _NoModulo()


"""
Residue ring element in Montgomery representation.
Multiplication is much faster than for [`ModUInt`](@ref), addition and subtraction are the same,
but division and conversion to regular integers is slower.

Supports `+`, `-`, `*`, `divrem`, `div`, `rem`, `^`, `<`, `<=`, `>`, `>=`,
`zero`, `one` and `isodd`.
Note that the division is a regular division, not multiplication by the inverse
(which is not guaranteed to exist for any `M`).

    MgModUInt{T, M}(x::Integer) where {T <: Unsigned, M}

Creates an `MgModUInt` object. `M` must have the type `T` and be an odd number.
"""
struct MgModUInt{T, M} <: AbstractModUInt{T, M}

    value :: T

    # This is the only method using `new` to ensure `M` has the type `T`
    # (since we cannot enforce it with Julia syntax)
    @inline function MgModUInt(x::T, m::T, ::_Verbatim) where T <: Unsigned
        new{T, m}(x)
    end

    @inline function MgModUInt{T, M}(x::T, ::_Verbatim) where {T <: Unsigned, M}
        new{T, M}(x)
    end

    @inline function MgModUInt{T, M}(x::T, ::_NoModulo) where {T <: Unsigned, M}
        MgModUInt{T, M}(to_montgomery(MgModUInt{T, M}, x), _verbatim)
    end

    @inline function MgModUInt{T, M}(x::Integer) where {T <: Unsigned, M}
        # Need to take the modulus before converting `x` to `T`,
        # in case `x` does not fit in `T`.
        MgModUInt{T, M}(convert(T, mod(x, M)), _no_modulo)
    end
end


"""
Caches the coefficient used for multiplication and conversion from Montgomery representation
based on the type.
"""
@inline @generated function montgomery_coeff(::Type{MgModUInt{T, M}}) where {T, M}
    res = get_montgomery_coeff(M)
    :( $res )
end


@inline function from_montgomery(x::MgModUInt{T, M}) where {T, M}
    from_montgomery(x.value, M, montgomery_coeff(MgModUInt{T, M}))
end


"""
Caches the coefficient used for conversion to Montgomery representation based on the type.
"""
@inline @generated function to_montgomery_coeff(::Type{MgModUInt{T, M}}) where {T, M}
    res = get_to_montgomery_coeff(M)
    :( $res )
end


@inline function to_montgomery(::Type{MgModUInt{T, M}}, x::T) where {T, M}
    m_coeff = montgomery_coeff(MgModUInt{T, M})
    to_m_coeff = to_montgomery_coeff(MgModUInt{T, M})
    to_montgomery(x, M, m_coeff, to_m_coeff)
end


# Needed in cases when convert(MgModUInt, MLUInt) is requested
# (including implicitly, e.g. in array element assignment)
# to prevent convert(::Type{<:Integer}, x::MLUInt) from firing.
@inline Base.convert(::Type{MgModUInt{T, M}}, x::V) where {T, M, V <: MLUInt} =
    MgModUInt{T, M}(convert(encompassing_type(V), x))
@inline Base.convert(::Type{MgModUInt{T, M}}, x::T) where {T <: MLUInt, M} =
    MgModUInt{T, M}(x)

@inline function Base.convert(::Type{ModUInt{T, M}}, x::MgModUInt{T, M}) where {T, M}
    ModUInt(from_montgomery(x), M, _verbatim)
end

@inline function Base.convert(::Type{MgModUInt{T, M}}, x::ModUInt{T, M}) where {T, M}
    MgModUInt{T, M}(x.value, M, _no_modulo)
end

@inline Base.convert(::Type{MgModUInt{T, M}}, x::MgModUInt{T, M}) where {T, M} = x

@inline function Base.convert(::Type{V}, x::MgModUInt{T, M}) where V <: Integer where {T, M}
    convert(V, from_montgomery(x))
end


@inline Base.promote_type(
    ::Type{MgModUInt{T, M}}, ::Type{MgModUInt{T, M}}) where {T, M} =
    MgModUInt{T, M}
@inline Base.promote_type(::Type{MgModUInt{T, M}}, ::Type{<:Integer}) where {T, M} =
    MgModUInt{T, M}
@inline Base.promote_type(::Type{<:Integer}, ::Type{MgModUInt{T, M}}) where {T, M} =
    MgModUInt{T, M}


# We need this to correctly process arithmetic operations on MgModUInt and Int
# (which is signed and the default in Julia for number literals)
# without defining specific methods for each operator.
@inline Base.signed(x::MgModUInt{T, M}) where {T, M} = x
@inline Base.unsigned(x::MgModUInt{T, M}) where {T, M} = x


# Unlike `one(x)`, `zero(x)` does not have a fallback `= zero(typeof(x))` in the standard library
# and uses conversion instead. So we are defining our own.
@inline Base.zero(::Type{MgModUInt{T, M}}) where {T, M} =
    MgModUInt(zero(T), M, _verbatim)
@inline Base.zero(::MgModUInt{T, M}) where {T, M} =
    zero(MgModUInt{T, M})


@inline Base.one(::Type{MgModUInt{T, M}}) where {T, M} =
    MgModUInt{T, M}(one(T))


@inline @generated function Base.:+(x::MgModUInt{T, M}, y::MgModUInt{T, M}) where {T, M}
    if addition_cant_overflow(M)
        addmod = :addmod_no_overflow
    else
        addmod = :addmod
    end
    :( MgModUInt($addmod(x.value, y.value, M), M, _verbatim) )
end


@inline function Base.:-(x::MgModUInt{T, M}, y::MgModUInt{T, M}) where {T, M}
    MgModUInt(submod(x.value, y.value, M), M, _verbatim)
end


@inline function Base.:-(x::MgModUInt{T, M}) where {T, M}
    # TODO: can be optimized
    zero(MgModUInt{T, M}) - x
end


@inline function Base.:*(x::MgModUInt{T, M}, y::MgModUInt{T, M}) where {T, M}
    res = mulmod_montgomery(x.value, y.value, M, montgomery_coeff(MgModUInt{T, M}))
    MgModUInt{T, M}(res, _verbatim)
end

#=
Montgomery multiplication gives `x, y -> x * y * R^(-1) mod M`.
So it works for one ModUInt and MgModUInt, returning an ModUInt:
`x, y * R -> x * y * R * R^(-1) mod M == x * y mod M`.
=#
@inline function Base.:*(x::ModUInt{T, M}, y::MgModUInt{T, M}) where {T, M}
    res = mulmod_montgomery(x.value, y.value, M, montgomery_coeff(MgModUInt{T, M}))
    ModUInt{T, M}(res, _verbatim)
end

@inline function Base.:*(x::MgModUInt{T, M}, y::ModUInt{T, M}) where {T, M}
    y * x
end


function Base.isodd(x::MgModUInt{T, M}) where {T, M}
    # TODO: optimize? Although currently it is not critical to the performance
    isodd(from_montgomery(x))
end


function Base.div(x::MgModUInt{T, M}, y::MgModUInt{T, M}) where {T, M}
    x_T = from_montgomery(x)
    y_T = from_montgomery(y)
    MgModUInt{T, M}(div(x_T, y_T))
end


function Base.div(x::MgModUInt{T, M}, y::Unsigned) where {T, M}
    if y >= M
        zero(MgModUInt{T, M})
    else
        x_T = from_montgomery(x)
        y_T = convert(T, y)
        MgModUInt{T, M}(div(x_T, y_T))
    end
end


# Apparently we cannot just define a method for `y::Integer`, since there is a
# `div(Unsigned, Union{...})` in Base, resulting in ambiguity.
function Base.div(x::MgModUInt{T, M}, y::Union{Int128, Int16, Int32, Int64, Int8}) where {T, M}
    y < 0 ? div(-x, unsigned(-y)) : div(x, unsigned(y))
end


function Base.divrem(x::MgModUInt{T, M}, y::MgModUInt{T, M}) where {T, M}
    x_T = from_montgomery(x)
    y_T = from_montgomery(y)
    d, r = divrem(x_T, y_T)
    MgModUInt{T, M}(d), MgModUInt{T, M}(r)
end


@inline function Base.inv(x::MgModUInt{T, M}) where {T, M}
    value = from_montgomery(x)
    MgModUInt{T, M}(invmod_(value, M), _no_modulo)
end


Base.string(x::MgModUInt{T, M}) where {T, M} =
    string(rr_value(change_representation(ModUInt, x))) * "RRM"


Base.show(io::IO, x::MgModUInt{T, M}) where {T, M} = print(io, string(x))


# Required for broadcasting


Base.length(x::MgModUInt{T, M}) where {T, M} = 1


Base.iterate(x::MgModUInt{T, M}) where {T, M} = (x, nothing)
Base.iterate(x::MgModUInt{T, M}, state) where {T, M} = nothing


encompassing_type(tp::Type{MgModUInt{T, M}}) where {T, M} = encompassing_type(T)
