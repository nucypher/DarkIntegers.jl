abstract type AbstractRRElem{T, M} <: Unsigned end


struct _Verbatim
end

const _verbatim = _Verbatim()


"""
Residue ring element, an unsigned integer with all operations performed modulo `M`.

Supports `+`, `-`, `*`, `divrem`, `div`, `rem`, `^`, `<`, `<=`, `>`, `>=`,
`zero`, `one` and `isodd`.
Note that the division is a regular division, not multiplication by the inverse
(which is not guaranteed to exist for any `M`).

    RRElem{T, M}(x::Integer) where {T <: Unsigned, M}

Creates an `RRElem` object. `M` must have the type `T`.
"""
struct RRElem{T, M} <: AbstractRRElem{T, M}
    value :: T

    # This is the only method using `new` to ensure `M` has the type `T`
    # (since we cannot enforce it with Julia syntax)
    @inline function RRElem(x::T, m::T, ::_Verbatim) where T <: Unsigned
        new{T, m}(x)
    end

    @inline function RRElem{T, M}(x::T, ::_Verbatim) where {T <: Unsigned, M}
        new{T, M}(x)
    end

    @inline function RRElem{T, M}(x::Integer) where {T <: Unsigned, M}
        RRElem{T, M}(convert(T, mod(x, M)), _verbatim)
    end
end

# Needed in cases when convert(RRElem, MPNumber) is requested
# (including implicitly, e.g. in array element assignment)
# to prevent convert(::Type{<:Integer}, x::MPNumber) from firing.
@inline Base.convert(::Type{RRElem{T, M}}, x::V) where {T, M, V <: MPNumber} =
    RRElem{T, M}(convert(encompassing_type(V), x))
@inline Base.convert(::Type{RRElem{T, M}}, x::T) where {T <: MPNumber, M} = RRElem{T, M}(x)

@inline Base.convert(::Type{RRElem{T, M}}, x::RRElem{T, M}) where {T, M} = x
@inline Base.convert(::Type{V}, x::RRElem{T, M}) where {V <: Integer, T, M} = convert(V, x.value)


@inline Base.promote_type(::Type{RRElem{T, M}}, ::Type{RRElem{T, M}}) where {T, M} = RRElem{T, M}
@inline Base.promote_type(::Type{RRElem{T, M}}, ::Type{<:Integer}) where {T, M} = RRElem{T, M}
@inline Base.promote_type(::Type{<:Integer}, ::Type{RRElem{T, M}}) where {T, M} = RRElem{T, M}


# We need this to correctly process arithmetic operations on RRElem and Int
# (which is signed and the default in Julia for number literals)
# without defining specific methods for each operator.
@inline Base.signed(x::RRElem{T, M}) where {T, M} = x
@inline Base.unsigned(x::RRElem{T, M}) where {T, M} = x


# Unlike `one(x)`, `zero(x)` does not have a fallback `= zero(typeof(x))` in the standard library
# and uses conversion instead. So we are defining our own.
@inline Base.zero(::Type{RRElem{T, M}}) where {T, M} = RRElem(zero(T), M, _verbatim)
@inline Base.zero(::RRElem{T, M}) where {T, M} = zero(RRElem{T, M})


@inline Base.one(::Type{RRElem{T, M}}) where {T, M} = RRElem(one(T), M, _verbatim)


@inline @generated function Base.:+(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M}
    if addition_cant_overflow(M)
        addmod = :addmod_no_overflow
    else
        addmod = :addmod
    end

    :( RRElem($addmod(x.value, y.value, M), M, _verbatim) )
end


@inline function Base.:-(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M}
    RRElem(submod(x.value, y.value, M), M, _verbatim)
end


@inline function Base.:-(x::RRElem{T, M}) where {T, M}
    # TODO: can be optimized
    zero(RRElem{T, M}) - x
end


@inline function Base.:*(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M}
    xt = x.value
    yt = y.value
    res = mulmod(xt, yt, M)
    RRElem(res, M, _verbatim)
end


@inline function Base.isodd(x::RRElem{T, M}) where {T, M}
    isodd(x.value)
end


@inline function Base.div(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M}
    RRElem(div(x.value, y.value), M, _verbatim)
end


@inline function Base.div(x::RRElem{T, M}, y::Unsigned) where {T, M}
    if y >= M
        zero(RRElem{T, M})
    else
        div(x, RRElem{T, M}(y))
    end
end


# Apparently we cannot just define a method for `y::Integer`, since there is a
# `div(Unsigned, Union{...})` in Base, resulting in ambiguity.
@inline function Base.div(x::RRElem{T, M}, y::Union{Int128, Int16, Int32, Int64, Int8}) where {T, M}
    y < 0 ? div(-x, unsigned(-y)) : div(x, unsigned(y))
end


@inline function Base.divrem(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M}
    d, r = divrem(x.value, y.value)
    RRElem(d, M, _verbatim), RRElem(r, M, _verbatim)
end


@inline function Base.inv(x::RRElem{T, M}) where {T, M}
    RRElem{T, M}(invmod_(x.value, M), _verbatim)
end


@inline Base.:<(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M} = x.value < y.value


@inline Base.:>(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M} = x.value > y.value


@inline Base.:<=(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M} = x.value <= y.value


@inline Base.:>=(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M} = x.value >= y.value


Base.string(x::RRElem{T, M}) where {T, M} = string(x.value) * "RR"


Base.show(io::IO, x::RRElem{T, M}) where {T, M} = print(io, string(x))


# Required for broadcasting


Base.length(x::RRElem{T, M}) where {T, M} = 1


Base.iterate(x::RRElem{T, M}) where {T, M} = (x, nothing)
Base.iterate(x::RRElem{T, M}, state) where {T, M} = nothing
