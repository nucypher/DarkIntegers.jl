struct _NoConversion
end

const _no_conversion = _NoConversion()


struct RRElemMontgomery{T, M} <: AbstractRRElem

    value :: T

    @inline function RRElemMontgomery(x::T, m::T, ::_NoConversion) where T <: Unsigned
        new{T, m}(x)
    end

    @inline function RRElemMontgomery(x::T, m::T) where T <: Unsigned
        RRElemMontgomery(to_montgomery(RRElemMontgomery{T, m}, x), m, _no_conversion)
    end

    @inline function RRElemMontgomery(x::RRElem{T, M}) where {T <: Unsigned, M}
        RRElemMontgomery(x.value, M)
    end

    @inline function RRElemMontgomery{T, M}(x::Integer) where {T <: Unsigned, M}
        if x < 0
            -_make_rr_elem_montgomery(unsigned(-x), M)
        else
            _make_rr_elem_montgomery(unsigned(x), M)
        end
    end

    @inline function RRElemMontgomery{T, M}(x::BigInt) where {T <: Unsigned, M}
        _make_rr_elem_montgomery(x, M)
    end
end


@inline function _make_rr_elem_montgomery(x::V, m::T) where {V <: Integer, T <: Unsigned}
    # Assumes that `x` is non-negative
    if bitsizeof(T) >= bitsizeof(V)
        RRElemMontgomery(mod(convert(T, x), m), m)
    else
        RRElemMontgomery(convert(T, mod(x, convert(V, m))), m)
    end
end


@inline @generated function montgomery_coeff(::Type{RRElemMontgomery{T, M}}) where {T, M}
    res = get_montgomery_coeff(M)
    :( $res )
end


@inline function from_montgomery(x::RRElemMontgomery{T, M}) where {T, M}
    from_montgomery(x.value, M, montgomery_coeff(RRElemMontgomery{T, M}))
end


@inline @generated function to_montgomery_coeff(::Type{RRElemMontgomery{T, M}}) where {T, M}
    res = get_to_montgomery_coeff(M)
    :( $res )
end


@inline function to_montgomery(::Type{RRElemMontgomery{T, M}}, x::T) where {T, M}
    to_montgomery(x, M, to_montgomery_coeff(RRElemMontgomery{T, M}))
end


@inline function Base.convert(::Type{RRElem{T, M}}, x::RRElemMontgomery{T, M}) where {T, M}
    RRElem(from_montgomery(x), M)
end

# TODO: this is used to prevent the convert(Integer, MPNumber) to activate.
# Is there a better way? Technically, this shouldn't be used at all - it's the constructor's job.
@inline Base.convert(::Type{RRElemMontgomery{T, M}}, x::MPNumber) where {T, M} =
    RRElemMontgomery(x, M)

@inline Base.convert(::Type{RRElemMontgomery{T, M}}, x::RRElemMontgomery{T, M}) where {T, M} = x

@inline function Base.convert(::Type{V}, x::RRElemMontgomery{T, M}) where V <: Integer where {T, M}
    convert(V, from_montgomery(x))
end


@inline function change_base_type(::Type{V}, x::RRElemMontgomery{T, M}) where {T, M, V <: Unsigned}
    RRElemMontgomery(convert(V, x.value), convert(V, M), _no_conversion)
end


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
    RRElemMontgomery(one(T), M)


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
    RRElemMontgomery(div(x_T, y_T), M)
end


function Base.div(x::RRElemMontgomery{T, M}, y::Unsigned) where {T, M}
    x_T = from_montgomery(x)
    y_T = convert(T, y) # TODO: assumes that `y` fits into RRElem
    RRElemMontgomery(div(x_T, y_T), M)
end


# Apparently we cannot just define a method for `y::Integer`, since there is a
# `div(Unsigned, Union{...})` in Base, resulting in ambiguity.
function Base.div(x::RRElemMontgomery{T, M}, y::Union{Int128, Int16, Int32, Int64, Int8}) where {T, M}
    y < 0 ? div(-x, unsigned(-y)) : div(x, unsigned(y))
end


function Base.divrem(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where {T, M}
    # TODO: optimize?
    x_T = from_montgomery(x)
    y_T = from_montgomery(y)
    d, r = divrem(x_T, y_T)
    RRElemMontgomery(d, M), RRElemMontgomery(r, M)
end


Base.string(x::RRElemMontgomery{T, M}) where {T, M} = string(x.value) * "RRM"


Base.show(io::IO, x::RRElemMontgomery{T, M}) where {T, M} = print(io, string(x))


# Required for broadcasting


Base.length(x::RRElemMontgomery{T, M}) where {T, M} = 1


Base.iterate(x::RRElemMontgomery{T, M}) where {T, M} = (x, nothing)
Base.iterate(x::RRElemMontgomery{T, M}, state) where {T, M} = nothing
