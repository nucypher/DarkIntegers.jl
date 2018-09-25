abstract type AbstractRRElem <: Unsigned end


struct RRElem{T, M} <: AbstractRRElem
    value :: T

    # TODO: or is it better to have `function RRElem{T, M}(x::T)`?
    @inline function RRElem(x::T, m::T) where T
        new{T, m}(x)
    end

    # If we're converting a negative number to a RRElem,
    # a conversion Integer -> MPNumber -> RRElem will give a wrong result.
    # So it has to be performed directly.
    @inline function RRElem{T, M}(x::Integer) where {T, M}
        # TODO: write an optimized version of this
        # TODO: if `x` is already a MPNumber, no need to convert to BigInt
        # @assert typeof(M) == T
        m_i = convert(BigInt, M)
        x_i = convert(BigInt, x)
        new{T, M}(mod(x_i, m_i))
    end
end


# TODO: this is used to prevent the convert(Integer, MPNumber) to activate.
# Is there a better way? Technically, this shouldn't be used at all - it's the constructor's job.
@inline Base.convert(::Type{RRElem{T, M}}, x::MPNumber) where {T, M} = RRElem(x, M)
@inline Base.convert(::Type{RRElem{T, M}}, x::RRElem{T, M}) where {T, M} = x
@inline Base.convert(::Type{<:Integer}, x::RRElem{T, M}) where {T, M} = convert(V, x.value)


@inline Base.promote_type(::Type{RRElem{T, M}}, ::Type{<:Integer}) where {T, M} = RRElem{T, M}
@inline Base.promote_type(::Type{<:Integer}, ::Type{RRElem{T, M}}) where {T, M} = RRElem{T, M}


# We need this to correctly process arithmetic operations on RRElem and Int
# (which is signed and the default in Julia for number literals)
# without defining specific methods for each operator.
@inline Base.signed(x::RRElem{T, M}) where {T, M} = x
@inline Base.unsigned(x::RRElem{T, M}) where {T, M} = x


@inline Base.zero(::Type{RRElem{T, M}}) where {T, M} = RRElem(zero(T), M)


@inline Base.one(::Type{RRElem{T, M}}) where {T, M} = RRElem(one(T), M)


@inline function Base.:+(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M}
    RRElem(addmod(x.value, y.value, M), M)
end


@inline function Base.:-(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M}
    RRElem(submod(x.value, y.value, M), M)
end


@inline function Base.:-(x::RRElem{T, M}) where {T, M}
    # TODO: can be optimized
    zero(RRElem{T, M}) - x
end


@inline function Base.:*(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M}
    xt = x.value
    yt = y.value
    res = mulmod_widemul(xt, yt, M)
    RRElem(res, M)
end


@inline function with_modulus(x::RRElem{T, M}, new_modulus::Unsigned) where {T, M}
    nm = convert(T, new_modulus)
    convert(RRElem{T, nm}, x)
end


@inline function modulus_reduction(x::RRElem{T, M}, new_modulus::Unsigned) where {T, M}
    nm = convert(T, new_modulus)
    # TODO: optimize
    xi = convert(BigInt, x)
    mi = convert(BigInt, M)

    # TODO: make it a purely integer algorithm
    convert(RRElem{T, nm}, round(BigInt, xi * new_modulus / mi))
end


@inline function Base.convert(::Type{RRElem{T, N}}, x::RRElem{T, M}) where {T, N, M}
    if N >= M
        RRElem(x.value, N)
    else
        # TODO: optimize
        RRElem(convert(T, convert(BigInt, x.value) % N))
    end
end


@inline function Base.isodd(x::RRElem{T, M}) where {T, M}
    isodd(x.value)
end


@inline function Base.div(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M}
    RRElem(div(x.value, y.value), M)
end


@inline function Base.div(x::RRElem{T, M}, y::Unsigned) where {T, M}
    # TODO: assumes that `y` fits into RRElem
    div(x, convert(RRElem{T, M}, y))
end


# Apparently we cannot just define a method for `y::Integer`, since there is a
# `div(Unsigned, Union{...})` in Base, resulting in ambiguity.
@inline function Base.div(x::RRElem{T, M}, y::Union{Int128, Int16, Int32, Int64, Int8}) where {T, M}
    y < 0 ? div(-x, unsigned(-y)) : div(x, unsigned(y))
end


@inline function Base.divrem(x::RRElem{T, M}, y::RRElem{T, M}) where {T, M}
    d, r = divrem(x.value, y.value)
    RRElem(d, M), RRElem(r, M)
end


@inline function Base.:^(x::T, y::Unsigned) where T <: AbstractRRElem
    # TODO: optimize
    res = one(T)
    @assert y >= 0
    for i in 1:y
        res *= x
    end
    res
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
