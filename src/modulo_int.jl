abstract type AbstractModUInt{T, M} <: Unsigned end


struct _Verbatim
end

const _verbatim = _Verbatim()


"""
An unsigned integer with all operations performed modulo `M`.

Supports `+`, `-`, `*`, `divrem`, `div`, `rem`, `^`, `<`, `<=`, `>`, `>=`,
`zero`, `one` and `isodd`.
Note that the division is a regular division, not multiplication by the inverse
(which is not guaranteed to exist for any `M`).

    ModUInt{T, M}(x::Integer) where {T <: Unsigned, M}

Creates an `ModUInt` object. `M` must have the type `T`.
"""
struct ModUInt{T <: Unsigned, M} <: AbstractModUInt{T, M}
    value :: T

    # This is the only method using `new` to ensure `M` has the type `T`
    # (since we cannot enforce it with Julia syntax)
    @inline function ModUInt(x::T, m::T, ::_Verbatim) where T
        new{T, m}(x)
    end

    @inline function ModUInt{T, M}(x::T, ::_Verbatim) where {T, M}
        ModUInt(x, M, _verbatim)
    end

    @inline function ModUInt{T, M}(x::Integer) where {T, M}
        ModUInt(convert(T, mod(x, M)), M, _verbatim)
    end
end


@inline Base.convert(::Type{ModUInt{T, M}}, x::T) where {T <: MLUInt, M} = ModUInt{T, M}(x)
@inline Base.convert(::Type{ModUInt{T, M}}, x::ModUInt{T, M}) where {T, M} = x
@inline Base.convert(::Type{ModUInt{T, M}}, x::Bool) where {T, M} =
    x ? one(ModUInt{T, M}) : zero(ModUInt{T, M})


@inline Base.promote_type(::Type{ModUInt{T, M}}, ::Type{ModUInt{T, M}}) where {T, M} = ModUInt{T, M}
@inline Base.promote_type(::Type{ModUInt{T, M}}, ::Type{<:Integer}) where {T, M} = ModUInt{T, M}
@inline Base.promote_type(::Type{<:Integer}, ::Type{ModUInt{T, M}}) where {T, M} = ModUInt{T, M}


# We need this to correctly process arithmetic operations on ModUInt and Int
# (which is signed and the default in Julia for number literals)
# without defining specific methods for each operator.
@inline Base.signed(x::ModUInt{T, M}) where {T, M} = x
@inline Base.unsigned(x::ModUInt{T, M}) where {T, M} = x


# Unlike `one(x)`, `zero(x)` does not have a fallback `= zero(typeof(x))` in the standard library
# and uses conversion instead. So we are defining our own.
@inline Base.zero(::Type{ModUInt{T, M}}) where {T, M} = ModUInt(zero(T), M, _verbatim)
@inline Base.zero(::ModUInt{T, M}) where {T, M} = zero(ModUInt{T, M})


@inline Base.one(::Type{ModUInt{T, M}}) where {T, M} = ModUInt(one(T), M, _verbatim)
@inline Base.one(::ModUInt{T, M}) where {T, M} = one(ModUInt{T, M})

@inline Base.oneunit(::Type{ModUInt{T, M}}) where {T, M} = one(ModUInt{T, M})
@inline Base.oneunit(x::ModUInt{T, M}) where {T, M} = oneunit(ModUInt{T, M})


@inline @generated function Base.:+(x::ModUInt{T, M}, y::ModUInt{T, M}) where {T, M}
    if addition_cant_overflow(M)
        addmod = :addmod_no_overflow
    else
        addmod = :addmod
    end

    :( ModUInt($addmod(x.value, y.value, M), M, _verbatim) )
end


@inline function Base.:-(x::ModUInt{T, M}, y::ModUInt{T, M}) where {T, M}
    ModUInt(submod(x.value, y.value, M), M, _verbatim)
end


@inline function Base.:-(x::ModUInt{T, M}) where {T, M}
    # TODO: (issue #20) can be optimized
    zero(ModUInt{T, M}) - x
end


@inline function Base.:*(x::ModUInt{T, M}, y::ModUInt{T, M}) where {T, M}
    xt = x.value
    yt = y.value
    res = mulmod(xt, yt, M)
    ModUInt(res, M, _verbatim)
end


@inline function Base.:>>(x::ModUInt{T, M}, shift::Int) where {T, M}
    if shift >= 0
        ModUInt(x.value >> shift, M, _verbatim)
    else
        ModUInt{T, M}(x.value >> shift)
    end
end


@inline function Base.:<<(x::ModUInt{T, M}, shift::Int) where {T, M}
    if shift >= 0
        ModUInt{T, M}(x.value << shift)
    else
        ModUInt(x.value << shift, M, _verbatim)
    end
end


@inline function Base.:&(x::ModUInt{T, M}, y::ModUInt{T, M}) where {T, M}
    # `&` of two numbers in [0, M) is guaranteed to produce the result in the same range,
    # so no need for reduction.
    ModUInt(x.value & y.value, M, _verbatim)
end


@inline function Base.isodd(x::ModUInt{T, M}) where {T, M}
    isodd(x.value)
end


@inline function Base.div(x::ModUInt{T, M}, y::ModUInt{T, M}) where {T, M}
    ModUInt(div(x.value, y.value), M, _verbatim)
end


@inline function Base.div(x::ModUInt{T, M}, y::Unsigned) where {T, M}
    if y >= M
        zero(ModUInt{T, M})
    else
        div(x, ModUInt{T, M}(y))
    end
end


# Apparently we cannot just define a method for `y::Integer`, since there is a
# `div(Unsigned, Union{...})` in Base, resulting in ambiguity.
@inline function Base.div(x::ModUInt{T, M}, y::Union{Int128, Int16, Int32, Int64, Int8}) where {T, M}
    signbit(y) ? div(-x, unsigned(-y)) : div(x, unsigned(y))
end


@inline function Base.divrem(x::ModUInt{T, M}, y::ModUInt{T, M}) where {T, M}
    d, r = divrem(x.value, y.value)
    ModUInt(d, M, _verbatim), ModUInt(r, M, _verbatim)
end


@inline Base.rem(x::ModUInt{T, M}, ::Type{ModUInt{T, M}}) where {T, M} = x

@inline Base.rem(x::ModUInt{T, M}, tp::Type{V}) where {T, M, V <: Integer} = rem(x.value, tp)


@inline function Base.inv(x::ModUInt{T, M}) where {T, M}
    ModUInt{T, M}(invmod_(x.value, M), _verbatim)
end


@inline Base.:<(x::ModUInt{T, M}, y::ModUInt{T, M}) where {T, M} = x.value < y.value


@inline Base.:>(x::ModUInt{T, M}, y::ModUInt{T, M}) where {T, M} = x.value > y.value


@inline Base.:<=(x::ModUInt{T, M}, y::ModUInt{T, M}) where {T, M} = x.value <= y.value


@inline Base.:>=(x::ModUInt{T, M}, y::ModUInt{T, M}) where {T, M} = x.value >= y.value


Base.string(x::ModUInt{T, M}) where {T, M} = string(value(x)) * "RR"


Base.show(io::IO, x::ModUInt{T, M}) where {T, M} = print(io, string(x))


Base.Broadcast.broadcastable(x::ModUInt) = (x,)


encompassing_type(tp::Type{ModUInt{T, M}}) where {T, M} = encompassing_type(T)


# Random number generation


function Base.rand(rng::AbstractRNG, ::Random.SamplerType{ModUInt{T, M}}) where {T, M}
    ModUInt{T, M}(rand(rng, zero(T):M-one(T)), _verbatim)
end


struct ModUIntSampler{Z, S} <: Random.Sampler{Z}
    sampler :: S

    function ModUIntSampler(RNG::Type{<:AbstractRNG}, r::UnitRange{ModUInt{T, M}}, n) where {T, M}
        sampler = Random.Sampler(RNG, value(r.start):value(r.stop), n)
        new{ModUInt{T, M}, typeof(sampler)}(sampler)
    end
end


function Random.Sampler(
        RNG::Type{<:AbstractRNG}, r::UnitRange{ModUInt{T, M}}, n::Union{Val{1}, Val{Inf}}) where {T, M}
    ModUIntSampler(RNG, r, n)
end


function Base.rand(rng::AbstractRNG, s::ModUIntSampler{Z, S}) where {Z, S}
    Z(rand(rng, s.sampler), _verbatim)
end


function num_bits(x::ModUInt{T, M}) where {T, M}
    num_bits(value(x))
end


function Base.trailing_zeros(x::ModUInt{T, M}) where {T, M}
    trailing_zeros(value(x))
end
