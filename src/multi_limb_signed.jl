#=
Multi-limb signed integers.
=#

struct MLInt{N, T <: Unsigned} <: Signed
    limbs :: NTuple{N, T}

    MLInt(x::NTuple{N, T}) where {N, T} = new{N, T}(x)
    MLInt{N, T}(x::MLInt{N, T}) where {N, T} = new{N, T}(x.limbs)
    MLInt{N, T}(x::NTuple{N, V}) where {N, T, V} = new{N, T}(convert.(T, x))
end


@inline Base.signed(x::MLInt{N, T}) where {N, T} = x


@inline Base.unsigned(x::MLInt{N, T}) where {N, T} = MLUInt{N, T}(x.limbs)


@inline function _most_significant_limb(x::MLInt{N, T}) where {N, T}
    s = signbit(x)
    for i in N:-1:1
        if (!s && !iszero(x[i])) || (s && x[i] != typemax(T))
            return i
        end
    end
    return 0
end


@inline function _unsafe_convert(::Type{MLInt{N, T}}, x::MLInt{M, T}) where {N, M, T}
    res = signed(_unsafe_convert(MLUInt{N, T}, unsigned(x)))
    if signbit(x)
        for i in min(N, M)+1:N
            res = setindex(res, typemax(T), i)
        end
    end
    res
end


@inline function _unsafe_convert(::Type{V}, x::MLInt{N, T}) where {V <: Integer, N, T}
    res = signed(_unsafe_convert(V, unsigned(x)))
    if signbit(x) && (V == BigInt || bitsizeof(MLInt{N, T}) < bitsizeof(V))
        res -= one(V) << bitsizeof(MLInt{N, T})
    end
    res
end


@inline function _unsafe_convert(::Type{MLInt{N, T}}, x::V) where {N, T, V <: Integer}
    if x == typemin(V)
        signed(typemax(MLUInt{N, T}) << (bitsizeof(V) - 1))
    else
        res = signed(_unsafe_convert(MLUInt{N, T}, unsigned(abs(x))))
        signbit(x) ? -res : res
    end
end


@inline function _unsafe_convert(::Type{MLInt{N, T}}, x::BigInt) where {N, T}
    res = signed(_unsafe_convert(MLUInt{N, T}, abs(x)))
    signbit(x) ? -res : res
end


@inline Base.convert(::Type{MLInt{N, T}}, x::MLInt{N, T}) where {N, T} = x

@inline function Base.convert(::Type{MLInt{N, T}}, x::MLInt{M, T}) where {N, M, T}
    if _most_significant_limb(x) > N
        throw(InexactError(:convert, MLInt{N, T}, x))
    end
    _unsafe_convert(MLInt{N, T}, x)
end

@inline function Base.convert(::Type{V}, x::MLInt{N, T}) where {V <: Signed, N, T}
    if V != BigInt
        nb = num_bits(x)
        bsv = bitsizeof(V)
        bsx = bitsizeof(x)

        # nb == bsx for MLInt can happen only if `x == typemin(MLInt{N, T})`
        is_typemin = nb == bsx

        # In case of `x == typemin(MLInt{N, T})` we just need num_bits(x) to fit into V
        # (because `nb` already includes the sign bit)
        if is_typemin && nb > bsv
            throw(InexactError(:convert, V, x))
        # if x != typemin() we need num_bits(x) + 1 bit for the sign to fit into V
        elseif !is_typemin && nb >= bsv
            throw(InexactError(:convert, V, x))
        end
    end
    _unsafe_convert(V, x)
end

@inline function Base.convert(::Type{V}, x::MLInt{N, T}) where {V <: Unsigned, N, T}
    # Similar behaviour to built-in Julia types - can't convert negative numbers to unsigned.
    if signbit(x)
        throw(InexactError(:convert, V, x))
    # Check if `x` fits in the target type
    elseif num_bits(x) >= bitsizeof(V)
        throw(InexactError(:convert, V, x))
    end
    _unsafe_convert(V, x)
end

@inline function Base.convert(::Type{MLInt{N, T}}, x::Integer) where {N, T}
    bs = bitsizeof(MLInt{N, T})
    nb = num_bits(x)
    s = signbit(x)

    # if `x == typemin(MLInt{N, T})` (need to check it without actual conversion, naturally)
    is_typemin = s && nb == bs && trailing_zeros(x) == bs - 1

    # If `x == typemin(MLInt{N, T})`, obviously it fits.
    # Otherwise we need num_bits(x) + 1 bit for the sign to fit into the target
    if !is_typemin && nb >= bs
        throw(InexactError(:convert, MLInt{N, T}, x))
    end

    _unsafe_convert(MLInt{N, T}, x)
end

@inline function Base.convert(::Type{MLInt{N, T}}, x::Bool) where {N, T}
    x ? one(MLInt{N, T}) : zero(MLInt{N, T})
end


@inline Base.promote_type(::Type{MLInt{N, T}}, ::Type{<:Signed}) where {N, T} = MLInt{N, T}
@inline Base.promote_type(::Type{<:Signed}, ::Type{MLInt{N, T}}) where {N, T} = MLInt{N, T}
@inline Base.promote_type(::Type{MLInt{N, T}}, ::Type{MLInt{M, T}}) where {N, M, T} =
    MLInt{max(M, N), T}
@inline Base.promote_type(::Type{MLInt{N, T}}, ::Type{MLInt{N, T}}) where {N, T} =
    MLInt{N, T}


@inline Base.signbit(x::MLInt{N, T}) where {N, T} = signbit(signed(x[N]))


@inline Base.zero(::Type{MLInt{N, T}}) where {N, T} = signed(zero(MLUInt{N, T}))


@inline Base.one(::Type{MLInt{N, T}}) where {N, T} = signed(one(MLUInt{N, T}))


@inline @generated function Base.typemin(::Type{MLInt{N, T}}) where {N, T}
    x = one(T) << (bitsizeof(T) - 1)
    exprs = vcat([:(zero(T)) for i in 1:N-1], [:($x)])
    quote
        MLInt(tuple($(exprs...)))
    end
end


@inline @generated function Base.typemax(::Type{MLInt{N, T}}) where {N, T}
    x = (one(T) << (bitsizeof(T) - 1)) - one(T)
    exprs = vcat([:(typemax(T)) for i in 1:N-1], [:($x)])
    quote
        MLInt(tuple($(exprs...)))
    end
end


Base.trailing_zeros(x::MLInt{N, T}) where {N, T} = trailing_zeros(unsigned(x))


@inline function Base.setindex(x::MLInt{N, T}, v::T, i::Integer) where {N, T}
    MLInt(setindex(x.limbs, v, i))
end


@inline Base.getindex(x::MLInt{N, T}, i::Integer) where {N, T} = x.limbs[i]


Base.string(x::MLInt{N, T}) where {N, T} = "{S" * string(x.limbs) * "}"


Base.show(io::IO, x::MLInt{N, T}) where {N, T} = print(io, string(x))


Base.:-(x::MLInt{N, T}) where {N, T} = ~x + one(MLInt{N, T})


Base.:+(x::MLInt{N, T}, y::MLInt{N, T}) where {N, T} = signed(unsigned(x) + unsigned(y))


Base.:~(x::MLInt{N, T}) where {N, T} = signed(~unsigned(x))


Base.:flipsign(x::T, y::T) where T <: Union{MLInt, MLUInt} = signbit(x) ? -x : x


