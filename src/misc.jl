
encompassing_type(tp::Type{<:Unsigned}) = tp

encompassing_type(::Type{UInt4}) = UInt8

function encompassing_type(tp::Type{MPNumber{N, T}}) where {N, T}
    total_size = cld(N * bitsizeof(T), 8)

    if total_size <= 1
        return UInt8
    elseif total_size <= 2
        return UInt16
    elseif total_size <= 4
        return UInt32
    elseif total_size <= 8
        return UInt64
    elseif total_size <= 16
        return UInt128
    else
        return BigInt
    end
end
