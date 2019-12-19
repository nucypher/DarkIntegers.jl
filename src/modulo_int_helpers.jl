"""
    raw_value(x::AbstractModUInt{T, M})

Returns the raw value of type `T` stored in the object (without any conversions).
"""
raw_value(x::ModUInt) = x.value
raw_value(x::MgModUInt) = x.value


"""
    raw_value(x::AbstractModUInt{T, M})

Returns the value of type `T` stored in the object
(converted out of any special representation if necessary).
"""
value(x::ModUInt) = raw_value(x)
value(x::MgModUInt) = from_montgomery(x)


Base.eltype(::Type{ModUInt{T, M}}) where {T, M} = T
Base.eltype(::Type{MgModUInt{T, M}}) where {T, M} = T
Base.eltype(::T) where T <: AbstractModUInt = eltype(T)


"""
    modulus(::Type{AbstractModUInt{T, M}})
    modulus(x::AbstractModUInt{T, M})

Returns the modulus `M` (of type `T`) for a modulo integer object or type.
"""
modulus(::Type{ModUInt{T, M}}) where {T, M} = M
modulus(::Type{MgModUInt{T, M}}) where {T, M} = M
modulus(::T) where T <: AbstractModUInt = modulus(T)
