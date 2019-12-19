raw_value(x::ModUInt) = x.value
raw_value(x::MgModUInt) = x.value


value(x::ModUInt) = raw_value(x)
value(x::MgModUInt) = from_montgomery(x)


Base.eltype(::Type{ModUInt{T, M}}) where {T, M} = T
Base.eltype(::Type{MgModUInt{T, M}}) where {T, M} = T
Base.eltype(::T) where T <: AbstractModUInt = eltype(T)


modulus(::Type{ModUInt{T, M}}) where {T, M} = M
modulus(::Type{MgModUInt{T, M}}) where {T, M} = M
modulus(::T) where T <: AbstractModUInt = modulus(T)
