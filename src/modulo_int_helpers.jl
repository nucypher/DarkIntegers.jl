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


"""
    change_length(new_length::Integer, p::Polynomial)

Change length of the polynomial, padding it with zeros.
The new length must be greater or equal to the current length.
"""
@inline function change_length(new_length::Integer, p::Polynomial{T}) where T
    old_length = length(p.coeffs)
    if new_length > old_length
        Polynomial([p.coeffs; zeros(T, new_length - length(p.coeffs))], p.modulus)
    elseif new_length < old_length
        Polynomial(p.coeffs[1:new_length], p.modulus)
    else
        p
    end
end
