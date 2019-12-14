raw_value(x::ModUInt) = x.value
raw_value(x::MgModUInt) = x.value


Base.eltype(::Type{ModUInt{T, M}}) where {T, M} = T
Base.eltype(::Type{MgModUInt{T, M}}) where {T, M} = T
Base.eltype(::T) where T <: AbstractModUInt = eltype(T)


modulus(::Type{ModUInt{T, M}}) where {T, M} = M
modulus(::Type{MgModUInt{T, M}}) where {T, M} = M
modulus(::T) where T <: AbstractModUInt = modulus(T)


raw_value_as_builtin(x::AbstractModUInt) = convert(encompassing_type(eltype(x)), raw_value(x))


value_as_builtin(x::AbstractModUInt) = convert(encompassing_type(eltype(x)), convert(eltype(x), x))


modulus_as_builtin(::Type{T}) where T <: AbstractModUInt =
    convert(encompassing_type(eltype(T)), modulus(T))
modulus_as_builtin(x::T) where T <: AbstractModUInt = modulus_as_builtin(T)


"""
    change_modulus(new_modulus::Unsigned, x::ModUInt{T, M}) where {T, M}

Change modulus to `new_modulus` (it will be converted to type `T`), keeping the value intact.
If the new modulus is smaller than the current one, the value will be taken modulo `new_modulus`.
"""
@inline function change_modulus(new_modulus::Unsigned, x::ModUInt{T, M}) where {T, M}
    nm = convert(T, new_modulus)
    if nm >= M
        ModUInt(x.value, nm, _verbatim)
    else
        ModUInt(mod(x.value, nm), nm, _verbatim)
    end
end

"""
    change_modulus(new_modulus::Unsigned, p::Polynomial{T}) where T <: AbstractModUInt

Apply `change_modulus()` to every coefficient of the polynomial.
"""
@inline function change_modulus(
        new_modulus::Unsigned, p::Polynomial{T}) where T <: AbstractModUInt
    # Convert the modulus in advance so that it is not converted for each element separately
    nm = convert(eltype(T), new_modulus)
    Polynomial(change_modulus.(nm, p.coeffs), p.negacyclic)
end


@doc """
    change_representation(
        new_repr::Union{ModUInt, MgModUInt},
        x::Union{ModUInt{T, M}, MgModUInt{T, M}}) where {T, M}

Change the representation of the given residue ring element to one of
[`ModUInt`](@ref), [`MgModUInt`](@ref).
""" change_representation()

@inline change_representation(::Type{ModUInt}, x::ModUInt{T, M}) where {T, M} = x
@inline change_representation(::Type{MgModUInt}, x::ModUInt{T, M}) where {T, M} =
    convert(MgModUInt{T, M}, x)
@inline change_representation(
    ::Type{MgModUInt}, x::MgModUInt{T, M}) where {T, M} = x
@inline change_representation(
    ::Type{ModUInt}, x::MgModUInt{T, M}) where {T, M} = convert(ModUInt{T, M}, x)

"""
    change_representation(new_repr, p::Polynomial{T}) where T <: AbstractModUInt

Apply `change_representation()` to every coefficient of the polynomial.
"""
@inline change_representation(new_repr, p::Polynomial{T}) where T <: AbstractModUInt =
    Polynomial(change_representation.(new_repr, p.coeffs), p.negacyclic)


@doc """
    change_base_type(::Type{V}, x::ModUInt{T, M}) where {T, M, V <: Unsigned}

Change the base type (`T`) of the given residue ring element to the type `V`.
The modulus `M` must fit into `V`.
""" change_base_type()

@inline function change_base_type(::Type{V}, x::ModUInt{T, M}) where {T, M, V <: Unsigned}
    @assert M <= typemax(V)
    ModUInt(convert(V, x.value), convert(V, M), _verbatim)
end

@inline function change_base_type(::Type{V}, x::MgModUInt{T, M}) where {T, M, V <: Unsigned}
    @assert M <= typemax(V)
    MgModUInt(convert(V, x.value), convert(V, M), _verbatim)
end

"""
    change_base_type(::Type{V}, p::Polynomial{T}) where {T <: AbstractModUInt, V <: Unsigned}

Apply `change_base_type()` to every coefficient of the polynomial.
"""
@inline function change_base_type(
        tp::Type{V}, p::Polynomial{T}) where {T <: AbstractModUInt, V <: Unsigned}
    Polynomial(change_base_type.(tp, p.coeffs), p.negacyclic)
end


@inline function _rescale(
        new_max::Unsigned, x::T, old_max::T, round_result::Bool) where T <: Unsigned
    nm = convert(T, new_max)
    hi, lo = mulhilo(x, nm)
    q, r = divremhilo(hi, lo, old_max)
    if round_result
        if r >= old_max ÷ 2 + (isodd(old_max) ? one(T) : zero(T))
            q += one(T)
            if q == nm
                q = zero(T)
            end
        end
    end
    q
end


"""
    rescale(new_max::Unsigned, x::ModUInt{T, M}, round_result::Bool)

Rescale `x` proportionally to the range `[0, new_max)` (where `new_max <= M`).
Equivalent to `floor(x * new_max / M)` or `round(...)`, depending on the value of `round_result`.
If `round_result` is `true`, and the value if equal to `new_max` after rounding, it is set to 0.
"""
@inline function rescale(
        new_max::Unsigned, x::ModUInt{T, M}, round_result::Bool) where {T, M}
    ModUInt(_rescale(new_max, x.value, M, round_result), M, _verbatim)
end


"""
    rescale(new_max::Unsigned, p::Polynomial{T}, round_result::Bool) where T <: ModUInt

Apply `rescale()` to every coefficient of the polynomial.
"""
@inline function rescale(
        new_max::Unsigned, p::Polynomial{T}, round_result::Bool) where T <: ModUInt
    Polynomial(rescale.(new_max, p.coeffs, round_result), p.negacyclic)
end


"""
    change_length(new_length::Integer, p::Polynomial)

Change length of the polynomial, padding it with zeros.
The new length must be greater or equal to the current length.
"""
@inline function change_length(new_length::Integer, p::Polynomial{T}) where T
    old_length = length(p.coeffs)
    if new_length > old_length
        Polynomial([p.coeffs; zeros(T, new_length - length(p.coeffs))], p.negacyclic)
    elseif new_length < old_length
        Polynomial(p.coeffs[1:new_length], p.negacyclic)
    else
        p
    end
end
