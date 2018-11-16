rr_value(x::RRElem) = x.value
rr_value(x::RRElemMontgomery) = x.value

rr_representation(::Type{RRElem{T, M}}) where {T, M} = RRElem
rr_representation(::Type{RRElemMontgomery{T, M}}) where {T, M} = RRElemMontgomery

rr_base_type(::Type{RRElem{T, M}}) where {T, M} = T
rr_base_type(::Type{RRElemMontgomery{T, M}}) where {T, M} = T

rr_modulus(::RRElem{T, M}) where {T, M} = M
rr_modulus(::RRElemMontgomery{T, M}) where {T, M} = M
rr_modulus(::Type{RRElem{T, M}}) where {T, M} = M
rr_modulus(::Type{RRElemMontgomery{T, M}}) where {T, M} = M


function rr_modulus_simple(tp::Type{<:AbstractRRElem})
    m = rr_modulus(tp)
    convert(encompassing_type(typeof(m)), m)
end


"""
    change_modulus(new_modulus::Unsigned, x::RRElem{T, M}) where {T, M}

Change modulus to `new_modulus` (it will be converted to type `T`), keeping the value intact.
If the new modulus is smaller than the current one, the value will be taken modulo `new_modulus`.
"""
@inline function change_modulus(new_modulus::Unsigned, x::RRElem{T, M}) where {T, M}
    nm = convert(T, new_modulus)
    if nm >= M
        RRElem(x.value, nm, _no_conversion)
    else
        RRElem(mod(x.value, nm), nm, _no_conversion)
    end
end

"""
    change_modulus(new_modulus::Unsigned, p::Polynomial{T}) where T <: AbstractRRElem

Apply `change_modulus()` to every coefficient of the polynomial.
"""
@inline function change_modulus(
        new_modulus::Unsigned, p::Polynomial{T}) where T <: AbstractRRElem
    # Convert the modulus in advance so that it is not converted for each element separately
    nm = convert(rr_base_type(T), new_modulus)
    Polynomial(change_modulus.(nm, p.coeffs), p.negacyclic)
end


@doc """
    change_representation(
        new_repr::Union{RRElem, RRElemMontgomery},
        x::Union{RRElem{T, M}, RRElemMontgomery{T, M}}) where {T, M}

Change the representation of the given residue ring element to one of
[`RRElem`](@ref), [`RRElemMontgomery`](@ref).
""" change_representation()

@inline change_representation(::Type{RRElem}, x::RRElem{T, M}) where {T, M} = x
@inline change_representation(::Type{RRElemMontgomery}, x::RRElem{T, M}) where {T, M} =
    convert(RRElemMontgomery{T, M}, x)
@inline change_representation(
    ::Type{RRElemMontgomery}, x::RRElemMontgomery{T, M}) where {T, M} = x
@inline change_representation(
    ::Type{RRElem}, x::RRElemMontgomery{T, M}) where {T, M} = convert(RRElem{T, M}, x)

"""
    change_representation(new_repr, p::Polynomial{T}) where T <: AbstractRRElem

Apply `change_representation()` to every coefficient of the polynomial.
"""
@inline change_representation(new_repr, p::Polynomial{T}) where T <: AbstractRRElem =
    Polynomial(change_representation.(new_repr, p.coeffs), p.negacyclic)


@doc """
    change_base_type(::Type{V}, x::RRElem{T, M}) where {T, M, V <: Unsigned}

Change the base type (`T`) of the given residue ring element to the type `V`.
The modulus `M` must fit into `V`.
""" change_base_type()

@inline function change_base_type(::Type{V}, x::RRElem{T, M}) where {T, M, V <: Unsigned}
    @assert M <= typemax(V)
    RRElem(convert(V, x.value), convert(V, M), _no_conversion)
end

@inline function change_base_type(::Type{V}, x::RRElemMontgomery{T, M}) where {T, M, V <: Unsigned}
    @assert M <= typemax(V)
    RRElemMontgomery(convert(V, x.value), convert(V, M), _no_conversion)
end

"""
    change_base_type(::Type{V}, p::Polynomial{T}) where {T <: AbstractRRElem, V <: Unsigned}

Apply `change_base_type()` to every coefficient of the polynomial.
"""
@inline function change_base_type(
        tp::Type{V}, p::Polynomial{T}) where {T <: AbstractRRElem, V <: Unsigned}
    Polynomial(change_base_type.(tp, p.coeffs), p.negacyclic)
end


@inline function _change_modulus_proportional(
        new_modulus::Unsigned, x::T, old_modulus::T) where T <: Unsigned
    hi, lo = mulhilo(x, convert(T, new_modulus))
    q, r = divremhilo(hi, lo, old_modulus)
    if r >= old_modulus ÷ 2 + (isodd(old_modulus) ? one(T) : zero(T))
        q += one(T)
    end
    q
end


"""
    change_modulus_proportional(new_modulus::Unsigned, x::RRElem{T, M}) where {T, M}

Change the modulus of the given residue ring element, adjusting the value proportionally.
That is, `x_new = round(T, x * new_modulus / M)` (with the proper handling of an overflow).
The new modulus must be lower or equal to `M`.
"""
@inline function change_modulus_proportional(new_modulus::Unsigned, x::RRElem{T, M}) where {T, M}
    RRElem(
        _change_modulus_proportional(new_modulus, x.value, M),
        convert(T, new_modulus),
        _no_conversion)
end

"""
    change_modulus_proportional(new_modulus::Unsigned, x::Polynomial{T}) where {T <: RRElem}

Apply `change_modulus_proportional()` to every coefficient of the polynomial.
"""
@inline function change_modulus_proportional(
        new_modulus::Unsigned, p::Polynomial{T}) where {T <: RRElem}
    Polynomial(
        change_modulus_proportional.(new_modulus, p.coeffs),
        p.negacyclic)
end


"""
    change_length(new_length::Integer, p::Polynomial)

Change length of the polynomial, padding it with zeros.
The new length must be greater or equal to the current length.
"""
@inline function change_length(new_length::Integer, p::Polynomial{T}) where T
    Polynomial([p.coeffs; zeros(T, new_length - length(p))], p.negacyclic)
end
