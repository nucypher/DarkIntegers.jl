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


function rr_value_simple(val::AbstractRRElem)
    v = rr_value(val)
    convert(encompassing_type(typeof(v)), v)
end


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
        RRElem(x.value, nm, _verbatim)
    else
        RRElem(mod(x.value, nm), nm, _verbatim)
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
    RRElem(convert(V, x.value), convert(V, M), _verbatim)
end

@inline function change_base_type(::Type{V}, x::RRElemMontgomery{T, M}) where {T, M, V <: Unsigned}
    @assert M <= typemax(V)
    RRElemMontgomery(convert(V, x.value), convert(V, M), _verbatim)
end

"""
    change_base_type(::Type{V}, p::Polynomial{T}) where {T <: AbstractRRElem, V <: Unsigned}

Apply `change_base_type()` to every coefficient of the polynomial.
"""
@inline function change_base_type(
        tp::Type{V}, p::Polynomial{T}) where {T <: AbstractRRElem, V <: Unsigned}
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
    rescale(new_max::Unsigned, x::RRElem{T, M}, round_result::Bool)

Rescale `x` proportionally to the range `[0, new_max)` (where `new_max <= M`).
Equivalent to `floor(x * new_max / M)` or `round(...)`, depending on the value of `round_result`.
If `round_result` is `true`, and the value if equal to `new_max` after rounding, it is set to 0.
"""
@inline function rescale(
        new_max::Unsigned, x::RRElem{T, M}, round_result::Bool) where {T, M}
    RRElem(_rescale(new_max, x.value, M, round_result), M, _verbatim)
end


"""
    rescale(new_max::Unsigned, p::Polynomial{T}, round_result::Bool) where T <: RRElem

Apply `rescale()` to every coefficient of the polynomial.
"""
@inline function rescale(
        new_max::Unsigned, p::Polynomial{T}, round_result::Bool) where T <: RRElem
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
