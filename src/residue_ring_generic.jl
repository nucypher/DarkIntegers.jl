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


@inline function change_modulus(new_modulus::Unsigned, x::RRElem{T, M}) where {T, M}
    nm = convert(T, new_modulus)
    if nm >= M
        RRElem(x.value, nm, _no_conversion)
    else
        RRElem(mod(x.value, nm), nm, _no_conversion)
    end
end

@inline function change_modulus(
        new_modulus::Unsigned, p::Polynomial{T}) where T <: AbstractRRElem
    # Convert the modulus in advance so that it is not converted for each element separately
    nm = convert(rr_base_type(T), new_modulus)
    Polynomial(change_modulus.(nm, p.coeffs), p.negacyclic)
end


@inline change_representation(::Type{RRElem}, x::RRElem{T, M}) where {T, M} = x
@inline change_representation(::Type{RRElemMontgomery}, x::RRElem{T, M}) where {T, M} =
    convert(RRElemMontgomery{T, M}, x)
@inline change_representation(
    ::Type{RRElemMontgomery}, x::RRElemMontgomery{T, M}) where {T, M} = x
@inline change_representation(
    ::Type{RRElem}, x::RRElemMontgomery{T, M}) where {T, M} = convert(RRElem{T, M}, x)
@inline change_representation(new_repr, p::Polynomial) =
    Polynomial(change_representation.(new_repr, p.coeffs), p.negacyclic)


@inline function change_base_type(::Type{V}, x::RRElem{T, M}) where {T, M, V <: Unsigned}
    @assert M <= typemax(V)
    RRElem(convert(V, x.value), convert(V, M), _no_conversion)
end

@inline function change_base_type(::Type{V}, x::RRElemMontgomery{T, M}) where {T, M, V <: Unsigned}
    @assert M <= typemax(V)
    RRElemMontgomery(convert(V, x.value), convert(V, M), _no_conversion)
end


@inline function change_modulus_proportional(
        new_modulus::Unsigned, x::T, old_modulus::T) where T <: Unsigned

    # TODO: optimize
    xi = convert(BigInt, x)
    mi = convert(BigInt, old_modulus)

    # TODO: make it a purely integer algorithm
    convert(T, round(BigInt, xi * new_modulus / mi))
end


@inline function change_modulus_proportional(new_modulus::Unsigned, x::RRElem{T, M}) where {T, M}
    RRElem(
        change_modulus_proportional(new_modulus, x.value, M),
        convert(T, new_modulus),
        _no_conversion)
end
