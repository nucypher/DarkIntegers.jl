# Taken from https://github.com/JuliaLang/julia/pull/30515
# as a workaround for https://github.com/JuliaLang/julia/issues/29971
# Replace by `invmod` when it is merged.

hastypemax(::Base.BitIntegerType) = true
hastypemax(::Type{T}) where {T} = applicable(typemax, T)

function invmod_(n::T, m::T) where T <: Integer
    g, x, y = gcdx(n, m)
    g != 1 && throw(DomainError((n, m), "Greatest common divisor is $g."))
    m == 0 && throw(DomainError(m, "`m` must not be 0."))
    # Note that m might be negative here.
    r = (T <: Unsigned && hastypemax(typeof(x)) && x > typemax(x)>>1) ? mod(x + m, m) : mod(x, m)
    # The postcondition is: mod(r * n, m) == mod(T(1), m) && div(r, m) == 0
    r
end


@inline halve(x::Integer) = x >> 1


@inline double(x::Integer) = x << 1


bitsizeof(tp) = sizeof(tp) << 3
