"""
Multi-precision unsigned integers.
"""

using Base: setindex


struct MPNumber{N, T <: Unsigned} <: Unsigned
    value :: NTuple{N, T}

    MPNumber(x::NTuple{N, T}) where {N, T} = new{N, T}(x)

    @inline function MPNumber{N, T}(x::Integer) where {N, T}
        res = zero(MPNumber{N, T})
        for i in 1:N
            res = setindex(res, T(x & typemax(T)), i)
            x >>= bitsizeof(T)
        end
        res
    end
end


@inline Base.convert(::Type{MPNumber{N, T}}, x::MPNumber{N, T}) where {N, T} = x


@inline function Base.convert(::Type{V}, x::MPNumber{N, T}) where {V <: Integer, N, T}
    res = zero(V)
    for i in 1:N
        res += convert(V, x.value[i]) << (bitsizeof(T) * (i - 1))
    end
    res
end


# Because MPNumber <: Unsigned, show() will be bypassed sometimes in favor of string()
Base.string(x::MPNumber{N, T}) where {N, T} = "{" * string(x.value) * "}"


Base.show(io::IO, x::MPNumber{N, T}) where {N, T} = print(io, string(x))


@inline @generated function Base.zero(::Type{MPNumber{N, T}}) where {N, T}
    exprs = [:(zero(T)) for i in 1:N]
    quote
        MPNumber(tuple($(exprs...)))
    end
end


@inline @generated function Base.one(::Type{MPNumber{N, T}}) where {N, T}
    exprs = [i == 1 ? :(one(T)) : :(zero(T)) for i in 1:N]
    quote
        MPNumber(tuple($(exprs...)))
    end
end


@inline Base.typemin(::Type{MPNumber{N, T}}) where {N, T} = zero(MPNumber{N, T})


@inline @generated function Base.typemax(::Type{MPNumber{N, T}}) where {N, T}
    exprs = [:(typemax(T)) for i in 1:N]
    quote
        MPNumber(tuple($(exprs...)))
    end
end


@inline function Base.setindex(x::MPNumber{N, T}, v::T, i::Integer) where {N, T}
    MPNumber(setindex(x.value, v, i))
end


@inline Base.getindex(x::MPNumber{N, T}, i::Integer) where {N, T} = x.value[i]


# Addition of unsigned numbers with carry
@inline function _addc(x::T, y::T) where T <: Unsigned
    r = x + y
    r, r < x
end


# Subtraction of unsigned numbers with carry
@inline function _subc(x::T, y::T) where T <: Unsigned
    r = x - y
    r, x < y
end


@inline function Base.:+(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}
    c = false
    out = zero(MPNumber{N, T})
    for i in 1:N
        r, new_c = _addc(x.value[i], y.value[i])
        r, new_c2 = _addc(r, T(c))
        out = setindex(out, r, i)

        # `v[i] + c` is at most the limb size,
        # So we will have to carry at most 1.
        c = new_c || new_c2
    end
    out
end


@inline function Base.:-(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}
    c = false
    out = zero(MPNumber{N, T})
    for i in 1:N
        r, new_c = _subc(x.value[i], y.value[i])
        r, new_c2 = _subc(r, T(c))
        out = setindex(out, r, i)

        # `v[i] + c` is at most the limb size,
        # So we will have to carry at most 1.
        c = new_c || new_c2
    end
    out
end


@inline function Base.:>=(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}
    for i in N:-1:1
        if x.value[i] == y.value[i]
            continue
        end
        return x.value[i] > y.value[i]
    end
    true
end


@inline function Base.:>(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}
    for i in N:-1:1
        if x.value[i] == y.value[i]
            continue
        end
        return x.value[i] > y.value[i]
    end
    false
end


@inline function Base.:<=(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}
    for i in N:-1:1
        if x.value[i] == y.value[i]
            continue
        end
        return x.value[i] < y.value[i]
    end
    true
end


@inline function Base.:<(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}
    for i in N:-1:1
        if x.value[i] == y.value[i]
            continue
        end
        return x.value[i] < y.value[i]
    end
    false
end


@inline function _most_significant_limb(x::MPNumber{N, T}) where {N, T}
    for i in N:-1:1
        if x.value[i] > 0
            return i
        end
    end
    return 0
end


@inline function Base.:*(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}
    # TODO: to protect from timing attacks we can assume n == t == N
    # This will also allow us to use a generated function and may turn out to be faster...
    n = _most_significant_limb(x) - 1
    t = _most_significant_limb(y) - 1
    w = zero(MPNumber{N, T})
    for i in 1:t+1
        c = zero(T)
        hi = zero(T)
        for j in 1:min(n+1, N+1-i) # i + j - 1 <= N -> j <= N + 1 - i
            hi, lo = mulhilo(x.value[j], y.value[i])
            hi, lo = addhilo(hi, lo, c)
            hi, lo = addhilo(hi, lo, w.value[i + j - 1])
            w = setindex(w, lo, i + j - 1)
            c = hi
        end
        if c > 0 && i + n + 1 <= N
            hi, lo = addhilo(zero(T), c, w.value[i + n + 1])
            w = setindex(w, lo, i + n + 1)
        end
    end
    w
end


@inline function Base.isodd(x::MPNumber{N, T}) where {N, T}
    isodd(x.value[1])
end


function Base.:^(x::MPNumber{N, T}, y::Integer) where {N, T}
    # TODO: optimize
    res = one(MPNumber{N, T})
    @assert y >= 0
    for i in 1:y
        res *= x
    end
    res
end


function _ge_shift(x::MPNumber{N, T}, y::MPNumber{N, T}, y_shift::Integer) where {N, T}
    # true if x >= y * b^y_shift
    for i in N:-1:1+y_shift
        if x.value[i] == y.value[i - y_shift]
            continue
        end
        return x.value[i] > y.value[i - y_shift]
    end
    true
end


function _shift_limbs(x::MPNumber{N, T}, shift::Integer) where {N, T}
    # x -> x * b^shift
    res = zero(MPNumber{N, T})
    for i in 1:shift
        res = setindex(res, zero(T), i)
    end
    for i in shift+1:N
        res = setindex(res, x[i - shift], i)
    end
    res
end


# Calculates y * (x0 + x1 * b) -> r0, r1, r2
function _mul_1d_2d(x0::T, x1::T, y::T) where T <: Unsigned
    # Adapted from the generic multiplication for radix integers
    hi, lo = mulhilo(x0, y)
    w0 = lo
    c = hi

    hi, lo = mulhilo(x1, y)
    hi, lo = addhilo(hi, lo, c)
    w1 = lo
    w2 = hi

    (w0, w1, w2)
end


# Calculates x - y * z
function _sub_mul(x::MPNumber{N, T}, y::T, z::MPNumber{N, T}) where {N, T}
    hi_carry1 = false
    hi_carry2 = false
    hi_carry3 = false
    hi = zero(T)
    out = zero(MPNumber{N, T})
    for j in 1:N
        new_hi, lo = mulhilo(y, z.value[j])

        r, new_hi_carry1 = _subc(x.value[j], lo)
        r, new_hi_carry2 = _subc(r, hi)
        r, new_hi_carry3 = _subc(r, T(hi_carry1 + hi_carry2 + hi_carry3))

        out = setindex(out, r, j)
        hi = new_hi
        hi_carry1 = new_hi_carry1
        hi_carry2 = new_hi_carry2
        hi_carry3 = new_hi_carry3
    end

    out_hi, c1 = _subc(zero(T), hi)
    out_hi, c2 = _subc(out_hi, T(hi_carry1 + hi_carry2 + hi_carry3))
    out, c1 || c2 # TODO: it seems that we can have at most 1 carried over to the n+2-th limb
end


function divrem_single_limb(x::MPNumber{N, T}, y::T) where {N, T}
    r = zero(T)
    q = zero(MPNumber{N, T})
    for j in N-1:-1:0
        d, r = divremhilo(r, x[j+1], y)
        q = setindex(q, d, j+1)
    end
    q, setindex(zero(MPNumber{N, T}), r, 1)
end


function Base.divrem(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}
    n = _most_significant_limb(x) - 1
    t = _most_significant_limb(y) - 1

    if n < t
        return zero(MPNumber{N, T}), x
    end

    q = zero(MPNumber{N, T})
    r = zero(MPNumber{N, T})

    if t == 0 && n == 0
        q1, r1 = divrem(x[1], y[1])
        return setindex(q, q1, 1), setindex(r, r1, 1)
    end

    if t == 0
        return divrem_single_limb(x, y[1])
    end

    # TODO: can be replaced by `<` and `shift_limbs`
    while _ge_shift(x, y, n - t)
        q = setindex(q, q[n - t + 1] + one(T), n - t + 1)
        x = x - _shift_limbs(y, n - t)
    end

    for i in n:-1:t+1
        if x[i+1] == y[t+1]
            q = setindex(q, typemax(T), i - t - 1 + 1)
        else
            q = setindex(q, divhilo(x[i+1], x[i-1+1], y[t+1]), i - t - 1 + 1)
        end

        while (MPNumber(_mul_1d_2d(y[t-1+1], y[t+1], q[i-t-1+1]))
                > MPNumber((x[i-2+1], x[i-1+1], x[i+1])))
            q = setindex(q, q[i - t - 1+1] - one(T), i - t - 1+1)
        end

        x, c = _sub_mul(x, q[i-t-1+1], _shift_limbs(y, i - t - 1))
        if c
            x = x + _shift_limbs(y, i - t - 1)
            q = setindex(q, q[i - t - 1+1] - one(T), i - t - 1+1)
        end
    end

    (q, x)
end


@inline function Base.div(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}
    # TODO: is there a faster way?
    d, r = divrem(x, y)
    d
end


bitsizeof(::Type{MPNumber{N, T}}) where {N, T} = bitsizeof(T) * N


# Required for broadcasting


Base.length(x::MPNumber{N, T}) where {N, T} = 1


Base.iterate(x::MPNumber{N, T}) where {N, T} = (x, nothing)
Base.iterate(x::MPNumber{N, T}, state) where {N, T} = nothing
