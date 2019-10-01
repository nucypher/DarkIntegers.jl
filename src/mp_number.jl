"""
Multi-precision unsigned integers.
"""

using Base: setindex


"""
Multi-precision unsigned integer type, with `N` limbs of type `T`
(which must be an unsigned integer type).

Supports `+`, `-`, `*`, `divrem`, `div`, `rem`, `^`, `<`, `<=`, `>`, `>=`,
`zero`, `one` and `isodd`.

    MPNumber{N, T}(x::Integer) where {N, T <: Unsigned}

Creates an `MPNumber` object. If `x` does not fit into `N` limbs of type `T`,
the excess bits will be ignored.
"""
struct MPNumber{N, T <: Unsigned} <: Unsigned
    value :: NTuple{N, T}

    MPNumber(x::NTuple{N, T}) where {N, T} = new{N, T}(x)
    MPNumber{N, T}(x::NTuple{N, T}) where {N, T} = new{N, T}(x)

    @inline function MPNumber{N, T}(x::Integer) where {N, T <: Unsigned}
        res = zero(MPNumber{N, T})
        for i in 1:N
            res = setindex(res, T(x & typemax(T)), i)
            x >>= bitsizeof(T)
        end
        res
    end

    @inline function MPNumber{N, T}(x::MPNumber{M, T}) where {N, M, T <: Unsigned}
        res = zero(MPNumber{N, T})
        for i in 1:min(N, M)
            res = setindex(res, x[i], i)
        end
        res
    end
end


# These are required to prevent the more general conversion to any integer from triggering.
@inline Base.convert(::Type{MPNumber{N, T}}, x::MPNumber{N, T}) where {N, T} = x
@inline Base.convert(::Type{MPNumber{N, T}}, x::MPNumber{M, T}) where {N, M, T} = MPNumber{N, T}(x)

@inline function Base.convert(::Type{V}, x::MPNumber{N, T}) where {V <: Integer, N, T}
    res = zero(V)
    for i in 1:N
        res += convert(V, x.value[i]) << (bitsizeof(T) * (i - 1))
    end
    res
end


@inline Base.promote_type(::Type{MPNumber{N, T}}, ::Type{<:Integer}) where {N, T} = MPNumber{N, T}
@inline Base.promote_type(::Type{<:Integer}, ::Type{MPNumber{N, T}}) where {N, T} = MPNumber{N, T}
@inline Base.promote_type(::Type{MPNumber{N, T}}, ::Type{MPNumber{M, T}}) where {N, M, T} =
    MPNumber{max(M, N), T}
@inline Base.promote_type(::Type{MPNumber{N, T}}, ::Type{MPNumber{N, T}}) where {N, T} =
    MPNumber{N, T}


# We need this to correctly process arithmetic operations on MPNumber and Int
# (which is signed and the default in Julia for number literals)
# without defining specific methods for each operator.
@inline Base.signed(x::MPNumber{N, T}) where {N, T} = x
@inline Base.unsigned(x::MPNumber{N, T}) where {N, T} = x


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


@inline Base.:-(x::MPNumber{N, T}) where {N, T} = zero(MPNumber{N, T}) - x


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


@inline function Base.:>>(x::MPNumber{N, T}, y::Int) where {N, T}
    res = zero(MPNumber{N, T})

    if y >= N * bitsizeof(T)
        return res
    end

    msl = _most_significant_limb(x)
    if msl == 0
        return res
    end

    full_limbs = y ÷ bitsizeof(T)
    remainder = y % bitsizeof(T)

    if remainder == 0
        @inbounds for i in 1:(msl - full_limbs)
            res = setindex(res, x[i + full_limbs], i)
        end
    else
        @inbounds for i in 1:(msl - full_limbs)
            lo = x[i + full_limbs] >> remainder
            if i < msl - full_limbs
                lo |= x[i + full_limbs + 1] << (bitsizeof(T) - remainder)
            end
            res = setindex(res, lo, i)
        end
    end

    res
end


@inline function Base.trailing_zeros(x::MPNumber{N, T}) where {N, T}
    bs = bitsizeof(T)
    res = 0
    @inbounds for i in 1:N
        t = trailing_zeros(x[i])
        res += t
        if t != bs
            break
        end
    end
    res
end


@inline function Base.isodd(x::MPNumber{N, T}) where {N, T}
    isodd(x.value[1])
end


@inline function divrem_single_limb(x::MPNumber{N, T}, y::T) where {N, T}
    r = zero(T)
    q = zero(MPNumber{N, T})
    for j in N-1:-1:0
        d, r = divremhilo(r, x[j+1], y)
        q = setindex(q, d, j+1)
    end
    q, setindex(zero(MPNumber{N, T}), r, 1)
end


# Shift a two-limb number `(hi, lo)` by `shift` bits to the left,
# returning the most significant half.
@inline function _lshift_through(hi::T, lo::T, shift::Integer) where T
    (hi << shift) | (lo >> (bitsizeof(T) - shift))
end


# Returns a pair of `x[x_start:x_start+z_end] - y * z[1:z_end]`
# and the overflow flag for the subtraction.
# Assumes that `x_start+z_end <= N+1` and `z_end <= N`.
@inline function _mul_sub_from_part(
        x::MPNumber{N, T}, x_start, y::T, z::MPNumber{N, T}, z_end) where {N, T}

    hi_carry1 = false
    hi_carry2 = false
    hi_carry3 = false
    hi = zero(T)
    for j in 1:z_end
        new_hi, lo = mulhilo(y, z[j])

        r, new_hi_carry1 = _subc(x[j + x_start - 1], lo)
        r, new_hi_carry2 = _subc(r, hi)
        r, new_hi_carry3 = _subc(r, T(hi_carry1 + hi_carry2 + hi_carry3))

        x = setindex(x, r, j + x_start - 1)

        hi = new_hi
        hi_carry1 = new_hi_carry1
        hi_carry2 = new_hi_carry2
        hi_carry3 = new_hi_carry3
    end

    last_idx = z_end + x_start
    x_hi = last_idx <= N ? x[last_idx] : zero(T)

    out_hi, c1 = _subc(x_hi, hi)
    out_hi, c2 = _subc(out_hi, T(hi_carry1 + hi_carry2 + hi_carry3))

    if last_idx <= N
        x = setindex(x, out_hi, last_idx)
    end

    x, c1 || c2 # TODO: it seems that we can have at most 1 carried over to the n+2-th limb
end


# Returns `x[x_start:x_start+y_end] + y[1:y_end]`.
# Assumes `x_start+y_end <= N+1` and `y_end <= N`.
@inline function _add_to_part(x::MPNumber{N, T}, x_start, y::MPNumber{N, T}, y_end) where {N, T}

    c = false
    for i in 1:y_end
        r, new_c = _addc(x[i+x_start-1], y[i])
        r, new_c2 = _addc(r, T(c))
        x = setindex(x, r, i+x_start-1)

        # `v[i] + c` is at most the limb size,
        # So we will have to carry at most 1.
        c = new_c || new_c2
    end

    last_idx = y_end + x_start
    if last_idx <= N
        x = setindex(x, x[last_idx] + T(c), last_idx)
    end

    x
end


@inline function Base.divrem(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}

    # Division algorithm from D. Knuth's "The Art of Computer Programming", vol. 2.

    t = _most_significant_limb(x)
    n = _most_significant_limb(y)

    if n > t || y > x
        return zero(MPNumber{N, T}), x
    end

    if n == 1
        return divrem_single_limb(x, y[1])
    end

    m = t - n
    q = zero(MPNumber{N, T})

    # For efficiency, we are normalizing the divisor so that the most significant bit
    # of its most significant limb is set.

    # For the normalization, we need the bitshift `l` such that `(y[n] << l) in [b/2, b)`,
    # where `b` is the radix (`typemax(T) + 1`).
    # Naturally, it is exactly `leading_zeros(v[n])`.
    l = leading_zeros(y[n])

    # The normalization is "virtual": we're not actually changing `x` and `y`,
    # but instead just calculate "shifted" versions of some of their limbs.

    # `y` has at least two limbs; if n = 2 we are taking the lowest limb to be `0`.
    y1 = y[n]
    y2 = y[n-1]
    y3 = n >= 3 ? y[n-2] : zero(T)

    y_p = _lshift_through(y1, y2, l)
    y_pp = _lshift_through(y2, y3, l)

    for j = m:-1:0

        # We need to find the quotient `q = x[1+j:1+j+n] / y`
        # It can be proved that `q_hat = (x[1+j+n] * b + x[j+n]) / y[n]`
        # is between `q` and `q + 2` (providing that the highest big of `y[n]` is set).
        # So we are calculating it, and then do some tests to adjust it.

        # `x` has at least two limbs; take the values out of bounds to be `0`.
        x0 = j < m ? x[n+j+1] : zero(T)
        x1 = x[n+j]
        x2 = x[n+j-1]
        x3 = n+j >= 3 ? x[n+j-2] : zero(T)

        x_p = _lshift_through(x0, x1, l)
        x_pp = _lshift_through(x1, x2, l)
        x_ppp = _lshift_through(x2, x3, l)

        # Calculate the initial approximation
        q_hat, r_hat, overflow = divremhilo(x_p, x_pp, y_p)
        if overflow
            q_hat = typemax(T)
            r_hat = x_pp
        end

        if overflow
            r_hat, overflow2 = _addc(r_hat, y_p)
        else
            overflow2 = false
        end

        if !overflow2
            # TODO: `q_hat` differs from the target value by at most 2,
            # so this loop can theoretically be unrolled.
            while true
                thi, tlo = mulhilo(q_hat, y_pp)
                if MPNumber((tlo, thi)) <= MPNumber((x_ppp, r_hat))
                    break
                end

                q_hat -= one(T)

                r_hat, c2 = _addc(r_hat, y_p)
                if c2
                    break
                end
            end
        end

        x, overflow3 = _mul_sub_from_part(x, j+1, q_hat, y, n)

        if overflow3
            q_hat -= one(T)
            x = _add_to_part(x, j+1, y, n)
        end

        q = setindex(q, q_hat, 1+j)
    end

    q, x
end


@inline function Base.div(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}
    # TODO: is there a faster way?
    q, r = divrem(x, y)
    q
end


@inline function Base.rem(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}
    # TODO: is there a faster way?
    q, r = divrem(x, y)
    r
end


bitsizeof(::Type{MPNumber{N, T}}) where {N, T} = bitsizeof(T) * N


# Required for broadcasting


Base.length(x::MPNumber{N, T}) where {N, T} = 1


Base.iterate(x::MPNumber{N, T}) where {N, T} = (x, nothing)
Base.iterate(x::MPNumber{N, T}, state) where {N, T} = nothing


# The following methods are needed for MPNumber to support mulmod_bitshift()


@inline function double(x::MPNumber{N, T}) where {N, T}
    msb = false
    for i in 1:N
        new_msb = leading_zeros(x[i]) == 0
        x = setindex(x, (x[i] << 1) | ifelse(msb, one(T), zero(T)), i)
        msb = new_msb
    end
    x
end


_msb(::Type{UInt4}) = UInt4(0x8)
_msb(::Type{UInt8}) = UInt8(0x80)
_msb(::Type{UInt16}) = UInt16(0x8000)
_msb(::Type{UInt32}) = UInt32(0x80000000)
_msb(::Type{UInt64}) = UInt64(0x8000000000000000)


@inline function halve(x::MPNumber{N, T}) where {N, T}
    lsb = false
    for i in N:-1:1
        new_lsb = trailing_zeros(x[i]) == 0
        x = setindex(x, (x[i] >> 1) | ifelse(lsb, _msb(T), zero(T)), i)
        lsb = new_lsb
    end
    x
end


# The following methods are needed for MPNumber to support mulmod_widemul()


Base.widen(::Type{MPNumber{N, T}}) where {N, T} = MPNumber{N*2, T}


function Base.widemul(x::MPNumber{N, T}, y::MPNumber{N, T}) where {N, T}
    n = _most_significant_limb(x) - 1
    t = _most_significant_limb(y) - 1
    w = zero(MPNumber{N*2, T})
    for i in 1:t+1
        c = zero(T)
        hi = zero(T)
        for j in 1:n+1
            hi, lo = mulhilo(x.value[j], y.value[i])
            hi, lo = addhilo(hi, lo, c)
            hi, lo = addhilo(hi, lo, w.value[i + j - 1])
            w = setindex(w, lo, i + j - 1)
            c = hi
        end
        w = setindex(w, c, i + n + 1)
    end
    w
end


@inline mulmod(x::MPNumber{N, T}, y::MPNumber{N, T}, modulus::MPNumber{N, T}) where {N, T} =
    mulmod_widemul(x, y, modulus)
