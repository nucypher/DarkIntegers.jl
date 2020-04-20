#=
Multi-limb unsigned integers.
=#


"""
Multi-precision unsigned integer type, with `N` limbs of type `T`
(which must be an unsigned integer type).

Supports `+`, `-`, `*`, `divrem`, `div`, `rem`, `mod`, `^`,
`==`, `!=`, `<`, `<=`, `>`, `>=`,
`zero`, `one`, `oneunit`, `iseven`, `isodd`, `typemin`, `typemax`, `iszero`,
`sizeof` (can be off if limbs have the type `UInt4`),
`bitsizeof`, `leading_zeros`, `trailing_zeros`, `eltype`, `abs`.

Also supports `num_bits`, `halve`, `double`, `encompassing_type`.

The objects can be created either with `convert()`, or as

    MLUInt(x::NTuple{N, T})
    MLUInt{N, T}(x::NTuple{N, V})
"""
struct MLUInt{N, T <: Unsigned} <: Unsigned
    limbs :: NTuple{N, T}

    MLUInt(x::NTuple{N, T}) where {N, T} = new{N, T}(x)
    MLUInt{N, T}(x::MLUInt{N, T}) where {N, T} = new{N, T}(x.limbs)
    MLUInt{N, T}(x::NTuple{N, V}) where {N, T, V} = new{N, T}(convert.(T, x))
end


@inline function _most_significant_limb(x::MLUInt{N, T}) where {N, T}
    for i in N:-1:1
        if !iszero(x[i])
            return i
        end
    end
    return 0
end


@inline function _unsafe_convert(::Type{MLUInt{N, T}}, x::MLUInt{M, T}) where {N, M, T}
    res = zero(MLUInt{N, T})
    for i in 1:min(N, M)
        res = setindex(res, x[i], i)
    end
    res
end

@inline function _unsafe_convert(::Type{V}, x::MLUInt{N, T}) where {V <: Integer, N, T}
    res = zero(V)
    for i in 1:N
        res |= convert(V, x[i]) << (bitsizeof(T) * (i - 1))
    end
    res
end

@inline function _unsafe_convert(::Type{MLUInt{N, T}}, x::Integer) where {N, T}
    res = zero(MLUInt{N, T})
    for i in 1:N
        res = setindex(res, convert(T, x & typemax(T)), i)
        x >>= bitsizeof(T)
        if iszero(x)
            break
        end
    end
    res
end


# These are required to prevent the more general conversion to any integer from triggering.
@inline Base.convert(::Type{MLUInt{N, T}}, x::MLUInt{N, T}) where {N, T} = x

@inline function Base.convert(::Type{MLUInt{N, T}}, x::MLUInt{M, T}) where {N, M, T}
    if _most_significant_limb(x) > N
        throw(InexactError(:convert, MLUInt{N, T}, x))
    end
    _unsafe_convert(MLUInt{N, T}, x)
end

@inline function Base.convert(::Type{V}, x::MLUInt{N, T}) where {V <: Signed, N, T}
    # `>=` since one bit of `V` will be reserved for the sign
    if V != BigInt && num_bits(x) >= bitsizeof(V)
        throw(InexactError(:convert, V, x))
    end
    _unsafe_convert(V, x)
end

@inline function Base.convert(::Type{V}, x::MLUInt{N, T}) where {V <: Unsigned, N, T}
    if num_bits(x) > bitsizeof(V)
        throw(InexactError(:convert, V, x))
    end
    _unsafe_convert(V, x)
end

@inline function Base.convert(::Type{MLUInt{N, T}}, x::Integer) where {N, T}
    # Similar behaviour to built-in Julia types - can't convert negative numbers to unsigned.
    if signbit(x)
        throw(InexactError(:convert, MLUInt{N, T}, x))
    # Check if `x` fits in the target type
    elseif num_bits(x) > bitsizeof(MLUInt{N, T})
        throw(InexactError(:convert, MLUInt{N, T}, x))
    end
    _unsafe_convert(MLUInt{N, T}, x)
end

@inline function Base.convert(::Type{MLUInt{N, T}}, x::Bool) where {N, T}
    x ? one(MLUInt{N, T}) : zero(MLUInt{N, T})
end


@inline Base.promote_type(::Type{MLUInt{N, T}}, ::Type{<:Unsigned}) where {N, T} = MLUInt{N, T}
@inline Base.promote_type(::Type{MLUInt{N, T}}, ::Type{<:Signed}) where {N, T} = MLUInt{N, T}
@inline Base.promote_type(::Type{<:Unsigned}, ::Type{MLUInt{N, T}}) where {N, T} = MLUInt{N, T}
@inline Base.promote_type(::Type{<:Signed}, ::Type{MLUInt{N, T}}) where {N, T} = MLUInt{N, T}
@inline Base.promote_type(::Type{MLUInt{N, T}}, ::Type{MLUInt{M, T}}) where {N, M, T} =
    MLUInt{max(M, N), T}
@inline Base.promote_type(::Type{MLUInt{N, T}}, ::Type{MLUInt{N, T}}) where {N, T} =
    MLUInt{N, T}


@inline Base.signed(x::MLUInt{N, T}) where {N, T} = MLInt{N, T}(x.limbs)


@inline Base.unsigned(x::MLUInt{N, T}) where {N, T} = x


@inline Base.signbit(::MLUInt{N, T}) where {N, T} = false


# Because MLUInt <: Unsigned, show() will be bypassed sometimes in favor of string()
Base.string(x::MLUInt{N, T}) where {N, T} = "{" * string(x.limbs) * "}"


Base.show(io::IO, x::MLUInt{N, T}) where {N, T} = print(io, string(x))


@inline @generated function Base.zero(::Type{MLUInt{N, T}}) where {N, T}
    exprs = [:(zero(T)) for i in 1:N]
    quote
        MLUInt(tuple($(exprs...)))
    end
end


@inline @generated function Base.one(::Type{MLUInt{N, T}}) where {N, T}
    exprs = [i == 1 ? :(one(T)) : :(zero(T)) for i in 1:N]
    quote
        MLUInt(tuple($(exprs...)))
    end
end


@inline Base.typemin(::Type{MLUInt{N, T}}) where {N, T} = zero(MLUInt{N, T})


@inline @generated function Base.typemax(::Type{MLUInt{N, T}}) where {N, T}
    exprs = [:(typemax(T)) for i in 1:N]
    quote
        MLUInt(tuple($(exprs...)))
    end
end


@inline function Base.setindex(x::MLUInt{N, T}, v::T, i::Integer) where {N, T}
    MLUInt(setindex(x.limbs, v, i))
end


@inline Base.getindex(x::MLUInt{N, T}, i::Integer) where {N, T} = x.limbs[i]



@inline function Base.:+(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    c = false
    out = zero(MLUInt{N, T})
    for i in 1:N
        r, new_c = _addc(x[i], y[i])
        r, new_c2 = _addc(r, T(c))
        out = setindex(out, r, i)

        # `v[i] + c` is at most the limb size,
        # So we will have to carry at most 1.
        c = new_c || new_c2
    end
    out
end


@inline function Base.:-(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    c = false
    out = zero(MLUInt{N, T})
    for i in 1:N
        r, new_c = _subc(x[i], y[i])
        r, new_c2 = _subc(r, T(c))
        out = setindex(out, r, i)

        # `v[i] + c` is at most the limb size,
        # So we will have to carry at most 1.
        c = new_c || new_c2
    end
    out
end


@inline Base.:-(x::MLUInt{N, T}) where {N, T} = zero(MLUInt{N, T}) - x


@inline function Base.:>=(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    for i in N:-1:1
        if x[i] == y[i]
            continue
        end
        return x[i] > y[i]
    end
    true
end


@inline function Base.:>(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    for i in N:-1:1
        if x[i] == y[i]
            continue
        end
        return x[i] > y[i]
    end
    false
end


@inline function Base.:<=(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    for i in N:-1:1
        if x[i] == y[i]
            continue
        end
        return x[i] < y[i]
    end
    true
end


@inline function Base.:<(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    for i in N:-1:1
        if x[i] == y[i]
            continue
        end
        return x[i] < y[i]
    end
    false
end


@inline function Base.:*(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    # TODO: (see issue #14) to protect from timing attacks we can assume n == t == N
    # This will also allow us to use a generated function and may turn out to be faster...
    n = _most_significant_limb(x) - 1
    t = _most_significant_limb(y) - 1
    w = zero(MLUInt{N, T})
    for i in 1:t+1
        c = zero(T)
        hi = zero(T)
        for j in 1:min(n+1, N+1-i) # i + j - 1 <= N -> j <= N + 1 - i
            hi, lo = mulhilo(x[j], y[i])
            hi, lo = addhilo(hi, lo, c)
            hi, lo = addhilo(hi, lo, w[i + j - 1])
            w = setindex(w, lo, i + j - 1)
            c = hi
        end
        if c > 0 && i + n + 1 <= N
            hi, lo = addhilo(zero(T), c, w[i + n + 1])
            w = setindex(w, lo, i + n + 1)
        end
    end
    w
end


@inline function Base.:>>(x::MLUInt{N, T}, shift::Int) where {N, T}

    if shift == 0 || iszero(x)
        return x
    end

    if signbit(shift)
        return x << (-shift)
    end

    res = zero(MLUInt{N, T})

    if shift >= bitsizeof(MLUInt{N, T})
        return res
    end

    lb = log_bitsizeof(T)
    full_shifts = shift >> lb
    shift = xor(shift, full_shifts << lb)

    msl = _most_significant_limb(x)

    if shift == 0
        @inbounds for i in 1:(msl - full_shifts)
            res = setindex(res, x[i + full_shifts], i)
        end
    else
        @inbounds for i in 1:(msl - full_shifts)
            lo = x[i + full_shifts] >> shift
            if i < msl - full_shifts
                lo |= x[i + full_shifts + 1] << (bitsizeof(T) - shift)
            end
            res = setindex(res, lo, i)
        end
    end

    res
end


function Base.:<<(x::MLUInt{N, T}, shift::Int) where {N, T}

    if shift == 0 || iszero(x)
        return x
    end

    if signbit(shift)
        return x >> (-shift)
    end

    res = zero(MLUInt{N, T})

    if shift >= bitsizeof(MLUInt{N, T})
        return res
    end

    lb = log_bitsizeof(T)
    full_shifts = shift >> lb
    shift = xor(shift, full_shifts << lb)
    c_shift = bitsizeof(T) - shift

    for i in full_shifts+1:N
        hi = x[i-full_shifts] << shift
        lo = i == full_shifts+1 ? zero(T) : (x[i-full_shifts-1] >> c_shift)
        res = setindex(res, hi | lo, i)
    end

    res
end


@inline function Base.trailing_zeros(x::MLUInt{N, T}) where {N, T}
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


@inline function Base.leading_zeros(x::MLUInt{N, T}) where {N, T}
    if iszero(x)
        bitsizeof(MLUInt{N, T})
    else
        msl = _most_significant_limb(x)
        leading_zeros(x[msl]) + (N - msl) * bitsizeof(T)
    end
end


@inline function Base.iszero(x::MLUInt{N, T}) where {N, T}
    @inbounds for i in 1:N
        if !iszero(x[i])
            return false
        end
    end
    return true
end


@inline Base.isodd(x::MLUInt{N, T}) where {N, T} = isodd(x[1])


@inline Base.iseven(x::MLUInt{N, T}) where {N, T} = iseven(x[1])


@inline function divrem_single_limb(x::MLUInt{N, T}, y::T) where {N, T}
    r = zero(T)
    q = zero(MLUInt{N, T})
    for j in N-1:-1:0
        d, r = divremhilo(r, x[j+1], y)
        q = setindex(q, d, j+1)
    end
    q, setindex(zero(MLUInt{N, T}), r, 1)
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
        x::MLUInt{N, T}, x_start, y::T, z::MLUInt{N, T}, z_end) where {N, T}

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

    # TODO: (issue #23) it seems that we can have at most 1 carried over to the n+2-th limb
    x, c1 || c2
end


# Returns `x[x_start:x_start+y_end] + y[1:y_end]`.
# Assumes `x_start+y_end <= N+1` and `y_end <= N`.
@inline function _add_to_part(x::MLUInt{N, T}, x_start, y::MLUInt{N, T}, y_end) where {N, T}

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


@inline function Base.divrem(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}

    # Division algorithm from D. Knuth's "The Art of Computer Programming", vol. 2.

    t = _most_significant_limb(x)
    n = _most_significant_limb(y)

    if n > t || y > x
        return zero(MLUInt{N, T}), x
    end

    if n == 1
        return divrem_single_limb(x, y[1])
    end

    m = t - n
    q = zero(MLUInt{N, T})

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
            # TODO: (issue #24) `q_hat` differs from the target value by at most 2,
            # so this loop can theoretically be unrolled.
            while true
                thi, tlo = mulhilo(q_hat, y_pp)
                if MLUInt((tlo, thi)) <= MLUInt((x_ppp, r_hat))
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


@inline function Base.div(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    # TODO: (issue #25) is there a faster way?
    q, r = divrem(x, y)
    q
end


@inline function Base.rem(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    # TODO: (issue #25) is there a faster way?
    q, r = divrem(x, y)
    r
end


@inline Base.rem(x::MLUInt{N, T}, ::Type{MLUInt{N, T}}) where {N, T} = x

@inline @generated function Base.rem(x::MLUInt{N, T}, ::Type{V}) where {N, T, V <: Integer}
    if bitsizeof(V) <= bitsizeof(T)
        :( x[1] % $V )
    elseif bitsizeof(V) % bitsizeof(T) == 0
        t = bitsizeof(V) ÷ bitsizeof(T)
        expr = :(x[1] % $V)
        for i in 2:t
            expr = :( $expr | ((x[$i] % V) << $((i-1) * bitsizeof(T))) )
        end
        expr
    else
        error("Truncating $(MLUInt{N, T}) to $V is not supported")
    end
end


Base.sizeof(::Type{MLUInt{N, T}}) where {N, T} = sizeof(T) * N


# Accounting for limb bitsizes < 8
bitsizeof(::Type{MLUInt{N, T}}) where {N, T} = bitsizeof(T) * N


Base.Broadcast.broadcastable(x::MLUInt) = (x,)


# The following methods are needed for MLUInt to support mulmod_bitshift()


@inline function double(x::MLUInt{N, T}) where {N, T}
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


@inline function halve(x::MLUInt{N, T}) where {N, T}
    lsb = false
    for i in N:-1:1
        new_lsb = trailing_zeros(x[i]) == 0
        x = setindex(x, (x[i] >> 1) | ifelse(lsb, _msb(T), zero(T)), i)
        lsb = new_lsb
    end
    x
end


# The following methods are needed for MLUInt to support mulmod_widemul()


Base.widen(::Type{MLUInt{N, T}}) where {N, T} = MLUInt{N*2, T}


function Base.widemul(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    n = _most_significant_limb(x) - 1
    t = _most_significant_limb(y) - 1
    w = zero(MLUInt{N*2, T})
    for i in 1:t+1
        c = zero(T)
        hi = zero(T)
        for j in 1:n+1
            hi, lo = mulhilo(x[j], y[i])
            hi, lo = addhilo(hi, lo, c)
            hi, lo = addhilo(hi, lo, w[i + j - 1])
            w = setindex(w, lo, i + j - 1)
            c = hi
        end
        w = setindex(w, c, i + n + 1)
    end
    w
end


@inline mulmod(x::MLUInt{N, T}, y::MLUInt{N, T}, modulus::MLUInt{N, T}) where {N, T} =
    mulmod_widemul(x, y, modulus)


function encompassing_type(tp::Type{MLUInt{N, T}}) where {N, T}
    total_size = cld(N * bitsizeof(T), 8)

    if total_size <= 1
        return UInt8
    elseif total_size <= 2
        return UInt16
    elseif total_size <= 4
        return UInt32
    elseif total_size <= 8
        return UInt64
    elseif total_size <= 16
        return UInt128
    else
        return BigInt
    end
end


Base.eltype(::Type{MLUInt{N, T}}) where {N, T} = T


function Base.:~(x::MLUInt{N, T}) where {N, T}
    for i in 1:N
        x = setindex(x, ~x[i], i)
    end
    x
end


@inline function Base.:&(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    res = zero(MLUInt{N, T})
    for i in 1:N
        res = setindex(res, x[i] & y[i], i)
    end
    res
end


@inline function Base.xor(x::MLUInt{N, T}, y::MLUInt{N, T}) where {N, T}
    res = zero(MLUInt{N, T})
    for i in 1:N
        res = setindex(res, xor(x[i], y[i]), i)
    end
    res
end
