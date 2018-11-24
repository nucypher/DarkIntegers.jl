using Primes: isprime


"""
Select the best available polynomial multiplication function
based on the element type an polynomial length.
"""
@inline @generated function _get_polynomial_mul_function(
        ::Type{T}, ::Val{N}, ::Val{NC}) where {T, N, NC}
    if T <: AbstractRRElem
        m = rr_modulus_simple(T)
        # Regular NTT needs (m - 1) to be a multiple of N,
        # tangent NTT (for negacyclic polynomials) needs it to be a multiple of 2N.
        factor = NC ? 2 * N : N
        if rem(m - 1, factor) == 0 && isprime(m)
            :( ntt_mul )
        else
            :( karatsuba_mul )
        end
    else
        :( karatsuba_mul )
    end
end


"""
Polynomials modulo `x^n-1` (cyclic) or `x^n+1` (negacyclic).
Supports any type that has arithmetic operators defined for it,
including [`RRElem`](@ref) and [`RRElemMontgomery`](@ref).

    Polynomial(coeffs::AbstractArray{T, 1}, negacyclic::Bool) where T

Create a polynomial given the array of coefficients
(the `i`-th coefficient corresponds to the `(i-1)`-th power).
"""
struct Polynomial{T}
    coeffs :: Array{T, 1}
    negacyclic :: Bool
    mul_function :: Function

    @inline function Polynomial(coeffs::AbstractArray{T, 1}, negacyclic::Bool) where T
        len = length(coeffs)
        mul_function = _get_polynomial_mul_function(T, Val(len), Val(negacyclic))
        new{T}(coeffs, negacyclic, mul_function)
    end

    @inline function Polynomial(coeffs::Array{T, 1}, negacyclic::Bool, mul_function) where T
        new{T}(coeffs, negacyclic, mul_function)
    end
end


# It has to match a polynomial with any size and modulus,
# So it can't be a `Polynomial` object.
struct ZeroPolynomial
end


@inline Base.zero(::Polynomial{T}) where T = ZeroPolynomial()


@inline Base.length(p::Polynomial{T}) where T = length(p.coeffs)


@inline function Base.:(==)(p1::Polynomial{T}, p2::Polynomial{T}) where T
    p1.negacyclic == p2.negacyclic && p1.coeffs == p2.coeffs
end


@inline function Base.:*(p1::Polynomial{T}, p2::Polynomial{T}) where T
    @assert p1.negacyclic == p2.negacyclic && length(p1) == length(p2)
    p1.mul_function(p1, p2)
end

@inline Base.:*(p1::Polynomial, p2::ZeroPolynomial) = ZeroPolynomial()
@inline Base.:*(p1::ZeroPolynomial, p2::Polynomial) = ZeroPolynomial()


@inline function Base.:*(p1::Polynomial{T}, p2::Integer) where T
    Polynomial(p1.coeffs .* convert(T, p2), p1.negacyclic, p1.mul_function)
end

@inline function Base.:*(p1::Polynomial{T}, p2::V) where T where V
    Polynomial(p1.coeffs .* convert(T, p2), p1.negacyclic, p1.mul_function)
end


@inline function Base.:*(p1::Integer, p2::Polynomial)
    p2 * p1
end

@inline function Base.:*(p1::Polynomial{T}, p2::T) where T
    Polynomial(p1.coeffs .* p2, p1.negacyclic, p1.mul_function)
end


@inline function Base.div(p1::Polynomial{T}, p2::Integer) where T
   Polynomial(div.(p1.coeffs, p2), p1.negacyclic, p1.mul_function)
end


@inline function Base.:+(p1::Polynomial{T}, p2::Polynomial{T}) where T
    Polynomial(p1.coeffs .+ p2.coeffs, p1.negacyclic, p1.mul_function)
end

@inline function Base.:+(p1::Polynomial{T}, p2::T) where T
    coeffs = copy(p1.coeffs)
    coeffs[1] += p2
    Polynomial(coeffs, p1.negacyclic, p1.mul_function)
end


@inline function Base.:-(p1::Polynomial{T}, p2::Polynomial{T}) where T
    Polynomial(p1.coeffs .- p2.coeffs, p1.negacyclic, p1.mul_function)
end

@inline function Base.:-(p1::Polynomial{T}, p2::Unsigned) where T
    Polynomial(p1.coeffs .- T(p2), p1.negacyclic, p1.mul_function)
end

@inline function Base.:-(p1::Polynomial{T}) where T
    Polynomial(.-p1.coeffs, p1.negacyclic, p1.mul_function)
end

@inline Base.:-(p1::Polynomial, p2::ZeroPolynomial) = p1


"""
    shift_polynomial(p::Polynomial{T}, shift::Integer) where T

Multiply the polynomial by `x^shift`. `shift` can be negative.
"""
@Base.propagate_inbounds @inline function shift_polynomial(
        p::Polynomial{T}, shift::Integer) where T

    if shift == 0
        p
    else
        n = length(p)
        cycle = isodd(fld(shift, n))
        shift = mod(shift, n)

        shift_first = p.negacyclic && (!cycle)
        shift_last = p.negacyclic && cycle

        new_coeffs = similar(p.coeffs)
        coeffs = p.coeffs
        for j in 1:shift
            new_coeffs[j] = shift_first ? -coeffs[n-shift+j] : coeffs[n-shift+j]
        end
        for j in shift+1:n
            new_coeffs[j] = shift_last ? -coeffs[j-shift] : coeffs[j-shift]
        end

        Polynomial(new_coeffs, p.negacyclic, p.mul_function)
    end
end


"""
Multiplies two polynomials of length `l` whose coefficients are stored in arrays
`p1c` and `p2c` starting from locations `p1_s` and `p2_s`, respectively.
The result is placed in the array `res` starting from the location `res_s`.
`res` is assumed to have enough space to store the full result of length `2l`.
"""
@inline function mul_naive(
        l::Int, res::Array{T, 1}, res_s::Int,
        p1c::Array{T, 1}, p1_s::Int, p2c::Array{T, 1}, p2_s::Int) where T
    @simd for j in 1:l
        for k in 1:l
            res[res_s+j+k-2] += p2c[p2_s+k-1] * p1c[p1_s+j-1]
        end
    end
end


"""
Recursive Karatsuba multiplication function.
Multiplies two polynomials of length `full_len` (must be a power of 2)
whose coefficients are stored in arrays
`p1c` and `p2c` starting from locations `p1_s` and `p2_s`, respectively.
The result is placed in the array `res` starting from the location `res_s`.
`res` is assumed to have enough space to store the full result of length `2 * full_len`.
`buf` and `buf2` are temporary buffers with enough space for `2 * full_len` elements
starting from locations `buf_s` and `buf2_s`, respectively.
"""
@inline function _karatsuba_mul(
        full_len::Int, res::Array{T, 1}, res_s::Int,
        p1c::Array{T, 1}, p1_s::Int,
        p2c::Array{T, 1}, p2_s::Int,
        buf::Array{T, 1}, buf_s::Int,
        buf2::Array{T, 1}, buf2_s::Int) where T

    if full_len <= 8
        mul_naive(full_len, res, res_s, p1c, p1_s, p2c, p2_s)
        return
    end

    half_len = div(full_len, 2)

    _karatsuba_mul(
        half_len, res, res_s, p1c, p1_s, p2c, p2_s, buf, buf_s+full_len, buf2, buf2_s)
    _karatsuba_mul(
        half_len, res, res_s+full_len, p1c, p1_s+half_len, p2c, p2_s+half_len, buf, buf_s+full_len,
        buf2, buf2_s)

    z = zero(eltype(p1c))
    @simd for i in buf_s:buf_s+full_len-1
        buf[i] = z
    end

    @simd for i in 1:half_len
        buf2[buf2_s+i-1] = p1c[p1_s+i-1] + p1c[p1_s+half_len+i-1]
        buf2[buf2_s+half_len+i-1] = p2c[p2_s+i-1] + p2c[p2_s+half_len+i-1]
    end

    _karatsuba_mul(half_len, buf, buf_s, buf2, buf2_s, buf2, buf2_s+half_len, buf, buf_s+full_len,
        buf2, buf2_s+full_len)

    @simd for i in 1:full_len
        buf2[buf2_s+i-1] = buf[buf_s+i-1] - res[res_s+i-1] - res[res_s+i+full_len-1]
    end

    @simd for i in 1:full_len
        res[res_s+i+half_len-1] += buf2[buf2_s+i-1]
    end
end


"""
Multiplies two polynomials using Karatsuba algorithm.
Assumes the polynomials have the same length and the same value of the `negacyclic` field.
"""
@Base.propagate_inbounds @inline function karatsuba_mul(
        p1::Polynomial{T}, p2::Polynomial{T}) where T

    full_len = length(p1)
    half_len = div(length(p1), 2)

    z = zero(T)
    r0 = similar(p1.coeffs)
    r1 = similar(p1.coeffs)
    r2 = similar(p1.coeffs)
    buf = similar(p1.coeffs)
    buf2 = similar(p1.coeffs)
    r3 = similar(p1.coeffs)

    @simd for i in 1:full_len
        r0[i] = z
        r1[i] = z
        r2[i] = z
    end

    _karatsuba_mul(half_len, r0, 1, p1.coeffs, 1, p2.coeffs, 1, buf, 1, buf2, 1)
    _karatsuba_mul(half_len, r2, 1, p1.coeffs, half_len+1, p2.coeffs, half_len+1, buf, 1, buf2, 1)

    @simd for i in 1:half_len
        r3[i] = p1.coeffs[i] + p1.coeffs[i+half_len]
        r3[i+half_len] = p2.coeffs[i] + p2.coeffs[i+half_len]
    end
    _karatsuba_mul(half_len, r1, 1, r3, 1, r3, half_len+1, buf, 1, buf2, 1)
    r1 .-= r2 .+ r0

    if p1.negacyclic
        @simd for i in 1:half_len
            r0[i+half_len] += r1[i]
            r0[i] -= r1[i+half_len]
        end
        @simd for i in 1:full_len
            r0[i] -= r2[i]
        end
    else
        @simd for i in 1:half_len
            r0[i+half_len] += r1[i]
            r0[i] += r1[i+half_len]
        end
        @simd for i in 1:full_len
            r0[i] += r2[i]
        end
    end

    Polynomial(r0, p1.negacyclic, p1.mul_function)
end


"""
Multiplies two polynomials using NTT.
Assumes the polynomials have the same length and the same value of the `negacyclic` field.
"""
@inline function ntt_mul(p1::Polynomial{T}, p2::Polynomial{T}) where T
    plan = get_ntt_plan(T, length(p1), p1.negacyclic)
    c1 = copy(p1.coeffs)
    ntt!(plan, c1, false)
    c2 = copy(p2.coeffs)
    ntt!(plan, c2, false)
    c1 .*= c2
    ntt!(plan, c1, true)
    Polynomial(c1, p1.negacyclic, p1.mul_function)
end


@inline function nussbaumer_mul_cyclic(x::Array{T, 1}, y::Array{T, 1}) where T
    # assuming that the polynomial length is power of 2
    n = trailing_zeros(length(x))

    if n == 1
        z1 = x[1] * y[1] + x[2] * y[2]
        z2 = (x[1] + x[2]) * (y[1] + y[2]) - z1
        return [z1, z2]
    end

    m = length(x) >> 1

    xx = copy(x)
    yy = copy(y)

    for k in 1:m
        t = x[m + k]
        xx[m + k] = x[k] - t
        xx[k] = x[k] + t

        t = y[m + k]
        yy[m + k] = y[k] - t
        yy[k] = y[k] + t
    end

    z = similar(x)
    z[1:m] = nussbaumer_mul_cyclic(xx[1:m], yy[1:m])
    z[m+1:2m] = nussbaumer_mul_negacyclic(xx[m+1:2m], yy[m+1:2m])

    for k in 1:m
        t = z[m + k]
        z[m + k] = (z[k] - t)
        z[k] = (z[k] + t)
    end

    # To avoid final rescaling `z` must be divided by 2 here.
    # Instead, in order for the accumulated scaling from cyclic multiplication for the length `2^n`
    # to be the same as for the negacyclic multiplication, we add a multiplication by 2
    # for certain values of `n`.
    if 2^trailing_zeros(n - 1) == n - 1
        z .*= 2
    end

    z
end


@inline function bitreverse32(x::UInt32)
    x = ((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1)
    x = ((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2)
    x = ((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4)
    x = ((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8)
    (x >> 16) | (x << 16)
end


function nussbaumer_mul_negacyclic(x::Array{T, 1}, y::Array{T, 1}) where T

    n = trailing_zeros(length(x))

    if n == 1
        t = x[1] * (y[1] + y[2])
        z1 = t - (x[1] + x[2]) * y[2]
        z2 = t + (x[2] - x[1]) * y[1]
        return [z1, z2]
    end

    m = 1 << (fld(n, 2))
    r = 1 << (cld(n, 2))

    X = Array{T}(undef, 2m, r)
    Y = Array{T}(undef, 2m, r)

    X[1:m,:] .= reshape(x, m, r)
    X[m+1:2m,:] .= reshape(x, m, r)
    Y[1:m,:] .= reshape(y, m, r)
    Y[m+1:2m,:] .= reshape(y, m, r)

    jmax = fld(n, 2)
    for j in fld(n, 2)-1:-1:0
        for st in 1:m
            s = ((st - 1) >> j) << (j + 1)
            t = (st - 1) & ((1 << j) - 1)

            sp = bitreverse32(UInt32(s)) >> (32 - fld(n, 2) - 2 - j) # Remove hardcoding
            sp = sp ÷ 2

            k = sp * r ÷ m

            cycle = isodd(fld(k, r))
            k = mod(k, r)

            shift_first = (!cycle) ? -1 : 1
            shift_last = cycle ? -1 : 1

            st_ = s + t + 1

            Xl = copy(X[st_ + (1 << j),:])
            X[st_+2^j,1:k] .= X[st_,1:k] .+ (-shift_first) .* Xl[r-k+1:r]
            X[st_+2^j,k+1:r] .= X[st_,k+1:r] .+ (-shift_last) .* Xl[1:r-k]
            X[st_,1:k] .= X[st_,1:k] .+ shift_first .* Xl[r-k+1:r]
            X[st_,k+1:r] .= X[st_,k+1:r] .+ shift_last .* Xl[1:r-k]

            Yl = copy(Y[st_ + (1 << j),:])
            Y[st_+2^j,1:k] .= Y[st_,1:k] .+ (-shift_first) .* Yl[r-k+1:r]
            Y[st_+2^j,k+1:r] .= Y[st_,k+1:r] .+ (-shift_last) .* Yl[1:r-k]
            Y[st_,1:k] .= Y[st_,1:k] .+ shift_first .* Yl[r-k+1:r]
            Y[st_,k+1:r] .= Y[st_,k+1:r] .+ shift_last .* Yl[1:r-k]

        end
    end

    Z = similar(X)
    for i = 1:2m
        Z[i,:] .= nussbaumer_mul_negacyclic(X[i,:], Y[i,:])
    end

    for j = 0:fld(n, 2)
        for st in 1:m
            s = ((st - 1) >> j) << (j + 1)
            t = (st - 1) & ((1 << j) - 1)

            sp = bitreverse32(UInt32(s)) >> (32 - fld(n, 2) - 2 - j) # Remove hardcoding
            sp = sp ÷ 2

            k = -(sp * r ÷ m)

            cycle = isodd(fld(k, r))
            k = mod(k, r)

            shift_first = (!cycle) ? -1 : 1
            shift_last = cycle ? -1 : 1

            st_ = s + t + 1

            Zl = copy(Z[st_ + (1 << j),:])
            Z[st_+2^j,1:k] .= shift_first .* (Z[st_,end-k+1:end] .- Zl[end-k+1:end])
            Z[st_+2^j,k+1:end] .= shift_last .* (Z[st_,1:r-k] .- Zl[1:r-k])
            Z[st_,:] .= (Z[st_,:] .+ Zl)

            # To avoid final rescaling `Z` must be divided by 2 here.
        end
    end

    z = similar(x)

    z[1:m] .= Z[1:m,1] .- Z[m+1:2m,r]
    for j = 2:r
        z[m*(j-1)+1:m*j] .= Z[1:m,j] .+ Z[m+1:2m,j-1]
    end

    z
end



function get_scale(n, negacyclic)
    if n < 2
        0
    else
        if negacyclic
            fld(n, 2) + 1 + get_scale(cld(n, 2), true)
        else
            1 + get_scale(n-1, true) + (2^trailing_zeros(n - 1) == n - 1)
        end
    end
end


"""
Polynomial multiplication algorithm by
H.J. Nussbaumer, "Fast Polynomial Transform Algorithms for Digital Convolution"
IEEE Transactions on Acoustics, Speech, and Signal Processing, 28(2), 205–215 (1980)
doi:10.1109/TASSP.1980.1163372

The algorithm in a more clear form is provided in Knuth's TAOCP Vol.2, Exercise 4.6.4-59.
"""
@inline function nussbaumer_mul(p1::Polynomial{T}, p2::Polynomial{T}) where T
    # Assuming that the polynomial length is power of 2
    n = trailing_zeros(length(p1))

    if p1.negacyclic
        coeffs = nussbaumer_mul_negacyclic(p1.coeffs, p2.coeffs)
    else
        coeffs = nussbaumer_mul_cyclic(p1.coeffs, p2.coeffs)
    end

    # Since we did not divide by 2 in the internal multiplication functions
    # (because it may be slow for some residue ring element representations),
    # we need to rescale here.
    scale_exp = get_scale(n, p1.negacyclic)
    scale = 1 << scale_exp
    if T <: AbstractRRElem
        # Assuming here that gcd(scale, modulus) == 1
        # Since scale is a power of 2, it is enough for the modulus to be odd.
        coeffs .*= invmod(scale, rr_modulus_simple(T))
    else
        coeffs .÷= scale
    end

    Polynomial(coeffs, p1.negacyclic)
end

