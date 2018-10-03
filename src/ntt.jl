using Primes: factor


rr_modulus(::Type{RRElem{T, M}}) where {T, M} = M
rr_modulus(::Type{RRElemMontgomery{T, M}}) where {T, M} = M


function rr_modulus_simple(tp::Type{<:AbstractRRElem})
    m = rr_modulus(tp)
    convert(encompassing_type(typeof(m)), m)
end


function ff_inverse(x::T) where T <: AbstractRRElem
    m = rr_modulus_simple(T)
    x^(m - 2)
end


function get_root(::Type{T}) where T <: AbstractRRElem
    modulus = rr_modulus_simple(T)
    factors = keys(factor(modulus - 1))
    for w in 2:modulus-1
        found = true
        gw = T(w)
        for q in factors
            if gw^(div(modulus - 1, q)) == one(T)
                found = false
                break
            end
        end
        if found
            return gw
        end
    end
    zero(T)
end


function get_twiddle_base(::Type{T}, N, inverse) where T <: AbstractRRElem
    m = rr_modulus_simple(T)
    w = get_root(T) # w^(modulus(T) - 1) = 1
    b = w^div(m - 1, N)
    if inverse
        ff_inverse(b)
    else
        b
    end
end


function get_inverse_coeff(::Type{T}, N) where T <: AbstractRRElem
    # To compare with the DFT
    m = rr_modulus_simple(T)
    T(m - div(m - 1, N))
end


function bitreverse(x, l)
    # Slow, but simple function. Needs some optimized algorithm in practice.
    parse(Int, "0b" * reverse(bitstring(x)[end-l+1:end]))
end


function fft_generic!(cdata, inverse)

    n = length(cdata)

    if n < 2 || n & (n-1) != 0
        error("n must be >=2 and a power of 2")
    end

    logn = round(Int, log2(n))

    for i in 1:n
        j = bitreverse(i-1, logn) + 1
        if j > i
            cdata[j], cdata[i] = cdata[i], cdata[j]
        end
    end

    W = get_twiddle_base(eltype(cdata), n, inverse)

    for stage in 1:logn
        mmax = 2^(stage-1)
        istep = mmax * 2

        for m = 1:mmax
            w = W^((m-1) * 2^(logn - stage))

            for i = m:istep:n
                j = i + mmax

                temp = w * cdata[j]

                cdata[j] = cdata[i] - temp
                cdata[i] += temp
            end
        end
    end

    if inverse
        cdata .= cdata .* get_inverse_coeff(eltype(cdata), n)
    end
end


function fft_generic(cdata, inverse)
    cdata = copy(cdata)
    fft_generic!(cdata, inverse)
    cdata
end


function tangent_ntt(a::Array{T, 1}, inverse::Bool) where T <: AbstractRRElem
    N = length(a)
    w = get_twiddle_base(T, 2 * N, false)
    idx = collect(0:N-1)

    if !inverse
        a = a .* w.^idx
    end

    res = fft_generic(a, inverse)

    if inverse
        res = res .* (w.^(mod.(2 * N .- idx, 2 * N)))
    end

    res
end


function ntt_mul(p1::Polynomial{T}, p2::Polynomial{T}) where T
    @assert p1.negacyclic && p2.negacyclic
    c1_tr = tangent_ntt(p1.coeffs, false)
    c2_tr = tangent_ntt(p2.coeffs, false)
    Polynomial(tangent_ntt(c1_tr .* c2_tr, true), true)
end
