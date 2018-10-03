using DarkIntegers:
    fft_generic, AbstractRRElem, rr_modulus, ff_inverse, get_root, get_twiddle_base,
    get_inverse_coeff


@testgroup "NTT" begin


function reference_dft(a::Array{T}, inverse) where T <: AbstractRRElem
    d = length(a)
    r = get_twiddle_base(T, d, inverse)

    R = Array{T}(undef, d, d)
    for i in 1:d
        for j in 1:d
            R[i,j] = r^((i-1) * (j-1))
        end
    end

    at = R * a

    if inverse
        get_inverse_coeff(T, d) .* at
    else
        at
    end
end


@testcase "correctness" for tp in (UInt64, MPNumber{2, UInt32}), rr_tp in (RRElem, RRElemMontgomery)
    len = 8

    m = UInt64(1) - UInt64(2)^32 # 2^64 - 2^32 + 1 - a prime number
    a = UInt64.(0:7) # rand(UInt64, len) .% m

    m_tp = tp(m)

    a_rr = rr_tp{tp, m_tp}.(a)
    af_ref = reference_dft(a_rr, false)
    a_back_ref = reference_dft(af_ref, true)

    af = fft_generic(a_rr, false)
    a_back = fft_generic(af, true)

    @test af == af_ref
    @test a_back == a_back_ref
    @test a_back == a
end


end
