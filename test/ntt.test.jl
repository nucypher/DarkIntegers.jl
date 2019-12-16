using DarkIntegers: AbstractModUInt, get_root_of_one, get_inverse_coeff, ntt


_known_generator_called = 0


const _known_generator_test_val = convert(MLUInt{2, UInt8}, 257)


function DarkIntegers.known_generator(::Val{_known_generator_test_val})
    global _known_generator_called
    _known_generator_called += 1
    convert(MLUInt{2, UInt8}, 3)
end


@testgroup "NTT" begin


function reference_dft(a::Array{T}, inverse) where T <: AbstractModUInt
    d = length(a)
    r = get_root_of_one(T, d, inverse)

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


@testcase "correctness" for tp in (UInt64, MLUInt{2, UInt32}), mod_tp in (ModUInt, MgModUInt)

    len = 8

    m = UInt64(1) - UInt64(2)^32 # 2^64 - 2^32 + 1 - a prime number
    a = UInt64.(0:len-1) # rand(UInt64, len) .% m

    m_tp = convert(tp, m)

    a_rr = mod_tp{tp, m_tp}.(a)
    af_ref = reference_dft(a_rr, false)
    a_back_ref = reference_dft(af_ref, true)

    af = ntt(a_rr, inverse=false, tangent=false)
    a_back = ntt(af, inverse=true, tangent=false)

    # The results will be in a different order, which is irrelevant for the NTT's main purpose
    # (multiplying polynomials). So we just check that all the elements are the same.
    af_bi = convert.(BigInt, af)
    af_ref_bi = convert.(BigInt, af_ref)
    @test sort(af_bi) == sort(af_ref_bi)

    @test a_back == a_back_ref
    @test a_back == a
end


@testcase "known_generator()" begin
    c1 = _known_generator_called
    g = DarkIntegers.known_generator(Val(_known_generator_test_val))
    c2 = _known_generator_called
    @test c2 - c1 == 1
    @test as_builtin(g) == 3
end


end
