using DarkIntegers: shift_polynomial, karatsuba_mul, ntt_mul, nussbaumer_mul


@testgroup "polynomials" begin


cyclicity = ([false, true] => ["cyclic", "negacyclic"])


@testcase "shift_polynomial()" for negacyclic in cyclicity

    coeffs = 0:9
    modulus = 21
    rtp = UInt8 # MPNumber{2, UInt8}
    mr = rtp(modulus)
    mtp = RRElemMontgomery{rtp, mr}

    p = Polynomial(mtp.(coeffs), negacyclic)
    s = negacyclic ? -1 : 1

    @test shift_polynomial(p, 4) == Polynomial(mtp.([s .* (6:9); 0:5]), negacyclic)
    @test shift_polynomial(p, 10) == Polynomial(mtp.([s .* (0:9);]), negacyclic)
    @test shift_polynomial(p, 14) == Polynomial(mtp.([6:9; s .* (0:5)]), negacyclic)
    @test shift_polynomial(p, 20) == Polynomial(mtp.([0:9;]), negacyclic)
    @test shift_polynomial(p, 24) == Polynomial(mtp.([s .* (6:9); 0:5]), negacyclic)

    @test shift_polynomial(p, -4) == Polynomial(mtp.([4:9; s .* (0:3)]), negacyclic)
    @test shift_polynomial(p, -10) == Polynomial(mtp.([s .* (0:9);]), negacyclic)
    @test shift_polynomial(p, -14) == Polynomial(mtp.([s .* (4:9); 0:3]), negacyclic)
    @test shift_polynomial(p, -20) == Polynomial(mtp.([0:9;]), negacyclic)
    @test shift_polynomial(p, -24) == Polynomial(mtp.([4:9; s .* (0:3)]), negacyclic)

end


@testcase tags=[:performance] "shift_polynomial(), performance" begin

    modulus = BigInt(1) << 80 + 1
    p1_ref = BigInt.(rand(UInt128, 64)) .% modulus

    rtp = MPNumber{2, UInt64}
    mr = rtp(modulus)
    mtp = RRElemMontgomery{rtp, mr}

    p1 = Polynomial(mtp.(p1_ref), true)

    trial = @benchmark shift_polynomial($p1, 123)
    @test_result benchmark_result(trial)

end


function reference_mul(p1::Polynomial{T}, p2::Polynomial{T}) where T
    res = Polynomial(zeros(T, length(p1)), p1.negacyclic, p1.mul_function)
    for (j, c) in enumerate(p1.coeffs)
        res = res + shift_polynomial(p2, j - 1) * c
    end
    res
end


@testcase "multiplication" for negacyclic in cyclicity

    # A prime slightly greater than 2^80, and (modulus - 1) is a multiple of 64 (required for NTT)
    modulus = BigInt(1440321777275241790996481)
    p1_ref = BigInt.(rand(UInt128, 64)) .% modulus
    p2_ref = BigInt.(rand(UInt128, 64)) .% modulus

    rtp = MPNumber{2, UInt64}
    mr = convert(rtp, modulus)
    mtp = RRElemMontgomery{rtp, mr}

    p1 = Polynomial(mtp.(p1_ref), negacyclic)
    p2 = Polynomial(mtp.(p2_ref), negacyclic)

    ref = reference_mul(p1, p2)
    test1 = karatsuba_mul(p1, p2)
    test2 = ntt_mul(p1, p2)
    test3 = nussbaumer_mul(p1, p2)

    @test ref == test1
    @test ref == test2
    @test ref == test3
end


choice_types = [
    # regular multiplication
    UInt16,
    # modulus is not prime - Karatsuba multiplication
    RRElemMontgomery{UInt8, UInt8(95)},
    # modulus is prime, but modulus-1 is not a multiple of 2*32 - use Karatsuba
    RRElemMontgomery{UInt8, UInt8(97)},
    # modulus is prime, modulus-1 is a multiple of 2*32 - can use NTT for polynomials of size 64
    RRElem{UInt8, UInt8(193)}
]


@testcase "multiplication choice" for tp in choice_types, negacyclic in cyclicity

    if tp <: AbstractRRElem
        max_val = Int(DarkIntegers.rr_modulus(tp))
    else
        max_val = 256
    end

    len = 32
    p1 = Polynomial(tp.(mod.(rand(Int, len), max_val)), negacyclic)
    p2 = Polynomial(tp.(mod.(rand(Int, len), max_val)), negacyclic)

    p = p1 * p2
    ref = reference_mul(p1, p2)

    @test ref == p
end


@testcase tags=[:performance] "multiplication, performance" begin

    negacyclic = true
    # A prime slightly greater than 2^59, and (modulus - 1) is a multiple of 2^17
    # (required for NTT with sizes up to 2^16 to work)
    modulus = BigInt(576460752308273153)
    p1_ref = BigInt.(rand(UInt64, 512)) .% modulus
    p2_ref = BigInt.(rand(UInt64, 512)) .% modulus

    rtp = UInt64
    mr = convert(rtp, modulus)
    mtp = RRElemMontgomery{rtp, mr}

    p1 = Polynomial(mtp.(p1_ref), negacyclic)
    p2 = Polynomial(mtp.(p2_ref), negacyclic)

    trial = @benchmark karatsuba_mul($p1, $p2)
    @test_result "Karatsuba: " * benchmark_result(trial)

    trial = @benchmark ntt_mul($p1, $p2)
    @test_result "NTT: " * benchmark_result(trial)

    trial = @benchmark nussbaumer_mul($p1, $p2)
    @test_result "Nussbaumer: " * benchmark_result(trial)

    trial = @benchmark $p1 * $p2
    @test_result "default: " * benchmark_result(trial)

end


end
