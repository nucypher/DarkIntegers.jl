using DarkIntegers: karatsuba_mul, ntt_mul, nussbaumer_mul


_known_isprime_called = 0


function DarkIntegers.known_isprime(::Val{UInt64(257)})
    global _known_isprime_called
    _known_isprime_called += 1
    true
end


@testgroup "polynomials" begin


cyclicity = ([false, true] => ["cyclic", "negacyclic"])


@testcase "mul_by_monomial()" for negacyclic in cyclicity

    coeffs = 0:9
    m = 21
    rtp = UInt8 # MLUInt{2, UInt8}
    mr = rtp(m)
    mtp = MgModUInt{rtp, mr}

    pm = negacyclic ? negacyclic_modulus : cyclic_modulus

    p = Polynomial(mtp.(coeffs), pm)
    s = negacyclic ? -1 : 1

    @test mul_by_monomial(p, 4) == Polynomial(mtp.([s .* (6:9); 0:5]), pm)
    @test mul_by_monomial(p, 10) == Polynomial(mtp.([s .* (0:9);]), pm)
    @test mul_by_monomial(p, 14) == Polynomial(mtp.([6:9; s .* (0:5)]), pm)
    @test mul_by_monomial(p, 20) == Polynomial(mtp.([0:9;]), pm)
    @test mul_by_monomial(p, 24) == Polynomial(mtp.([s .* (6:9); 0:5]), pm)

    @test mul_by_monomial(p, -4) == Polynomial(mtp.([4:9; s .* (0:3)]), pm)
    @test mul_by_monomial(p, -10) == Polynomial(mtp.([s .* (0:9);]), pm)
    @test mul_by_monomial(p, -14) == Polynomial(mtp.([s .* (4:9); 0:3]), pm)
    @test mul_by_monomial(p, -20) == Polynomial(mtp.([0:9;]), pm)
    @test mul_by_monomial(p, -24) == Polynomial(mtp.([4:9; s .* (0:3)]), pm)

end


@testcase tags=[:performance] "mul_by_monomial(), performance" begin

    m = BigInt(1) << 80 + 1
    p1_ref = BigInt.(rand(UInt128, 64)) .% m

    rtp = MLUInt{2, UInt64}
    mr = convert(rtp, m)
    mtp = MgModUInt{rtp, mr}

    p1 = Polynomial(convert.(mtp, p1_ref), negacyclic_modulus)

    trial = @benchmark mul_by_monomial($p1, 123)
    @test_result benchmark_result(trial)

end


function reference_mul(p1::Polynomial{T, N}, p2::Polynomial{T, N}) where {T, N}
    res = zero(Polynomial{T, N})
    for (j, c) in enumerate(p1.coeffs)
        res = res + mul_by_monomial(p2, j - 1) * c
    end
    res
end


@testcase "multiplication" for negacyclic in cyclicity

    # A prime slightly greater than 2^80, and (modulus - 1) is a multiple of 64 (required for NTT)
    m = BigInt(1440321777275241790996481)
    p1_ref = BigInt.(rand(UInt128, 64)) .% m
    p2_ref = BigInt.(rand(UInt128, 64)) .% m

    rtp = MLUInt{2, UInt64}
    mr = convert(rtp, m)
    mtp = MgModUInt{rtp, mr}

    pm = negacyclic ? negacyclic_modulus : cyclic_modulus
    p1 = Polynomial(convert.(mtp, p1_ref), pm)
    p2 = Polynomial(convert.(mtp, p2_ref), pm)

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
    MgModUInt{UInt8, UInt8(95)},
    # modulus is prime, but modulus-1 is not a multiple of 2*32 - use Karatsuba
    MgModUInt{UInt8, UInt8(97)},
    # modulus is prime, modulus-1 is a multiple of 2*32 - can use NTT for polynomials of size 64
    ModUInt{UInt8, UInt8(193)}
]


@testcase "multiplication choice" for tp in choice_types, negacyclic in cyclicity

    if tp <: AbstractModUInt
        max_val = convert(Int, modulus(tp))
    else
        max_val = 256
    end

    len = 32
    pm = negacyclic ? negacyclic_modulus : cyclic_modulus
    p1 = Polynomial(convert.(tp, mod.(rand(Int, len), max_val)), pm)
    p2 = Polynomial(convert.(tp, mod.(rand(Int, len), max_val)), pm)

    p = p1 * p2
    ref = reference_mul(p1, p2)

    @test ref == p
end


@testcase tags=[:performance] "multiplication, performance" begin

    pm = negacyclic_modulus

    # A prime slightly greater than 2^59, and (modulus - 1) is a multiple of 2^17
    # (required for NTT with sizes up to 2^16 to work)
    m = BigInt(576460752308273153)
    p1_ref = BigInt.(rand(UInt64, 512)) .% m
    p2_ref = BigInt.(rand(UInt64, 512)) .% m

    rtp = UInt64
    mr = convert(rtp, m)
    mtp = MgModUInt{rtp, mr}

    p1 = Polynomial(convert.(mtp, p1_ref), pm)
    p2 = Polynomial(convert.(mtp, p2_ref), pm)

    trial = @benchmark karatsuba_mul($p1, $p2)
    @test_result "Karatsuba: " * benchmark_result(trial)

    trial = @benchmark ntt_mul($p1, $p2)
    @test_result "NTT: " * benchmark_result(trial)

    trial = @benchmark nussbaumer_mul($p1, $p2)
    @test_result "Nussbaumer: " * benchmark_result(trial)

    trial = @benchmark $p1 * $p2
    @test_result "default: " * benchmark_result(trial)

end


@testcase "broadcasting" begin
    coeffs1 = [1, 0, 0, 0]
    coeffs2 = [1, 2, 3, 4]

    p1 = Polynomial(coeffs1, negacyclic_modulus)
    p2 = Polynomial(coeffs2, negacyclic_modulus)

    @test (p1 + 1).coeffs == [2, 0, 0, 0] # check that simple addition works as expected
    @test [p1, p2] .+ 1 == [p1 + 1, p2 + 1] # broadcasting an array of polynomials
    @test broadcast_into_polynomial(+, p1, 1).coeffs == [2, 1, 1, 1] # broadcasted addition

    # Broadcasting of non-polynomials creates a new polynomial with default modulus
    p3 = broadcast_into_polynomial(+, [4, 5, 6, 7], 1)
    @test p3.coeffs == [5, 6, 7, 8]
    @test p3.modulus == cyclic_modulus

    # Modulus mismatch is not allowed
    @test_throws Exception broadcast_into_polynomial(+, p1, p3)

    # Check that the multiplication function is replaced when the coefficient type changes
    @test p1.mul_function == DarkIntegers.karatsuba_mul
    tp = ModUInt{UInt64, UInt64(17)}
    func(x) = convert(tp, x)
    p4 = broadcast_into_polynomial(func, p1)
    @test eltype(p4) == tp
    @test p4.mul_function == DarkIntegers.ntt_mul

    # Broadcast into the old container
    # Converts the values back into the type of `p1`
    broadcast_into_polynomial!(value, p1, p4)
    @test eltype(p1) == Int
    @test p1.mul_function == DarkIntegers.karatsuba_mul

    # Modulus mismatch with the destination is not allowed
    @test_throws Exception broadcast_into_polynomial!(value, p1, p3)
end


@testcase "known_isprime()" begin
    tp = ModUInt{UInt64, UInt64(257)}
    coeffs = tp[1, 0, 0, 0]
    c1 = _known_isprime_called
    p = Polynomial(coeffs, negacyclic_modulus)
    c2 = _known_isprime_called
    @test c2 - c1 == 1
    @test p.mul_function == DarkIntegers.ntt_mul
end


@testcase "type stability" begin
    tp = Polynomial{Int, 4}
    @test Base.promote_op(zero, tp) == tp
    @test Base.promote_op(+, tp, tp) == tp
    @test Base.promote_op(-, tp, tp) == tp
    @test Base.promote_op(-, tp) == tp
    @test Base.promote_op(*, tp, tp) == tp
end


end
