using DarkIntegers


@testgroup "residue ring elements, Montgomery representation" begin


@testcase "creation/conversion/promotion" begin
    T = MPNumber{2, UInt8}
    modulus = T(177)
    x = RRElemMontgomery(T(100), modulus)
    y = RRElemMontgomery(T(90), modulus)

    @test x + y == 13
    @test x + 1 == 101
    @test 1 + x == 101

    @test x - y == 10
    @test x - 1 == 99
    @test 1 - x == 78

    @test x * y == 150
    @test x * 3 == 123
    @test 110 * x == 26

    @test x ÷ y == 1
    @test x ÷ 3 == 33
    @test 110 ÷ x == 1

end


@testcase tags=[:performance] "*, performance" for rng in fixed_rng

    modulus = UInt128(2)^80 + 1
    x = rand(rng, UInt128(1):modulus-1)
    y = rand(rng, UInt128(1):modulus-1)

    x_m1 = RRElemMontgomery(x, modulus)
    y_m1 = RRElemMontgomery(y, modulus)
    trial = @benchmark $x_m1 * $y_m1
    @test_result "UInt128: " * benchmark_result(trial)

    mptp = MPNumber{2, UInt64}
    x_mp = mptp(x)
    y_mp = mptp(y)
    m_mp = mptp(modulus)
    x_m2 = RRElemMontgomery(x_mp, m_mp)
    y_m2 = RRElemMontgomery(y_mp, m_mp)
    trial = @benchmark $x_m2 * $y_m2
    @test_result "2xUInt64: " * benchmark_result(trial)

    mptp = MPNumber{3, UInt32}
    x_mp = mptp(x)
    y_mp = mptp(y)
    m_mp = mptp(modulus)
    x_m3 = RRElemMontgomery(x_mp, m_mp)
    y_m3 = RRElemMontgomery(y_mp, m_mp)
    trial = @benchmark $x_m3 * $y_m3
    @test_result "3xUInt32: " * benchmark_result(trial)

end


end
