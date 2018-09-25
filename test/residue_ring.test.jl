using DarkIntegers
using DarkIntegers: mulmod, mulmod_bitshift, mulmod_modhilo, mulmod_widemul


@testgroup "residue ring elements" begin


@testcase "creation/conversion/promotion" begin
    T = MPNumber{2, UInt8}
    modulus = T(177)
    x = RRElem(T(100), modulus)
    y = RRElem(T(90), modulus)

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
    x = rand(rng, UInt128) % modulus
    y = rand(rng, UInt128) % modulus

    trial = @benchmark mulmod($x, $y, $modulus)
    @test_result "UInt128: " * benchmark_result(trial)

    mptp = MPNumber{2, UInt64}
    x_mp = mptp(x)
    y_mp = mptp(y)
    modulus_mp = mptp(modulus)

    trial = @benchmark mulmod_bitshift($x_mp, $y_mp, $modulus_mp)
    @test_result "bitshift: " * benchmark_result(trial)

    trial = @benchmark mulmod_widemul($x_mp, $y_mp, $modulus_mp)
    @test_result "widen: " * benchmark_result(trial)
end


end
