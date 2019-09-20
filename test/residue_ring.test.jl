using DarkIntegers
using DarkIntegers: _verbatim, mulmod, mulmod_bitshift, mulmod_remhilo, mulmod_widemul


@testgroup "residue ring elements" begin


@testcase "construction" begin
    T = UInt16
    modulus = T(177)
    val = T(200)

    # Check that even a value greater than the modulus is not modified
    # when no conversion is requested.

    x = RRElem(val, modulus, _verbatim)
    @test rr_value(x) == val

    x = RRElem{T, modulus}(val, _verbatim)
    @test rr_value(x) == val

    # Check that a value greater than the modulus is converted correctly

    x = RRElem{T, modulus}(val)
    @test rr_value(x) == mod(val, modulus)

    big_val = Int64(2^50)
    x = RRElem{T, modulus}(big_val)
    @test rr_value(x) == convert(T, mod(big_val, modulus))

    big_val = Int64(-2^50)
    x = RRElem{T, modulus}(big_val)
    @test rr_value(x) == convert(T, mod(big_val, modulus))

end


@testcase "conversion" begin
    mp_tp = MPNumber{2, UInt8}
    modulus = mp_tp(177)

    rr_tp = RRElem{mp_tp, modulus}

    @test convert(rr_tp, mp_tp(1)) == rr_tp(1)
    @test convert(rr_tp, MPNumber{3, UInt16}(1)) == rr_tp(1)
    @test convert(rr_tp, mp_tp(1)) == rr_tp(1)
    @test convert(Int, rr_tp(1)) == 1
end


@testcase "promotion" begin
    T = MPNumber{2, UInt8}
    modulus = T(177)
    x = RRElem{T, modulus}(T(100))
    y = RRElem{T, modulus}(T(90))

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

    # check that negative integers are processed correctly:
    # first the modulus is taken, and only then they are converted to the target unsigned type.
    @test RRElem{T, modulus}(-10) == 167
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


@testcase "inv" for rng in fixed_rng
    modulus = UInt64(251)
    tp = RRElem{UInt64, modulus}

    for i in 1:100
        x = rand(rng, 1:250)
        x_tp = convert(tp, x)
        ix = inv(x_tp)
        if x * ix != one(tp)
            @test_fail "Incorrect result for $x: $(convert(Int, ix))"
        end
    end
end


end
