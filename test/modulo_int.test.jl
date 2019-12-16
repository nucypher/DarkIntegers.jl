using DarkIntegers
using DarkIntegers: _verbatim, mulmod, mulmod_bitshift, mulmod_remhilo, mulmod_widemul


@testgroup "modulo integers" begin


@testcase "construction" begin
    T = UInt16
    m = convert(T, 177)
    val = convert(T, 200)

    # Check that even a value greater than the modulus is not modified
    # when no conversion is requested.

    x = ModUInt(val, m, _verbatim)
    @test convert(T, x) == val

    x = ModUInt{T, m}(val, _verbatim)
    @test convert(T, x) == val

    # Check that a value greater than the modulus is converted correctly

    x = ModUInt{T, m}(val)
    @test convert(T, x) == mod(val, m)

    big_val = Int64(2^50)
    x = ModUInt{T, m}(big_val)
    @test convert(T, x) == convert(T, mod(big_val, m))

    big_val = Int64(-2^50)
    x = ModUInt{T, m}(big_val)
    @test convert(T, x) == convert(T, mod(big_val, m))

end


@testcase "conversion" begin
    mp_tp = MLUInt{2, UInt8}
    m = convert(mp_tp, 177)

    mod_tp = ModUInt{mp_tp, m}

    @test convert(mod_tp, convert(mp_tp, 1)) == mod_tp(1)
    @test convert(mod_tp, convert(MLUInt{3, UInt16}, 1)) == mod_tp(1)
    @test convert(mod_tp, convert(mp_tp, 1)) == mod_tp(1)
    @test convert(Int, mod_tp(1)) == 1
end


@testcase "promotion" begin
    T = MLUInt{2, UInt8}
    modulus = convert(T, 177)
    x = ModUInt{T, modulus}(convert(T, 100))
    y = ModUInt{T, modulus}(convert(T, 90))

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

    # check that big integers are processed correctly:
    # first the modulus is taken, and only then they are converted to the target unsigned type.
    @test ModUInt{T, modulus}(187) == 10
end


@testcase tags=[:performance] "*, performance" for rng in fixed_rng

    modulus = UInt128(2)^80 + 1
    x = rand(rng, UInt128) % modulus
    y = rand(rng, UInt128) % modulus

    trial = @benchmark mulmod($x, $y, $modulus)
    @test_result "UInt128: " * benchmark_result(trial)

    mptp = MLUInt{2, UInt64}
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
    tp = ModUInt{UInt64, modulus}

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
