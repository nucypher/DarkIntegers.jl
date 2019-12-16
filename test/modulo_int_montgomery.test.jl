using DarkIntegers


@testgroup "modulo integers, Montgomery representation" begin


@testcase "construction" begin
    T = UInt16
    m = T(177)
    val = T(200)

    # Check that even a value greater than the modulus is not modified
    # when no conversion is requested.

    x = MgModUInt(val, m, _verbatim)
    @test raw_value(x) == val

    x = MgModUInt{T, m}(val, _verbatim)
    @test raw_value(x) == val

    # Check that a value greater than the modulus is converted correctly

    x = MgModUInt{T, m}(val)
    @test convert(T, x) != val
    @test convert(T, x) == mod(val, m)

    big_val = Int64(2^50)
    x = MgModUInt{T, m}(big_val)
    @test convert(T, x) == mod(big_val, m)

    big_val = Int64(-2^50)
    x = MgModUInt{T, m}(big_val)
    @test convert(T, x) == mod(big_val, m)
end


@testcase "conversion" begin
    mp_tp = MLUInt{2, UInt8}
    m = convert(mp_tp, 177)

    mod_tp = MgModUInt{mp_tp, m}

    @test convert(mod_tp, convert(mp_tp, 1)) == mod_tp(1)
    @test convert(mod_tp, convert(MLUInt{3, UInt8}, 1)) == mod_tp(1)
    @test convert(mod_tp, convert(mp_tp, 1)) == mod_tp(1)
    @test convert(Int, mod_tp(1)) == 1
end


@testcase "promotion" begin
    T = MLUInt{2, UInt8}
    modulus = convert(T, 177)
    x = MgModUInt{T, modulus}(100)
    y = MgModUInt{T, modulus}(90)

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
    @test MgModUInt{T, modulus}(187) == 10
end


@testcase tags=[:performance] "*, performance" for rng in fixed_rng

    modulus = UInt128(2)^80 + 1
    x = rand(rng, UInt128(1):modulus-1)
    y = rand(rng, UInt128(1):modulus-1)

    x_m1 = MgModUInt{UInt128, modulus}(x)
    y_m1 = MgModUInt{UInt128, modulus}(y)
    trial = @benchmark $x_m1 * $y_m1
    @test_result "UInt128: " * benchmark_result(trial)

    mptp1 = MLUInt{2, UInt64}
    x_mp = mptp1(x)
    y_mp = mptp1(y)
    m_mp = mptp1(modulus)
    x_m2 = MgModUInt{mptp1, m_mp}(x_mp)
    y_m2 = MgModUInt{mptp1, m_mp}(y_mp)
    trial = @benchmark $x_m2 * $y_m2
    @test_result "2xUInt64: " * benchmark_result(trial)

    mptp2 = MLUInt{3, UInt32}
    x_mp = mptp2(x)
    y_mp = mptp2(y)
    m_mp = mptp2(modulus)
    x_m3 = MgModUInt{mptp2, m_mp}(x_mp)
    y_m3 = MgModUInt{mptp2, m_mp}(y_mp)
    trial = @benchmark $x_m3 * $y_m3
    @test_result "3xUInt32: " * benchmark_result(trial)

end


@testcase "inv" for rng in fixed_rng
    modulus = UInt64(251)
    tp = MgModUInt{UInt64, modulus}

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
