using DarkIntegers
using DarkIntegers: _verbatim, rr_value


@testgroup "residue ring elements, Montgomery representation" begin


@testcase "construction" begin
    T = UInt16
    modulus = T(177)
    val = T(200)

    # Check that even a value greater than the modulus is not modified
    # when no conversion is requested.

    x = RRElemMontgomery(val, modulus, _verbatim)
    @test rr_value(x) == val

    x = RRElemMontgomery{T, modulus}(val, _verbatim)
    @test rr_value(x) == val

    # Check that a value greater than the modulus is converted correctly

    x = RRElemMontgomery{T, modulus}(val)
    @test rr_value(x) != val
    @test convert(T, x) == mod(val, modulus)

    big_val = Int64(2^50)
    x = RRElemMontgomery{T, modulus}(big_val)
    @test convert(T, x) == mod(big_val, modulus)

    big_val = Int64(-2^50)
    x = RRElemMontgomery{T, modulus}(big_val)
    @test convert(T, x) == mod(big_val, modulus)
end


@testcase "conversion" begin
    mp_tp = MPNumber{2, UInt8}
    modulus = mp_tp(177)

    rr_tp = RRElemMontgomery{mp_tp, modulus}

    @test convert(rr_tp, mp_tp(1)) == rr_tp(1)
    @test convert(rr_tp, MPNumber{3, UInt16}(1)) == rr_tp(1)
    @test convert(rr_tp, mp_tp(1)) == rr_tp(1)
    @test convert(Int, rr_tp(1)) == 1
end


@testcase "promotion" begin
    T = MPNumber{2, UInt8}
    modulus = T(177)
    x = RRElemMontgomery{T, modulus}(100)
    y = RRElemMontgomery{T, modulus}(90)

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
    @test RRElemMontgomery{T, modulus}(-10) == 167
end


@testcase tags=[:performance] "*, performance" for rng in fixed_rng

    modulus = UInt128(2)^80 + 1
    x = rand(rng, UInt128(1):modulus-1)
    y = rand(rng, UInt128(1):modulus-1)

    x_m1 = RRElemMontgomery{UInt128, modulus}(x)
    y_m1 = RRElemMontgomery{UInt128, modulus}(y)
    trial = @benchmark $x_m1 * $y_m1
    @test_result "UInt128: " * benchmark_result(trial)

    mptp1 = MPNumber{2, UInt64}
    x_mp = mptp1(x)
    y_mp = mptp1(y)
    m_mp = mptp1(modulus)
    x_m2 = RRElemMontgomery{mptp1, m_mp}(x_mp)
    y_m2 = RRElemMontgomery{mptp1, m_mp}(y_mp)
    trial = @benchmark $x_m2 * $y_m2
    @test_result "2xUInt64: " * benchmark_result(trial)

    mptp2 = MPNumber{3, UInt32}
    x_mp = mptp2(x)
    y_mp = mptp2(y)
    m_mp = mptp2(modulus)
    x_m3 = RRElemMontgomery{mptp2, m_mp}(x_mp)
    y_m3 = RRElemMontgomery{mptp2, m_mp}(y_mp)
    trial = @benchmark $x_m3 * $y_m3
    @test_result "3xUInt32: " * benchmark_result(trial)

end


@testcase "inv" for rng in fixed_rng
    modulus = UInt64(251)
    tp = RRElemMontgomery{UInt64, modulus}

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
