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
    @test value(x) != val
    @test value(x) == mod(val, m)

    big_val = Int64(2^50)
    x = MgModUInt{T, m}(big_val)
    @test value(x) == mod(big_val, m)

    big_val = Int64(-2^50)
    x = MgModUInt{T, m}(big_val)
    @test value(x) == mod(big_val, m)


    mp_tp = MLUInt{2, UInt8}
    m = convert(mp_tp, 177)

    mod_tp = MgModUInt{mp_tp, m}

    @test as_builtin(value(mod_tp(convert(mp_tp, 1)))) == 1
    @test as_builtin(value(mod_tp(convert(MLUInt{3, UInt8}, 1)))) == 1
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


@testcase tags=[:performance] "*, performance" for rng in fixed_rng(123)

    modulus = UInt128(2)^80 + 1
    x = rand(rng, UInt128(1):modulus-1)
    y = rand(rng, UInt128(1):modulus-1)

    x_m1 = MgModUInt{UInt128, modulus}(x)
    y_m1 = MgModUInt{UInt128, modulus}(y)
    trial = @benchmark $x_m1 * $y_m1
    @test_result "UInt128: " * benchmark_result(trial)

    mptp1 = MLUInt{2, UInt64}
    x_mp = convert(mptp1, x)
    y_mp = convert(mptp1, y)
    m_mp = convert(mptp1, modulus)
    x_m2 = MgModUInt{mptp1, m_mp}(x_mp)
    y_m2 = MgModUInt{mptp1, m_mp}(y_mp)
    trial = @benchmark $x_m2 * $y_m2
    @test_result "2xUInt64: " * benchmark_result(trial)

    mptp2 = MLUInt{3, UInt32}
    x_mp = convert(mptp2, x)
    y_mp = convert(mptp2, y)
    m_mp = convert(mptp2, modulus)
    x_m3 = MgModUInt{mptp2, m_mp}(x_mp)
    y_m3 = MgModUInt{mptp2, m_mp}(y_mp)
    trial = @benchmark $x_m3 * $y_m3
    @test_result "3xUInt32: " * benchmark_result(trial)

end


@testcase "inv" for rng in fixed_rng(123)
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


@testcase "type stability" begin
    tp = MgModUInt{UInt64, UInt64(5)}
    @test Base.promote_op(zero, tp) == tp
    @test Base.promote_op(one, tp) == tp
    @test Base.promote_op(+, tp, tp) == tp
    @test Base.promote_op(-, tp, tp) == tp
    @test Base.promote_op(-, tp) == tp
    @test Base.promote_op(*, tp, tp) == tp
end


@testcase "rand(type)" for rng in fixed_rng(123)
    m = convert(MLUInt{3, UInt8}, 2^22 - 27)
    tp = MgModUInt{MLUInt{3, UInt8}, m}
    res = rand(rng, tp, 10000)
    @test eltype(res) == tp
    vals = convert.(Int, value.(res))
    @test minimum(vals) <= 10000
    @test maximum(vals) >= 2^22 - 27 - 10000
    @test maximum(vals) < 2^22 - 27
end


@testcase "rand(range)" for rng in fixed_rng(123)

    m = convert(MLUInt{3, UInt8}, 2^22 - 27)
    tp = MgModUInt{MLUInt{3, UInt8}, m}

    start = 10
    len = 256 * 8 - 123

    a = convert(tp, start)
    b = convert(tp, start + len - 1)

    res = rand(rng, a:b, 10000)
    @test eltype(res) == tp

    res_int = convert.(Int, value.(res))
    @test all(res_int .>= start)
    @test all(res_int .<= start + len - 1)
end


@testcase tags=[:performance] "rand() performance" for rng in fixed_rng(123)

    m = convert(MLUInt{4, UInt64}, big(2)^255-19)
    tp = MgModUInt{MLUInt{4, UInt64}, m}

    trial = @benchmark rand($rng, $tp)
    @test_result "type: " * benchmark_result(trial)

    a = rand(rng, tp)
    b = rand(rng, tp)
    r = min(a, b):max(a, b)
    trial = @benchmark rand($rng, $r)
    @test_result "range: " * benchmark_result(trial)
end


end
