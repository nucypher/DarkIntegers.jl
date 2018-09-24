using DarkIntegers
using DarkIntegers: UInt4, _sub_mul


@testgroup "multiprecision arithmetic" begin


@testcase "+" begin
    check_function_random(MPNumber{4, UInt64}, +, +, 2)
end


@testcase tags=[:exhaustive] "+, exhaustive" begin
    check_function_exhaustive(MPNumber{3, UInt4}, +, +, 2)
end


@testcase "-" begin
    check_function_random(MPNumber{4, UInt64}, -, -, 2)
end


@testcase tags=[:exhaustive] "-, exhaustive" begin
    check_function_exhaustive(MPNumber{3, UInt4}, -, -, 2)
end


@testcase "*" begin
    check_function_random(MPNumber{4, UInt64}, *, *, 2)
end


@testcase tags=[:exhaustive] "*, exhaustive" begin
    check_function_exhaustive(MPNumber{3, UInt4}, *, *, 2)
end


@testcase "div" begin
    check_function_random(
        MPNumber{4, UInt64}, div, div, 2;
        args_filter_predicate=(x, y) -> y > 0)
end


@testcase tags=[:exhaustive] "div, exhaustive" begin
    check_function_exhaustive(
        MPNumber{3, UInt4}, div, div, 2;
        args_filter_predicate=(x, y) -> y > 0)
end


@testcase "divrem" begin
    check_function_random(
        MPNumber{4, UInt64}, divrem, divrem, 2;
        args_filter_predicate=(x, y) -> y > 0)
end


@testcase tags=[:exhaustive] "divrem, exhaustive" begin
    check_function_exhaustive(
        MPNumber{3, UInt4}, divrem, divrem, 2;
        args_filter_predicate=(x, y) -> y > 0)
end


@testcase tags=[:performance] "divrem performance, single limb divisor" for rng in fixed_rng
    mptp = MPNumber{2, UInt64}

    x = rand(rng, UInt128(1):typemax(UInt128))
    y = rand(rng, UInt128(1):UInt128(typemax(UInt64)))

    trial = @benchmark divrem($x, $y)
    @test_result "UInt128: " * benchmark_result(trial)

    x_mp = convert(mptp, x)
    y_mp = convert(mptp, y)

    trial = @benchmark divrem($x_mp, $y_mp)
    @test_result "2xUInt64: " * benchmark_result(trial)
end


@testcase tags=[:performance] "divrem performance, multiple limb divisor" for rng in fixed_rng
    mptp = MPNumber{8, UInt16}

    x = rand(rng, UInt128(1):typemax(UInt128))
    y = rand(rng, UInt128(1):typemax(UInt128))

    x_mp = convert(mptp, x)
    y_mp = convert(mptp, y)

    trial = @benchmark divrem($x_mp, $y_mp)
    @test_result "8xUInt16: " * benchmark_result(trial)
end


@testcase "_sub_mul" begin

    tp = UInt8
    len = 4
    rmod = 2^(bitsizeof(tp) * len)

    for i in 1:1000

        xr = MPNumber(tuple(rand(tp, len)...))
        yr = rand(tp)
        zr = MPNumber(tuple(rand(tp, len)...))

        x = convert(Int, xr)
        y = Int(yr)
        z = convert(Int, zr)

        resr, carry = _sub_mul(xr, yr, zr)
        test = convert(Int, resr)
        ref = x - y * z

        if ref < 0
            ref_carry = true
            ref = mod(ref, rmod)
        else
            ref_carry = false
        end

        if (test, carry) != (ref, ref_carry)
            @test_fail (
                "Incorrect result for $((x, y, z)): " *
                "got $((test, carry)), expected $((ref, ref_carry))")
            return
        end
    end
end


end
