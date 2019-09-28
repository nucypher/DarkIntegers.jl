using DarkIntegers
using DarkIntegers: UInt4, mulmod, mulmod_bitshift, mulmod_widemul


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


@testcase "creation/conversion/promotion" begin
    xi = 12345
    yi = 400

    tp = MPNumber{2, UInt8}
    modulus = 2^bitsizeof(tp)

    x = MPNumber{2, UInt8}(xi)
    y = MPNumber{2, UInt8}(yi)

    @test x + y == mod(xi + yi, modulus)
    @test x + 1 == mod(xi + 1, modulus)
    @test 1 + x == mod(1 + xi, modulus)

    @test x - y == mod(xi - yi, modulus)
    @test x - 1 == mod(xi - 1, modulus)
    @test 1 - x == mod(1 - xi, modulus)

    @test x * y == mod(xi * yi, modulus)
    @test x * 3 == mod(xi * 3, modulus)
    @test 110 * x == mod(110 * xi, modulus)

    @test x ÷ y == mod(xi ÷ yi, modulus)
    @test x ÷ 3 == mod(xi ÷ 3, modulus)
    @test 110 ÷ x == mod(110 ÷ xi, modulus)
end


mulmod_funcs = [mulmod_bitshift, mulmod_widemul]
mulmod_names = ["bitshift", "widemul"]


function mulmod_ref(x::T, y::T, modulus::T) where T <: Unsigned
    T2 = widen(T)
    T(mod(T2(x) * T2(y), modulus))
end


modulo_args_filter(x, y, modulus) = modulus >= 2 && x < modulus && y < modulus


@testcase "mulmod" for func in (mulmod_funcs => mulmod_names)
    check_function_random(
        MPNumber{4, UInt8}, func, mulmod_ref, 3;
        args_filter_predicate=modulo_args_filter)
end


@testcase tags=[:performance] "mulmod performance" for rng in fixed_rng
    mptp = MPNumber{2, UInt64}

    modulus = UInt128(2)^80 + 1
    x = rand(rng, UInt128(1):modulus-1)
    y = rand(rng, UInt128(1):modulus-1)

    trial = @benchmark mulmod($x, $y, $modulus)
    @test_result "UInt128: " * benchmark_result(trial)

    x_mp = mptp(x)
    y_mp = mptp(y)
    m_mp = mptp(modulus)

    trial = @benchmark mulmod_bitshift($x_mp, $y_mp, $m_mp)
    @test_result "2xUInt64, bitshift: " * benchmark_result(trial)

    trial = @benchmark mulmod_widemul($x_mp, $y_mp, $m_mp)
    @test_result "2xUInt64, widemul: " * benchmark_result(trial)
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

    function make_args(rng)
        x = rand(rng, UInt128(1):typemax(UInt128))
        y = rand(rng, UInt128(1):typemax(UInt128))
        convert(mptp, x), convert(mptp, y)
    end

    trial = benchmark_distribution(rng, divrem, make_args, 2)
    @test_result "8xUInt16: " * benchmark_distribution_result(trial)
end


@testcase "lshift" for rng in fixed_rng
    for limbs in 1:8
        for msl in 1:limbs
            for shift in [0:8:64; 2:8:72]
                tp = MPNumber{limbs, UInt8}
                x = rand(rng, UInt64) & ((UInt64(1) << (msl * 8)) - 1)
                x_mp = convert(tp, x)
                res_mp = x_mp >> shift
                res = convert(UInt64, res_mp)
                ref = x >> shift
                if res != ref
                    @test_fail "Incorrect result for x=$x_mp and shift=$shift"
                end
            end
        end
    end
end


@testcase tags=[:performance] "lshift performance" for rng in fixed_rng
    x = MPNumber{8, UInt32}(reverse((
        0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFE,
        0xBAAEDCE6, 0xAF48A03B, 0xBFD25E8C, 0xD0364141)))

    mptp = MPNumber{8, UInt32}

    function make_args(rng)
        x = rand(rng, one(BigInt):((one(BigInt) << 256) - one(BigInt)))
        y = rand(rng, 1:65)
        convert(mptp, x), y
    end

    trial = benchmark_distribution(rng, >>, make_args, 2)
    @test_result "8xUInt32: " * benchmark_distribution_result(trial)
end


end
