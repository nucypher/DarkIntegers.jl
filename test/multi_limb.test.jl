using DarkIntegers
using DarkIntegers: UInt4, mulmod_bitshift, mulmod_widemul, _unsafe_convert


@testgroup "multi-limb integers" begin


function check_convert(target, value, unsafe::Bool)
    if unsafe
        _unsafe_convert(target, value)
    else
        convert(target, value)
    end
end


unsafe_fx = [true, false] => ["checked", "unchecked"]


@testcase "conversion, MLUInt to MLUInt" for unsafe in unsafe_fx
    @test check_convert(MLUInt{2, UInt64}, MLUInt{2, UInt64}((123, 456)), unsafe) ==
        MLUInt{2, UInt64}((123, 456))
    @test check_convert(MLUInt{3, UInt64}, MLUInt{2, UInt64}((123, 456)), unsafe) ==
        MLUInt{3, UInt64}((123, 456, 0))
    @test check_convert(MLUInt{1, UInt64}, MLUInt{2, UInt64}((123, 0)), unsafe) ==
        MLUInt{1, UInt64}((123,))

    if !unsafe
        @test_throws InexactError convert(MLUInt{1, UInt64}, MLUInt{2, UInt64}((123, 456)))
    end
end


@testcase "conversion, MLUInt to Signed" for unsafe in unsafe_fx
    @test check_convert(BigInt, MLUInt{2, UInt64}((123, 456)), unsafe) ==
        big(123) + big(456) << 64
    @test check_convert(Int64, MLUInt{2, UInt16}((123, 456)), unsafe) == 123 + 456 << 16
    @test check_convert(Int64, MLUInt{2, UInt32}((0, 0)), unsafe) == 0
    @test check_convert(Int64, MLUInt{2, UInt32}((0xffffffff, 0x7fffffff)), unsafe) ==
        typemax(Int64)

    if !unsafe
        @test_throws InexactError convert(Int64, MLUInt{2, UInt32}((0xffffffff, 0xffffffff)))
        @test_throws InexactError convert(Int32, MLUInt{3, UInt16}((0xffff, 0xffff, 0x0001)))
    end
end


@testcase "conversion, MLUInt to Unsigned" for unsafe in unsafe_fx
    @test check_convert(UInt64, MLUInt{2, UInt16}((123, 456)), unsafe) == 123 + 456 << 16
    @test check_convert(UInt64, MLUInt{2, UInt32}((0, 0)), unsafe) == 0
    @test check_convert(UInt64, MLUInt{2, UInt32}((0xffffffff, 0xffffffff)), unsafe) == typemax(UInt64)

    if !unsafe
        @test_throws InexactError convert(UInt32, MLUInt{3, UInt16}((0xffff, 0xffff, 0x0001)))
    end
end


@testcase "conversion, Integer to MLUInt" for unsafe in unsafe_fx
    @test check_convert(MLUInt{2, UInt32}, 0, unsafe) == MLUInt{2, UInt32}((0, 0))
    @test check_convert(MLUInt{2, UInt32}, typemax(UInt64), unsafe) ==
        MLUInt{2, UInt32}((0xffffffff, 0xffffffff))

    if !unsafe
        @test_throws InexactError convert(MLUInt{2, UInt64}, -1)
        @test_throws InexactError convert(MLUInt{2, UInt16}, typemax(UInt64))
    end
end


@testcase "conversion, Bool to MLUInt" begin
    @test convert(MLUInt{2, UInt32}, false) == zero(MLUInt{2, UInt32})
    @test convert(MLUInt{2, UInt32}, true) == one(MLUInt{2, UInt32})
end


@testcase "comparisons" for op in [<, >, <=, >=]
    check_function_random(MLUInt{3, UInt8}, op, op, 2)
end


@testcase tags=[:exhaustive] "comparisons, exhaustive" for op in [<, >, <=, >=]
    check_function_exhaustive(MLUInt{2, UInt4}, op, op, 2)
end


@testcase "utility functions" begin
    @test zero(MLUInt{2, UInt64}) == MLUInt{2, UInt64}((0, 0))
    @test one(MLUInt{2, UInt64}) == MLUInt{2, UInt64}((1, 0))
    @test oneunit(MLUInt{2, UInt64}) == MLUInt{2, UInt64}((1, 0))
    @test typemin(MLUInt{2, UInt64}) == MLUInt{2, UInt64}((0, 0))
    @test typemax(MLUInt{2, UInt64}) == MLUInt{2, UInt64}((typemax(UInt64), typemax(UInt64)))

    @test zero(MLUInt{2, UInt4}) == MLUInt{2, UInt4}((0, 0))
    @test one(MLUInt{2, UInt4}) == MLUInt{2, UInt4}((1, 0))
    @test oneunit(MLUInt{2, UInt4}) == MLUInt{2, UInt4}((1, 0))
    @test typemin(MLUInt{2, UInt4}) == MLUInt{2, UInt4}((0, 0))
    @test typemax(MLUInt{2, UInt4}) == MLUInt{2, UInt4}((typemax(UInt4), typemax(UInt4)))

    @test iseven(MLUInt{2, UInt64}((2, 1)))
    @test !iseven(MLUInt{2, UInt64}((1, 1)))
    @test !isodd(MLUInt{2, UInt64}((2, 1)))
    @test isodd(MLUInt{2, UInt64}((1, 1)))
    @test iszero(zero(MLUInt{2, UInt64}))
    @test !iszero(one(MLUInt{2, UInt64}))

    @test leading_zeros(MLUInt{3, UInt64}((0, 1234, 0))) == leading_zeros(UInt64(1234)) + 64
    @test trailing_zeros(MLUInt{3, UInt64}((0, 1234, 0))) == trailing_zeros(UInt64(1234)) + 64

    @test eltype(MLUInt{3, UInt64}) == UInt64
    @test eltype(one(MLUInt{3, UInt64})) == UInt64

    @test sizeof(MLUInt{3, UInt64}) == 3 * 8
    @test bitsizeof(MLUInt{3, UInt64}) == 3 * 64
    @test sizeof(MLUInt{3, UInt4}) == 3
    @test bitsizeof(MLUInt{3, UInt4}) == 3 * 4

    @test abs(MLUInt{2, UInt64}((2, 1))) == MLUInt{2, UInt64}((2, 1))
end


@testcase "+" begin
    check_function_random(MLUInt{4, UInt64}, +, +, 2)
end


@testcase tags=[:exhaustive] "+, exhaustive" begin
    check_function_exhaustive(MLUInt{3, UInt4}, +, +, 2)
end


@testcase "-" begin
    check_function_random(MLUInt{4, UInt64}, -, -, 2)
end


@testcase tags=[:exhaustive] "-, exhaustive" begin
    check_function_exhaustive(MLUInt{3, UInt4}, -, -, 2)
end


@testcase "*" begin
    check_function_random(MLUInt{4, UInt64}, *, *, 2)
end


@testcase tags=[:exhaustive] "*, exhaustive" begin
    check_function_exhaustive(MLUInt{3, UInt4}, *, *, 2)
end


@testcase "^" for power in [0, 1, 2, 3, 4, 5]
    func(x) = x^power
    check_function_random(MLUInt{3, UInt8}, func, func, 1)
end


@testcase "div" begin
    check_function_random(
        MLUInt{4, UInt64}, div, div, 2;
        args_filter_predicate=(x, y) -> y > 0)
end


@testcase tags=[:exhaustive] "div, exhaustive" begin
    check_function_exhaustive(
        MLUInt{3, UInt4}, div, div, 2;
        args_filter_predicate=(x, y) -> y > 0)
end


@testcase "divrem" begin
    check_function_random(
        MLUInt{4, UInt64}, divrem, divrem, 2;
        args_filter_predicate=(x, y) -> y > 0)
end


@testcase tags=[:exhaustive] "divrem, exhaustive" begin
    check_function_exhaustive(
        MLUInt{3, UInt4}, divrem, divrem, 2;
        args_filter_predicate=(x, y) -> y > 0)
end


@testcase "promotion" begin
    xi = 12345
    yi = 400

    tp = MLUInt{2, UInt8}
    modulus = 2^bitsizeof(tp)

    x = convert(tp, xi)
    y = convert(tp, yi)

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
    convert(T, mod(T2(x) * T2(y), modulus))
end


modulo_args_filter(x, y, modulus) = modulus >= 2 && x < modulus && y < modulus


@testcase "mulmod" for func in (mulmod_funcs => mulmod_names)
    check_function_random(
        MLUInt{4, UInt8}, func, mulmod_ref, 3;
        args_filter_predicate=modulo_args_filter)
end


@testcase tags=[:performance] "mulmod performance" for rng in fixed_rng
    mptp = MLUInt{2, UInt64}

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
    mptp = MLUInt{2, UInt64}

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
    mptp = MLUInt{8, UInt16}

    function make_args(rng)
        x = rand(rng, UInt128(1):typemax(UInt128))
        y = rand(rng, UInt128(1):typemax(UInt128))
        convert(mptp, x), convert(mptp, y)
    end

    trial = benchmark_distribution(rng, divrem, make_args, 2)
    @test_result "8xUInt16: " * benchmark_distribution_result(trial)
end


@testcase "shift" for rng in fixed_rng, shift_dir in ["left", "right"]
    op = shift_dir == "left" ? (<<) : (>>)
    for limbs in 1:8
        for msl in 1:limbs
            for shift in [0:8:64; 2:8:72]
                tp = MLUInt{limbs, UInt8}
                x = rand(rng, UInt64) & ((UInt64(1) << (msl * 8)) - 1)
                x_mp = convert(tp, x)
                res_mp = op(x_mp, shift)
                res = convert(UInt64, res_mp)
                ref = op(x, shift) & ((one(UInt64) << (limbs * 8)) - 1)
                if res != ref
                    @test_fail "Incorrect result for x=$x_mp and shift=$shift: expected $ref, got $res"
                    return
                end
            end
        end
    end
end


@testcase tags=[:performance] "lshift performance" for rng in fixed_rng
    x = MLUInt{8, UInt32}(reverse((
        0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFE,
        0xBAAEDCE6, 0xAF48A03B, 0xBFD25E8C, 0xD0364141)))

    mptp = MLUInt{8, UInt32}

    function make_args(rng)
        x = rand(rng, one(BigInt):((one(BigInt) << 256) - one(BigInt)))
        y = rand(rng, 1:65)
        convert(mptp, x), y
    end

    trial = benchmark_distribution(rng, >>, make_args, 2)
    @test_result "8xUInt32: " * benchmark_distribution_result(trial)
end


end
