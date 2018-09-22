using DarkIntegers: UInt4, addmod, submod, mulmod, mulmod_bitshift, mulmod_modhilo, mulmod_widen


@testgroup "single limb modulo arithmetic" begin


@testcase "addmod" begin
    tp = UInt64
    for i in 1:10000
        modulus = rand(tp(2):typemax(tp))
        x = rand(tp) % modulus
        y = rand(tp) % modulus
        ref = mod(BigInt(x) + BigInt(y), BigInt(modulus))
        test = addmod(x, y, modulus)
        if ref != test
            @test_fail "Incorrect result for $x + $y: got $test, expected $ref"
            return
        end
    end
end


@testcase tags=[:exhaustive] "addmod, exhaustive" begin
    tp = UInt4
    for modulus in UInt4(2):typemax(tp)
        for x in zero(UInt4):modulus-one(UInt4), y in zero(UInt4):modulus-one(UInt4)
            ref = mod(BigInt(x) + BigInt(y), BigInt(modulus))
            test = addmod(x, y, modulus)
            if ref != test
                @test_fail "Incorrect result for $x + $y: got $test, expected $ref"
                return
            end
        end
    end
end


@testcase "submod" begin
    tp = UInt64
    for i in 1:10000
        modulus = rand(tp(2):typemax(tp))
        x = rand(tp) % modulus
        y = rand(tp) % modulus
        ref = mod(BigInt(x) - BigInt(y), BigInt(modulus))
        test = submod(x, y, modulus)
        if ref != test
            @test_fail "Incorrect result for $x - $y: got $test, expected $ref"
            return
        end
    end
end


@testcase tags=[:exhaustive] "submod, exhaustive" begin
    tp = UInt4
    for modulus in UInt4(2):typemax(tp)
        for x in zero(UInt4):modulus-one(UInt4), y in zero(UInt4):modulus-one(UInt4)
            ref = mod(BigInt(x) - BigInt(y), BigInt(modulus))
            test = submod(x, y, modulus)
            if ref != test
                @test_fail "Incorrect result for $x - $y: got $test, expected $ref"
                return
            end
        end
    end
end


mulmod_funcs = [mulmod_bitshift, mulmod_modhilo, mulmod_widen, mulmod]
mulmod_names = ["bitshift", "modhilo", "widen", "mulmod"]


@testcase "mulmod" for func in (mulmod_funcs => mulmod_names)
    tp = UInt64
    for i in 1:10000
        modulus = rand(tp(2):typemax(tp))
        x = rand(tp) % modulus
        y = rand(tp) % modulus
        ref = mod(BigInt(x) * BigInt(y), BigInt(modulus))
        test = func(x, y, modulus)
        if ref != test
            @test_fail "Incorrect result for $x * $y: got $test, expected $ref"
            return
        end
    end
end


@testcase tags=[:exhaustive] "mulmod, exhaustive" for func in (mulmod_funcs => mulmod_names)
    tp = UInt4
    for modulus in UInt4(2):typemax(tp)
        for x in zero(UInt4):modulus-one(UInt4), y in zero(UInt4):modulus-one(UInt4)
            ref = mod(BigInt(x) * BigInt(y), BigInt(modulus))
            test = func(x, y, modulus)
            if ref != test
                @test_fail "Incorrect result for $x * $y: got $test, expected $ref"
                return
            end
        end
    end
end


@testcase tags=[:performance] "mulmod, performance" for rng in fixed_rng, tp in builtin_uint_types

    modulus = rand(rng, tp(2):typemax(tp))
    x = rand(rng, tp) % modulus
    y = rand(rng, tp) % modulus

    trial = @benchmark mulmod_bitshift($x, $y, $modulus)
    @test_result "bitshift: " * benchmark_result(trial)

    trial = @benchmark mulmod_modhilo($x, $y, $modulus)
    @test_result "modhilo: " * benchmark_result(trial)

    trial = @benchmark mulmod_widen($x, $y, $modulus)
    @test_result "widen: " * benchmark_result(trial)
end


end
