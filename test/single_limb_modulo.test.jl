using DarkIntegers: UInt4, addmod, submod, mulmod, mulmod_bitshift, mulmod_modhilo, mulmod_widen


@testgroup "single limb modulo arithmetic" begin


addmod_ref(x, y, modulus) = mod(BigInt(x) + BigInt(y), BigInt(modulus))
addmod_args_filter(x, y, modulus) = modulus >= 2 && x < modulus && y < modulus


@testcase "addmod" begin
    check_function_random(addmod, addmod_ref, 3, addmod_args_filter)
end


@testcase tags=[:exhaustive] "addmod, exhaustive" begin
    check_function_exhaustive(addmod, addmod_ref, 3, addmod_args_filter)
end


submod_ref(x, y, modulus) = mod(BigInt(x) - BigInt(y), BigInt(modulus))
submod_args_filter(x, y, modulus) = modulus >= 2 && x < modulus && y < modulus


@testcase "submod" begin
    check_function_random(submod, submod_ref, 3, submod_args_filter)
end


@testcase tags=[:exhaustive] "submod, exhaustive" begin
    check_function_exhaustive(submod, submod_ref, 3, submod_args_filter)
end


mulmod_funcs = [mulmod_bitshift, mulmod_modhilo, mulmod_widen, mulmod]
mulmod_names = ["bitshift", "modhilo", "widen", "mulmod"]


mulmod_ref(x, y, modulus) = mod(BigInt(x) * BigInt(y), BigInt(modulus))
mulmod_args_filter(x, y, modulus) = modulus >= 2 && x < modulus && y < modulus


@testcase "mulmod" for func in (mulmod_funcs => mulmod_names)
    check_function_random(func, mulmod_ref, 3, mulmod_args_filter)
end


@testcase tags=[:exhaustive] "mulmod, exhaustive" for func in (mulmod_funcs => mulmod_names)
    check_function_exhaustive(func, mulmod_ref, 3, mulmod_args_filter)
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
