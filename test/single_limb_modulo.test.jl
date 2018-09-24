using DarkIntegers: UInt4, addmod, submod, mulmod, mulmod_bitshift, mulmod_modhilo, mulmod_widen


@testgroup "single limb modulo arithmetic" begin


modulo_args_filter(x, y, modulus) = modulus >= 2 && x < modulus && y < modulus


function addmod_ref(x::T, y::T, modulus::T) where T <: Unsigned
    T2 = widen(T)
    T(mod(T2(x) + T2(y), modulus))
end


@testcase "addmod" begin
    check_function_random(
        UInt64, addmod, addmod_ref, 3;
        args_filter_predicate=modulo_args_filter)
end


@testcase tags=[:exhaustive] "addmod, exhaustive" begin
    check_function_exhaustive(
        UInt4, addmod, addmod_ref, 3;
        args_filter_predicate=modulo_args_filter)
end


function submod_ref(x::T, y::T, modulus::T) where T <: Unsigned
    T2 = widen(T)
    T(mod(T2(x) + modulus - y, modulus))
end


@testcase "submod" begin
    check_function_random(
        UInt64, submod, submod_ref, 3;
        args_filter_predicate=modulo_args_filter)
end


@testcase tags=[:exhaustive] "submod, exhaustive" begin
    check_function_exhaustive(
        UInt4, submod, submod_ref, 3;
        args_filter_predicate=modulo_args_filter)
end


mulmod_funcs = [mulmod_bitshift, mulmod_modhilo, mulmod_widen, mulmod]
mulmod_names = ["bitshift", "modhilo", "widen", "mulmod"]


function mulmod_ref(x::T, y::T, modulus::T) where T <: Unsigned
    T2 = widen(T)
    T(mod(T2(x) * T2(y), modulus))
end


@testcase "mulmod" for func in (mulmod_funcs => mulmod_names)
    check_function_random(
        UInt64, func, mulmod_ref, 3;
        args_filter_predicate=modulo_args_filter)
end


@testcase tags=[:exhaustive] "mulmod, exhaustive" for func in (mulmod_funcs => mulmod_names)
    check_function_exhaustive(
        UInt4, func, mulmod_ref, 3;
        args_filter_predicate=modulo_args_filter)
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
