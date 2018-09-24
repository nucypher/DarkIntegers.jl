using DarkIntegers:
    UInt4, addhilo, mulhilo, mulhilo_widemul, mulhilo_same_type, bitsizeof,
    divhilo, modhilo, divremhilo


@testgroup "single limb arithmetic" begin


function addhilo_ref(x_hi::T, x_lo::T, y::T) where T <: Unsigned
    res = (BigInt(x_lo) + BigInt(x_hi) * BigInt(1) << bitsizeof(T)) + BigInt(y)
    T(res >> bitsizeof(T)), T(res & typemax(T))
end


@testcase "addhilo" begin
    check_function_random(addhilo, addhilo_ref, 3)
end


@testcase tags=[:exhaustive] "addhilo, exhaustive" begin
    check_function_exhaustive(addhilo, addhilo_ref, 3)
end


function mulhilo_ref(x::T, y::T) where T <: Unsigned
    res = BigInt(x) * BigInt(y)
    T(res >> bitsizeof(T)), T(res & typemax(T))
end


mulhilo_funcs = [mulhilo_same_type, mulhilo_widemul, mulhilo]
mulhilo_names = ["same_type", "widemul", "mulhilo"]


@testcase "mulhilo" for func in (mulhilo_funcs => mulhilo_names)
    check_function_random(func, mulhilo_ref, 2)
end


@testcase tags=[:exhaustive] "mulhilo, exhaustive" for func in (mulhilo_funcs => mulhilo_names)
    check_function_exhaustive(func, mulhilo_ref, 2)
end


@testcase tags=[:performance] "mulhilo, performance" for rng in fixed_rng, tp in builtin_uint_types
    x = rand(rng, tp)
    y = rand(rng, tp)

    trial = @benchmark mulhilo_same_type($x, $y)
    @test_result "same_type: " * benchmark_result(trial)

    trial = @benchmark mulhilo_widemul($x, $y)
    @test_result "widemul: " * benchmark_result(trial)
end


division_args_filter(x_hi, x_lo, y) = y > 0 && x_hi < y


divhilo_ref(x_hi::T, x_lo::T, y::T) where T <: Unsigned =
    T(div((BigInt(x_lo) + BigInt(x_hi) * BigInt(1) << bitsizeof(T)), BigInt(y)))


@testcase "divhilo" begin
    check_function_random(divhilo, divhilo_ref, 3, division_args_filter)
end


@testcase tags=[:exhaustive] "divhilo, exhaustive" begin
    check_function_exhaustive(divhilo, divhilo_ref, 3, division_args_filter)
end


modhilo_ref(x_hi::T, x_lo::T, y::T) where T <: Unsigned =
    T(mod((BigInt(x_lo) + BigInt(x_hi) * BigInt(1) << bitsizeof(T)), BigInt(y)))


@testcase "modhilo" begin
    check_function_random(modhilo, modhilo_ref, 3, division_args_filter)
end


@testcase tags=[:exhaustive] "modhilo, exhaustive" begin
    check_function_exhaustive(modhilo, modhilo_ref, 3, division_args_filter)
end


divremhilo_ref(x_hi::T, x_lo::T, y::T) where T <: Unsigned =
    T.(divrem((BigInt(x_lo) + BigInt(x_hi) * BigInt(1) << bitsizeof(T)), BigInt(y)))


@testcase "divremhilo" begin
    check_function_random(divremhilo, divremhilo_ref, 3, division_args_filter)
end


@testcase tags=[:exhaustive] "divremhilo, exhaustive" begin
    check_function_exhaustive(divremhilo, divremhilo_ref, 3, division_args_filter)
end


end
