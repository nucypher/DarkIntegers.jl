using DarkIntegers:
    UInt4, addhilo, mulhilo, mulhilo_widemul, mulhilo_same_type, bitsizeof,
    divhilo, remhilo, divremhilo, divremhilo_widen, divremhilo_same_type


@testgroup "single limb arithmetic" begin


function addhilo_ref(bitsize, x_hi::T, x_lo::T, y::T) where T <: Unsigned
    T2 = widen(T)
    res = (T2(x_lo) + T2(x_hi) << bitsize) + T2(y)
    T(res >> bitsize), T(res & (one(T) << bitsize - 1))
end


@testcase "addhilo" begin
    check_function_random(UInt64, addhilo, addhilo_ref, 3; ref_needs_bitsize=true)
end


@testcase tags=[:exhaustive] "addhilo, exhaustive" begin
    check_function_exhaustive(UInt4, addhilo, addhilo_ref, 3; ref_needs_bitsize=true)
end


function mulhilo_ref(bitsize, x::T, y::T) where T <: Unsigned
    T2 = widen(T)
    res = T2(x) * T2(y)
    T(res >> bitsize), T(res & (one(T) << bitsize - 1))
end


mulhilo_funcs = [mulhilo_same_type, mulhilo_widemul, mulhilo]
mulhilo_names = ["same_type", "widemul", "mulhilo"]


@testcase "mulhilo" for func in (mulhilo_funcs => mulhilo_names)
    check_function_random(UInt64, func, mulhilo_ref, 2; ref_needs_bitsize=true)
end


@testcase tags=[:exhaustive] "mulhilo, exhaustive" for func in (mulhilo_funcs => mulhilo_names)
    check_function_exhaustive(UInt4, func, mulhilo_ref, 2; ref_needs_bitsize=true)
end


@testcase tags=[:performance] "mulhilo, performance" for rng in fixed_rng, tp in builtin_uint_types
    x = rand(rng, tp)
    y = rand(rng, tp)

    trial = @benchmark mulhilo_same_type($x, $y)
    @test_result "same_type: " * benchmark_result(trial)

    trial = @benchmark mulhilo_widemul($x, $y)
    @test_result "widemul: " * benchmark_result(trial)
end


division_args_filter(x_hi, x_lo, y) = y > 0


function divremhilo_ref(bitsize, x_hi::T, x_lo::T, y::T) where T <: Unsigned
    T2 = widen(T)
    x = T2(x_lo) + T2(x_hi) << bitsize
    q, r = divrem(x, y)
    T(q & typemax(T)), T(r), q >= (one(T2) << bitsize)
end


function divhilo_ref(bitsize, x_hi::T, x_lo::T, y::T) where T <: Unsigned
    q, r, o = divremhilo_ref(bitsize, x_hi, x_lo, y)
    q, o
end


function remhilo_ref(bitsize, x_hi::T, x_lo::T, y::T) where T <: Unsigned
    q, r, o = divremhilo_ref(bitsize, x_hi, x_lo, y)
    r
end


divremhilo_funcs = [divremhilo_same_type, divremhilo_widen, divremhilo]
divremhilo_names = ["same_type", "widen", "divremhilo"]


@testcase "divhilo" begin
    check_function_random(
        UInt64, divhilo, divhilo_ref, 3;
        args_filter_predicate=division_args_filter, ref_needs_bitsize=true)
end


@testcase tags=[:exhaustive] "divhilo, exhaustive" begin
    check_function_exhaustive(
        UInt4, divhilo, divhilo_ref, 3;
        args_filter_predicate=division_args_filter, ref_needs_bitsize=true)
end


@testcase "remhilo" begin
    check_function_random(
        UInt64, remhilo, remhilo_ref, 3;
        args_filter_predicate=division_args_filter, ref_needs_bitsize=true)
end


@testcase tags=[:exhaustive] "remhilo, exhaustive" begin
    check_function_exhaustive(
        UInt4, remhilo, remhilo_ref, 3;
        args_filter_predicate=division_args_filter, ref_needs_bitsize=true)
end


@testcase "divremhilo" for func in (divremhilo_funcs => divremhilo_names)
    check_function_random(
        UInt64, func, divremhilo_ref, 3;
        args_filter_predicate=division_args_filter, ref_needs_bitsize=true)
end


@testcase tags=[:exhaustive] "divremhilo, exhaustive" for func in (divremhilo_funcs => divremhilo_names)
    check_function_exhaustive(
        UInt4, func, divremhilo_ref, 3;
        args_filter_predicate=division_args_filter, ref_needs_bitsize=true)
end


@testcase tags=[:performance] "divremhilo, performance" for rng in fixed_rng, tp in builtin_uint_types
    x_hi = rand(rng, zero(tp):typemax(tp))
    x_lo = rand(rng, zero(tp):typemax(tp))
    y = rand(rng, one(tp):typemax(tp))

    trial = @benchmark divremhilo_same_type($x_hi, $x_lo, $y)
    @test_result "same_type: " * benchmark_result(trial)

    trial = @benchmark divremhilo_widen($x_hi, $x_lo, $y)
    @test_result "widen: " * benchmark_result(trial)

end


end
