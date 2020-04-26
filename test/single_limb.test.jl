using DarkIntegers:
    UInt4, addhilo, mulhilo, mulhilo_widemul, mulhilo_same_type, bitsizeof,
    divhilo, remhilo, divremhilo, divremhilo_widen, divremhilo_same_type,
    addcarry, subborrow, muladdcarry


@testgroup "single limb arithmetic" begin


function addcarry_ref(x::T, y::T, z::T=zero(T)) where T
    b = big(x) + big(y) + big(z)
    (b >> bitsizeof(T)) % T, b % T
end


@testcase "addcarry" for tp in (UInt64, UInt128), arity in (2, 3)
    check_function_random(tp, addcarry, addcarry_ref, arity)
end


@testcase tags=[:performance] "addcarry, performance" for rng in fixed_rng, tp in builtin_uint_types
    trial = @benchmark addcarry(x, y, z) setup=begin
        x = rand($rng, $tp)
        y = rand($rng, $tp)
        z = rand($rng, $tp)
    end
    @test_result benchmark_result(trial)
end


function subborrow_ref(x::T, y::T, borrow::T) where T
    b = big(x) - big(y) - (big(borrow >> (bitsizeof(T) - 1)))
    if b >= 0
        zero(T), b % T
    else
        bb = b + (one(BigInt) << bitsizeof(T))
        typemax(T), bb % T
    end
end


@testcase "subborrow" for tp in (UInt64, UInt128)
    # Two checks to test both possible values of `borrow`
    check_function_random(
        tp,
        (x, y) -> subborrow(x, y, zero(UInt64)),
        (x, y) -> subborrow_ref(x, y, zero(UInt64)),
        2)
    check_function_random(
        tp,
        (x, y) -> subborrow(x, y, typemax(UInt64)),
        (x, y) -> subborrow_ref(x, y, typemax(UInt64)),
        2)

    check_function_random(tp, subborrow(x, y), subborrow_ref(x, y), 2)
end


@testcase tags=[:performance] "subborrow, performance" for rng in fixed_rng, tp in builtin_uint_types
    trial = @benchmark subborrow(x, y, borrow) setup=begin
        x = rand($rng, $tp)
        y = rand($rng, $tp)
        borrow = rand($rng, Bool) ? typemax($tp) : zero($tp)
    end
    @test_result benchmark_result(trial)
end


function muladdcarry_ref(x::T, y::T, z::T, w::T) where T
    b = big(x) + big(y) * big(z) + big(w)
    (b >> bitsizeof(T)) % T, b % T
end


@testcase "muladdcarry" for tp in (UInt64, UInt128)
    check_function_random(tp, muladdcarry, muladdcarry_ref, 4)
end


@testcase tags=[:performance] "muladdcarry, performance" for rng in fixed_rng, tp in builtin_uint_types
    trial = @benchmark muladdcarry(x, y, z, w) setup=begin
        x = rand($rng, $tp)
        y = rand($rng, $tp)
        z = rand($rng, $tp)
        w = rand($rng, $tp)
    end
    @test_result benchmark_result(trial)
end


function addhilo_ref(bitsize, x_hi::T, x_lo::T, y::T) where T <: Unsigned
    T2 = widen(T)
    res = (T2(x_lo) + T2(x_hi) << bitsize) + T2(y)
    convert(T, res >> bitsize), convert(T, res & (one(T) << bitsize - 1))
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
    convert(T, res >> bitsize), convert(T, res & (one(T) << bitsize - 1))
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
    convert(T, q & typemax(T)), convert(T, r), q >= (one(T2) << bitsize)
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
