using DarkIntegers:
    UInt4, addhilo, mulhilo, mulhilo_widemul, mulhilo_same_type, bitsizeof,
    divhilo, modhilo


@testgroup "single limb arithmetic" begin


function addhilo_ref(x_hi::T, x_lo::T, y::T) where T <: Unsigned
    res = (BigInt(x_lo) + BigInt(x_hi) * BigInt(1) << bitsizeof(T)) + BigInt(y)
    T(res >> bitsizeof(T)), T(res & typemax(T))
end


@testcase "addhilo" begin
    for i in 1:1000
        x_hi = rand(UInt64)
        x_lo = rand(UInt64)
        y = rand(UInt64)
        ref = addhilo_ref(x_hi, x_lo, y)
        test = addhilo(x_hi, x_lo, y)
        if ref != test
            @test_fail "Incorrect result for ($x_hi, $x_lo) + $y: got $test, expected $ref"
            return
        end
    end
end


@testcase tags=[:exhaustive] "addhilo, exhaustive" begin
    tp = UInt4
    tp_range = typemin(tp):typemax(tp)
    for x_hi in tp_range, x_lo in tp_range, y in tp_range
        ref = addhilo_ref(x_hi, x_lo, y)
        test = addhilo(x_hi, x_lo, y)
        if ref != test
            @test_fail "Incorrect result for ($x_hi, $x_lo) + $y: got $test, expected $ref"
            return
        end
    end
end


function mulhilo_ref(x::T, y::T) where T <: Unsigned
    res = BigInt(x) * BigInt(y)
    T(res >> bitsizeof(T)), T(res & typemax(T))
end


mulhilo_funcs = [mulhilo_same_type, mulhilo_widemul, mulhilo]
mulhilo_names = ["same_type", "widemul", "mulhilo"]


@testcase "mulhilo" for func in (mulhilo_funcs => mulhilo_names)
    for i in 1:1000
        x = rand(UInt64)
        y = rand(UInt64)
        ref = mulhilo_ref(x, y)
        test = func(x, y)
        if ref != test
            @test_fail "Incorrect result for $x * $y: got $test, expected $ref"
            return
        end
    end
end


@testcase tags=[:exhaustive] "mulhilo, exhaustive" for func in (mulhilo_funcs => mulhilo_names)
    tp = UInt4
    for x in typemin(tp):typemax(tp), y in typemin(tp):typemax(tp)
        ref = mulhilo_ref(x, y)
        test = func(x, y)
        if ref != test
            @test_fail "Incorrect result for $x * $y: got $test, expected $ref"
            return
        end
    end
end


@testcase tags=[:performance] "mulhilo, performance" for rng in fixed_rng, tp in builtin_uint_types
    x = rand(rng, tp)
    y = rand(rng, tp)

    trial = @benchmark mulhilo_same_type($x, $y)
    @test_result "same_type: " * benchmark_result(trial)

    trial = @benchmark mulhilo_widemul($x, $y)
    @test_result "widemul: " * benchmark_result(trial)
end


function divhilo_ref(x_hi::T, x_lo::T, y::T) where T <: Unsigned
    res = div((BigInt(x_lo) + BigInt(x_hi) * BigInt(1) << bitsizeof(T)), BigInt(y))
    T(res)
end


@testcase "divhilo" begin
    tp = UInt64
    for i in 1:1000
        y = rand(one(tp):typemax(tp))
        x_hi = rand(zero(tp):y-one(tp))
        x_lo = rand(tp)
        ref = divhilo_ref(x_hi, x_lo, y)
        test = divhilo(x_hi, x_lo, y)
        if ref != test
            @test_fail "Incorrect result for div(($x_hi, $x_lo), $y): got $test, expected $ref"
            return
        end
    end
end


@testcase tags=[:exhaustive] "divhilo, exhaustive" begin
    tp = UInt4
    for y in one(tp):typemax(tp)
        for x_hi in zero(tp):y-one(tp), x_lo in typemin(tp):typemax(tp)
            ref = divhilo_ref(x_hi, x_lo, y)
            test = divhilo(x_hi, x_lo, y)
            if ref != test
                @test_fail "Incorrect result for div(($x_hi, $x_lo), $y): got $test, expected $ref"
                return
            end
        end
    end
end


function modhilo_ref(x_hi::T, x_lo::T, y::T) where T <: Unsigned
    res = mod((BigInt(x_lo) + BigInt(x_hi) * BigInt(1) << bitsizeof(T)), BigInt(y))
    T(res)
end


@testcase "modhilo" begin
    tp = UInt64
    for i in 1:1000
        y = rand(one(tp):typemax(tp))
        x_hi = rand(zero(tp):y-one(tp))
        x_lo = rand(tp)
        ref = modhilo_ref(x_hi, x_lo, y)
        test = modhilo(x_hi, x_lo, y)
        if ref != test
            @test_fail "Incorrect result for mod(($x_hi, $x_lo), $y): got $test, expected $ref"
            return
        end
    end
end


@testcase tags=[:exhaustive] "modhilo, exhaustive" begin
    tp = UInt4
    for y in one(tp):typemax(tp)
        for x_hi in zero(tp):y-one(tp), x_lo in typemin(tp):typemax(tp)
            ref = modhilo_ref(x_hi, x_lo, y)
            test = modhilo(x_hi, x_lo, y)
            if ref != test
                @test_fail "Incorrect result for mod(($x_hi, $x_lo), $y): got $test, expected $ref"
                return
            end
        end
    end
end


end
