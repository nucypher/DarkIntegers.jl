using Base.Iterators: product
using Random
using BenchmarkTools: prettytime, prettymemory, @benchmark

using DarkIntegers
using DarkIntegers: UInt4, bitsizeof, encompassing_type


function benchmark_result(trial)
    time_str = prettytime(minimum(trial.times))

    if trial.allocs > 0
        mem_str = prettymemory(trial.memory)
        alloc_str = ", $mem_str ($(trial.allocs) allocs)"
    else
        alloc_str = ""
    end

    time_str * alloc_str
end


const builtin_uint_types = (UInt8, UInt16, UInt32, UInt64, UInt128)


const fixed_rng = @local_fixture begin
    seed = 123
    rng = MersenneTwister(seed)
    @produce rng "seed=$seed"
end


function check_function(tp, test_func, ref_func, args, ref_needs_bitsize)
    tp_args = convert.(tp, args)
    if ref_needs_bitsize
        ref = ref_func(bitsizeof(tp), args...)
    else
        ref = ref_func(args...)
    end
    tp_test = test_func(tp_args...)
    tp_ref = convert.(tp, ref)
    if tp_ref != tp_test
        @test_fail "Incorrect result for $(tuple(tp_args...)): got $tp_test, expected $tp_ref"
        return false
    end
    true
end


default_args_filter_predicate(args...) = true


function check_function_random(
        tp, test_func, ref_func, arity::Int;
        args_filter_predicate=default_args_filter_predicate,
        ref_needs_bitsize::Bool=false)

    ctr = 0
    itp = encompassing_type(tp)
    tmin = convert(itp, typemin(tp))
    tmax = convert(itp, typemax(tp))
    while ctr < 1000
        args = rand(tmin:tmax, arity)
        if !args_filter_predicate(args...)
            continue
        end
        if !check_function(tp, test_func, ref_func, args, ref_needs_bitsize)
            return
        end
        ctr += 1
    end
end


function check_function_exhaustive(
        tp, test_func, ref_func, arity::Int;
        args_filter_predicate=default_args_filter_predicate,
        ref_needs_bitsize::Bool=false)

    itp = encompassing_type(tp)
    tmin = convert(itp, typemin(tp))
    tmax = convert(itp, typemax(tp))
    for args in product([tmin:tmax for i in 1:arity]...)
        if !args_filter_predicate(args...)
            continue
        end
        if !check_function(tp, test_func, ref_func, args, ref_needs_bitsize)
            return false
        end
    end
    true
end
