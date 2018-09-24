using Base.Iterators: product
using Random
using BenchmarkTools: prettytime, prettymemory, @benchmark


function benchmark_result(trial)
    time_str = prettytime(minimum(trial.times))

    if trial.allocs > 0
        mem_str = prettymemory(trial.memory)
        alloc_str = "$mem_str ($(trial.allocs) allocs)"
    else
        alloc_str = "no allocs"
    end


    "$time_str, $alloc_str"
end


const builtin_uint_types = (UInt8, UInt16, UInt32, UInt64, UInt128)


const fixed_rng = @local_fixture begin
    seed = 123
    rng = MersenneTwister(seed)
    @produce rng "seed=$seed"
end


function check_function(test_func, ref_func, args)
    ref = ref_func(args...)
    test = test_func(args...)
    if ref != test
        @test_fail "Incorrect result for $(tuple(args...)): got $test, expected $ref"
        return false
    end
    true
end


function check_function_random(test_func, ref_func, arity::Int, args_filter_predicate=nothing)
    if args_filter_predicate === nothing
        args_filter_predicate = (args...) -> true
    end
    ctr = 0
    while ctr < 1000
        args = rand(UInt64, arity)
        if !args_filter_predicate(args...)
            continue
        end
        if !check_function(test_func, ref_func, args)
            return
        end
        ctr += 1
    end
end


function check_function_exhaustive(test_func, ref_func, arity::Int, args_filter_predicate=nothing)
    if args_filter_predicate === nothing
        args_filter_predicate = (args...) -> true
    end
    tp_range = zero(UInt4):typemax(UInt4)
    for args in product([tp_range for i in 1:arity]...)
        if !args_filter_predicate(args...)
            continue
        end
        if !check_function(test_func, ref_func, args)
            return false
        end
    end
    true
end
