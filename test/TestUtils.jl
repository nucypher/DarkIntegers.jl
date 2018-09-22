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
