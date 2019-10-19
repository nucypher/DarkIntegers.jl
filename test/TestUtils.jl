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


struct DistributionTrial
    times :: Array{Float64, 1} # in nanoseconds
end


@generated function _benchmark_distribution(
        rng, test_func, make_args, ::Val{nargs}, samples, iterations) where nargs

    # When measuring fast functions, splicing of arguments produced by `make_args()`
    # into `test_func()` introduces a noticeable delay (tens of nanoseconds).
    # This function is parametrized on the number of arguments,
    # generating the splicing in compile time.

    argnames = [gensym("arg") for i in 1:nargs]

    quote
        ts = Array{Float64, 1}(undef, samples)
        for i in 1:samples
            ($(argnames...),) = make_args(rng)
            tmin = typemax(Float64)
            for j in 1:10
                t = @elapsed test_func($(argnames...))
                if t < tmin
                    tmin = t
                end
            end
            ts[i] = tmin
        end
        DistributionTrial(ts * 1e9)
    end

end


# A trampoline for `_benchmark_distribution()` to avoid exposing the `Val(nargs)` argument.
function benchmark_distribution(
        rng, test_func, make_args, nargs::Int; samples=10000, iterations=10)
    _benchmark_distribution(
        rng, test_func, make_args, Val(nargs), samples, iterations)
end


mean(vals::AbstractArray) = sum(vals) / length(vals)


std(vals::AbstractArray) = sqrt(sum((vals .- mean(vals)).^2 / (length(vals) - 1)))


function benchmark_distribution_result(trial::DistributionTrial)
    mean_str = prettytime(mean(trial.times))
    std_str = prettytime(std(trial.times))
    "$mean_str (±$std_str)"
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
        if !args_filter_predicate(tp.(args)...)
            continue
        end
        if !check_function(tp, test_func, ref_func, args, ref_needs_bitsize)
            return false
        end
    end
    true
end
