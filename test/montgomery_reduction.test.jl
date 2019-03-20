using DarkIntegers
using DarkIntegers:
    UInt4, get_montgomery_coeff, mulmod_montgomery, to_montgomery, from_montgomery,
    get_to_montgomery_coeff


# Can't have it inside @testgroup, because then @benchmark does not see it.
function batched_mulmod(res, x, y, modulus, m_prime)
    @inbounds @simd for i in 1:length(x)
        res[i] = mulmod_montgomery(x[i], y[i], modulus, m_prime)
    end
end


@testgroup "Montgomery reduction" begin


limb_type(::Type{MPNumber{N, T}}) where {N, T} = T
limb_type(tp) = tp


@testcase "Montgomery coefficient" for tp in (UInt32, MPNumber{4, UInt8})

    for i in 1:100

        # in Montgomery representation the modulus needs to be odd
        m = rand(UInt128(0):UInt128(2)^(bitsizeof(tp) - 1)-1) * 2 + 1

        m_tp = tp(m)
        m_prime_tp = get_montgomery_coeff(m_tp)

        expected_tp = limb_type(tp)
        if !(typeof(m_prime_tp) == expected_tp)
            @test_fail "Expected type $expected_tp, got $(typeof(m_prime_tp))"
        end

        m_prime = convert(UInt128, m_prime_tp)
        b = convert(UInt128, typemax(limb_type(tp))) + 1

        if (m_prime * m) % b != 1
            @test_fail "Incorrect Montgomery coefficient for $m: $m_prime"
        end
    end
end


function mulmod_montgomery_ref(T_bitsize, x::T, y::T, modulus::T) where T
    T2 = widen(T)

    x2 = T2(x)
    y2 = T2(y)
    m2 = T2(modulus)

    R = one(T2) << T_bitsize
    inv_R = invmod(R, T2(modulus))

    # applying moduli separately so that intermediate results fit in T2
    (((x2 * y2) % m2) * inv_R) % m2
end


function mulmod_montgomery_test(x::T, y::T, modulus::T) where T
    mulmod_montgomery(x, y, modulus, get_montgomery_coeff(modulus))
end


function mulmod_montgomery_predicate(x, y, modulus)
    isodd(modulus) && x < modulus && y < modulus
end


@testcase "mulmod_montgomery()" for tp in (UInt32, MPNumber{4, UInt8})
    check_function_random(
        tp, mulmod_montgomery_test, mulmod_montgomery_ref, 3;
        args_filter_predicate=mulmod_montgomery_predicate,
        ref_needs_bitsize=true)
end


@testcase tags=[:exhaustive] "mulmod_montgomery(), exhaustive" for tp in (UInt4, MPNumber{2, UInt4})
    check_function_exhaustive(
        tp, mulmod_montgomery_test, mulmod_montgomery_ref, 3;
        args_filter_predicate=mulmod_montgomery_predicate,
        ref_needs_bitsize=true)
end


@testcase tags=[:performance] "mulmod_montgomery(), performance" for rng in fixed_rng

    # Test batched performance to check how well the operations will be vectorized

    batch = 1000

    modulus = UInt128(2)^80 + 1
    x = rand(rng, UInt128(1):modulus-1, batch)
    y = rand(rng, UInt128(1):modulus-1, batch)
    m_prime = get_montgomery_coeff(modulus)
    res = similar(x)

    trial = @benchmark batched_mulmod($res, $x, $y, $modulus, $m_prime)
    @test_result "UInt128: " * benchmark_result(trial)

    mp2tp = MPNumber{2, UInt64}
    x_mp2 = mp2tp.(x)
    y_mp2 = mp2tp.(y)
    m_mp2 = mp2tp(modulus)
    m_prime_mp2 = get_montgomery_coeff(m_mp2)
    res_mp2 = similar(x_mp2)

    trial = @benchmark batched_mulmod($res_mp2, $x_mp2, $y_mp2, $m_mp2, $m_prime_mp2)
    @test_result "2xUInt64: " * benchmark_result(trial)

    mp3tp = MPNumber{3, UInt32}
    x_mp3 = mp3tp.(x)
    y_mp3 = mp3tp.(y)
    m_mp3 = mp3tp(modulus)
    m_prime_mp3 = get_montgomery_coeff(m_mp3)
    res_mp3 = similar(x_mp3)

    trial = @benchmark batched_mulmod($res_mp3, $x_mp3, $y_mp3, $m_mp3, $m_prime_mp3)
    @test_result "3xUInt32: " * benchmark_result(trial)

end


function to_montgomery_ref(bitsize, x::T, modulus::T) where T
    T2 = widen(T)
    x2 = T2(x)
    m2 = T2(modulus)
    R = T2(1) << bitsize
    T((x2 * R) % m2)
end


function to_montgomery_test(x::T, modulus::T) where T
    to_montgomery(x, modulus, get_montgomery_coeff(modulus), get_to_montgomery_coeff(modulus))
end


function to_from_montgomery_predicate(x, modulus)
    isodd(modulus) && x < modulus
end


@testcase "to_montgomery()" for tp in (UInt32, MPNumber{4, UInt8})
    check_function_random(
        tp, to_montgomery_test, to_montgomery_ref, 2;
        args_filter_predicate=to_from_montgomery_predicate,
        ref_needs_bitsize=true)
end


@testcase tags=[:exhaustive] "to_montgomery(), exhaustive" for tp in (UInt4, MPNumber{2, UInt4})
    check_function_random(
        tp, to_montgomery_test, to_montgomery_ref, 2;
        args_filter_predicate=to_from_montgomery_predicate,
        ref_needs_bitsize=true)
end


@testcase tags=[:performance] "to_montgomery(), performance" for rng in fixed_rng

    modulus = UInt128(2)^80 + 1
    x = rand(rng, UInt128(1):modulus-1)

    mc = get_montgomery_coeff(x)
    c = get_to_montgomery_coeff(x)
    trial = @benchmark to_montgomery($x, $modulus, $mc, $c)
    @test_result "UInt128: " * benchmark_result(trial)

    mptp1 = MPNumber{2, UInt64}
    x_mp1 = mptp1(x)
    m_mp1 = mptp1(modulus)
    mc_mp1 = get_montgomery_coeff(m_mp1)
    c_mp1 = get_to_montgomery_coeff(m_mp1)

    trial = @benchmark to_montgomery($x_mp1, $m_mp1, $mc_mp1, $c_mp1)
    @test_result "2xUInt64: " * benchmark_result(trial)

    mptp2 = MPNumber{3, UInt32}
    x_mp2 = mptp2(x)
    m_mp2 = mptp2(modulus)
    mc_mp2 = get_montgomery_coeff(m_mp2)
    c_mp2 = get_to_montgomery_coeff(m_mp2)

    trial = @benchmark to_montgomery($x_mp2, $m_mp2, $mc_mp2, $c_mp2)
    @test_result "3xUInt32: " * benchmark_result(trial)

end


function from_montgomery_ref(bitsize, x::T, modulus::T) where T
    T2 = widen(T)
    x2 = T2(x)
    m2 = T2(modulus)
    R = T2(1) << bitsize
    iR = invmod(R, m2)
    T((x2 * iR) % m2)
end


function from_montgomery_test(x::T, modulus::T) where T
    from_montgomery(x, modulus, get_montgomery_coeff(modulus))
end


@testcase "from_montgomery()" for tp in (UInt32, MPNumber{4, UInt8})
    check_function_random(
        tp, from_montgomery_test, from_montgomery_ref, 2;
        args_filter_predicate=to_from_montgomery_predicate,
        ref_needs_bitsize=true)
end


@testcase tags=[:exhaustive] "from_montgomery(), exhaustive" for tp in (UInt4, MPNumber{2, UInt4})
    check_function_random(
        tp, from_montgomery_test, from_montgomery_ref, 2;
        args_filter_predicate=to_from_montgomery_predicate,
        ref_needs_bitsize=true)
end


@testcase tags=[:performance] "from_montgomery(), performance" for rng in fixed_rng

    modulus = UInt128(2)^80 + 1
    x = rand(rng, UInt128(1):modulus-1)

    c = get_montgomery_coeff(x)
    trial = @benchmark from_montgomery($x, $modulus, $c)
    @test_result "UInt128: " * benchmark_result(trial)

    mptp1 = MPNumber{2, UInt64}
    x_mp1 = mptp1(x)
    m_mp1 = mptp1(modulus)
    c_mp1 = get_montgomery_coeff(m_mp1)

    trial = @benchmark from_montgomery($x_mp1, $m_mp1, $c_mp1)
    @test_result "2xUInt64: " * benchmark_result(trial)

    mptp2 = MPNumber{3, UInt32}
    x_mp2 = mptp2(x)
    m_mp2 = mptp2(modulus)
    c_mp2 = get_montgomery_coeff(m_mp2)

    trial = @benchmark from_montgomery($x_mp2, $m_mp2, $c_mp2)
    @test_result "3xUInt32: " * benchmark_result(trial)

end


end
