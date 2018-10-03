using DarkIntegers:
    shift_polynomial,
    reference_mul, fast_reference_mul, karatsuba_mul


@testgroup "polynomials" begin


@testcase "shift_polynomial()" begin

    coeffs = 0:9
    modulus = 21
    rtp = UInt8 # MPNumber{2, UInt8}
    mr = rtp(modulus)
    mtp = RRElemMontgomery{rtp, mr}

    p = Polynomial(mtp, coeffs, true)

    @test shift_polynomial(p, 4) == Polynomial(mtp, [-(6:9); 0:5], true)
    @test shift_polynomial(p, 10) == Polynomial(mtp, [-(0:9);], true)
    @test shift_polynomial(p, 14) == Polynomial(mtp, [6:9; -(0:5)], true)
    @test shift_polynomial(p, 20) == Polynomial(mtp, [0:9;], true)
    @test shift_polynomial(p, 24) == Polynomial(mtp, [-(6:9); 0:5], true)

    @test shift_polynomial(p, -4) == Polynomial(mtp, [4:9; -(0:3)], true)
    @test shift_polynomial(p, -10) == Polynomial(mtp, [-(0:9);], true)
    @test shift_polynomial(p, -14) == Polynomial(mtp, [-(4:9); 0:3], true)
    @test shift_polynomial(p, -20) == Polynomial(mtp, [0:9;], true)
    @test shift_polynomial(p, -24) == Polynomial(mtp, [4:9; -(0:3)], true)

end


@testcase tags=[:performance] "shift_polynomial(), performance" begin

    modulus = BigInt(1) << 80 + 1
    p1_ref = BigInt.(rand(UInt128, 64)) .% modulus

    rtp = MPNumber{2, UInt64}
    mr = rtp(modulus)
    mtp = RRElemMontgomery{rtp, mr}

    p1 = Polynomial(mtp, p1_ref, true)

    trial = @benchmark shift_polynomial($p1, 123)
    @test_result benchmark_result(trial)

end


# TODO: it's completely generic, so the main `shift` implementation
# can be just made to accept anything array-like
function reference_poly_shift(
        p::Array{BigInt, 1}, modulus::BigInt, negacyclic::Bool, shift::Integer)

    if shift == 0
        p
    else
        cycle = mod(fld(shift, length(p)), 2)
        shift = mod(shift, length(p))

        if cycle == 1 && negacyclic
            coeffs = modulus .- p
        else
            coeffs = p
        end

        new_coeffs = circshift(coeffs, shift)

        if negacyclic
            new_coeffs[1:shift] .= modulus .- new_coeffs[1:shift]
        end
        new_coeffs
    end
end


function reference_poly_mul(
        p1::Array{BigInt, 1}, p2::Array{BigInt, 1}, modulus::BigInt, negacyclic::Bool)

    res = zeros(eltype(p1), length(p1))
    for (j, c) in enumerate(p1)
        res = res + reference_poly_shift(p2, modulus, negacyclic, j - 1) * c
    end
    res .% modulus
end


@testcase "multiplication" begin

    negacyclic = true
    modulus = BigInt(1) << 80 + 1
    p1_ref = BigInt.(rand(UInt128, 64)) .% modulus
    p2_ref = BigInt.(rand(UInt128, 64)) .% modulus

    rtp = MPNumber{2, UInt64}
    mr = convert(rtp, modulus)
    mtp = RRElemMontgomery{rtp, mr}

    p1 = Polynomial(mtp, p1_ref, negacyclic)
    p2 = Polynomial(mtp, p2_ref, negacyclic)

    ref = reference_poly_mul(p1_ref, p2_ref, modulus, negacyclic)
    test1 = reference_mul(p1, p2)
    test2 = fast_reference_mul(p1, p2)
    test3 = karatsuba_mul(p1, p2)

    @test ref == convert.(BigInt, test1.coeffs)
    @test ref == convert.(BigInt, test2.coeffs)
    @test ref == convert.(BigInt, test3.coeffs)
end


@testcase tags=[:performance] "multiplication, performance" begin

    negacyclic = true
    modulus = BigInt(1) << 80 + 1
    p1_ref = BigInt.(rand(UInt128, 64)) .% modulus
    p2_ref = BigInt.(rand(UInt128, 64)) .% modulus

    rtp = MPNumber{2, UInt64}
    mr = convert(rtp, modulus)
    mtp = RRElemMontgomery{rtp, mr}

    p1 = Polynomial(mtp, p1_ref, negacyclic)
    p2 = Polynomial(mtp, p2_ref, negacyclic)

    trial = @benchmark fast_reference_mul($p1, $p2)
    @test_result "naive: " * benchmark_result(trial)

    trial = @benchmark karatsuba_mul($p1, $p2)
    @test_result "Karatsuba: " * benchmark_result(trial)

end


end
