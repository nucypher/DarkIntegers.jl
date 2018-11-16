using DarkIntegers
using DarkIntegers: _change_modulus_proportional


@testgroup "modification" begin


@inline function change_modulus_proportional_ref(
        new_modulus::Unsigned, x::T, old_modulus::T) where T <: Unsigned
    xi = convert(BigInt, x)
    mi = convert(BigInt, old_modulus)
    convert(T, round(BigInt, xi * new_modulus / mi))
end


@testcase(
"change_modulus_proportional",
for odd_new_modulus in ([false, true] => ["new modulus is even", "new modulus is odd"])

    tp = UInt64
    old_modulus_i = 2^12+1
    old_modulus = tp(old_modulus_i)
    new_modulus = unsigned(odd_new_modulus ? 2^4+1 : 2^4)

    for i in 0:old_modulus_i-1

        res = _change_modulus_proportional(new_modulus, tp(i), old_modulus)
        ref = change_modulus_proportional_ref(new_modulus, tp(i), old_modulus)

        if res != ref
            @critical @test_fail(
                "Reducing $i from modulus $old_modulus to $new_modulus: got $res, expected $ref")
        end
    end

end)

end
