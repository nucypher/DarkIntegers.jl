using DarkIntegers
using DarkIntegers: _rescale


@testgroup "modification" begin


@inline function rescale_ref(
        new_max::Unsigned, x::T, old_max::T, round_result::Bool) where T <: Unsigned
    xi = convert(BigInt, x)
    mi = convert(BigInt, old_max)
    float_res = xi * new_max / mi
    if round_result
        res = round(BigInt, float_res)
        if res == new_max
            res = zero(BigInt)
        end
    else
        res = floor(BigInt, float_res)
    end
    convert(T, res)
end


@testcase(
"rescale",
for odd_new_max in ([false, true] => ["new max is even", "new max is odd"]),
    round_result in ([false, true] => ["floor", "round"])

    tp = UInt64
    old_max_i = 2^12+1
    old_max = tp(old_max_i)
    new_max = unsigned(odd_new_max ? 2^4+1 : 2^4)

    for i in 0:old_max_i-1

        res = _rescale(new_max, tp(i), old_max, round_result)
        ref = rescale_ref(new_max, tp(i), old_max, round_result)

        if res != ref
            @critical @test_fail(
                "Rescaling $i from $old_max to $new_max: got $res, expected $ref")
        end
    end

end)

end
