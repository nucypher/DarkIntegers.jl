using DarkIntegers
using DarkIntegers: _unsafe_convert, Int4


@testgroup "multi-limb integers, signed" begin


function check_convert(target, value, unsafe::Bool)
    if unsafe
        _unsafe_convert(target, value)
    else
        convert(target, value)
    end
end


unsafe_fx = [true, false] => ["checked", "unchecked"]


@testcase "conversion, MLInt to MLInt" for unsafe in unsafe_fx
    @test check_convert(MLInt{2, UInt64}, MLInt{2, UInt64}((123, 456)), unsafe) ==
        MLInt{2, UInt64}((123, 456))
    @test check_convert(MLInt{3, UInt64}, MLInt{2, UInt64}((123, 456)), unsafe) ==
        MLInt{3, UInt64}((123, 456, 0))
    @test check_convert(MLInt{1, UInt64}, MLInt{2, UInt64}((123, 0)), unsafe) ==
        MLInt{1, UInt64}((123,))
    # Check that high limbs are filled with `typemax()` when converting a negative number.
    @test check_convert(MLInt{3, UInt16}, MLInt{2, UInt16}((0x4567, 0x8123)), unsafe) ==
        MLInt{3, UInt16}((0x4567, 0x8123, 0xffff))

    if !unsafe
        @test_throws InexactError convert(MLInt{1, UInt64}, MLInt{2, UInt64}((123, 456)))
    end
end


@testcase "conversion, MLInt to Signed" for unsafe in unsafe_fx
    @test check_convert(BigInt, MLInt{2, UInt64}((123, 456)), unsafe) ==
        big(123) + big(456) << 64
    # Test conversion of a negative number
    @test check_convert(BigInt, MLInt{2, UInt16}((0x4567, 0x8123)), unsafe) ==
        big(signed(0x81234567))
    @test check_convert(Int64, MLInt{2, UInt16}((123, 456)), unsafe) == 123 + 456 << 16
    @test check_convert(Int64, MLInt{2, UInt32}((0, 0)), unsafe) == 0
    @test check_convert(Int64, MLInt{2, UInt32}((0xffffffff, 0x7fffffff)), unsafe) == typemax(Int64)
    @test check_convert(Int64, MLInt{2, UInt32}((0x00000000, 0x80000000)), unsafe) == typemin(Int64)

    if !unsafe
        @test_throws InexactError convert(Int32, MLInt{3, UInt16}((0xffff, 0xffff, 0x0001)))
        @test_throws InexactError convert(Int32, MLInt{2, UInt32}((0x00000000, 0x80000000)))
    end
end


@testcase "conversion, MLInt to Unsigned" for unsafe in unsafe_fx
    @test check_convert(UInt64, MLInt{2, UInt16}((123, 456)), unsafe) == 123 + 456 << 16
    @test check_convert(UInt64, MLInt{2, UInt32}((0, 0)), unsafe) == 0
    @test check_convert(UInt64, MLInt{2, UInt32}((0xffffffff, 0x7fffffff)), unsafe) ==
        unsigned(typemax(Int64))

    if !unsafe
        @test_throws InexactError convert(UInt64, MLInt{2, UInt16}((0x4567, 0x8123)))
        @test_throws InexactError convert(UInt32, MLInt{3, UInt16}((0xffff, 0xffff, 0x0001)))
    end
end


@testcase "conversion, Integer to MLInt" for unsafe in unsafe_fx
    @test check_convert(MLInt{2, UInt32}, big(0), unsafe) == zero(MLInt{2, UInt32})
    @test check_convert(MLInt{2, UInt32}, big(typemin(Int64)), unsafe) == typemin(MLInt{2, UInt32})
    @test check_convert(MLInt{2, UInt32}, big(typemax(Int64)), unsafe) == typemax(MLInt{2, UInt32})
    @test check_convert(MLInt{2, UInt32}, 0, unsafe) == zero(MLInt{2, UInt32})
    @test check_convert(MLInt{2, UInt32}, typemin(Int64), unsafe) == typemin(MLInt{2, UInt32})
    @test check_convert(MLInt{2, UInt32}, typemin(Int64) + 1, unsafe) ==
        MLInt{2, UInt32}((0x00000001, 0x80000000))

    if !unsafe
        @test_throws InexactError convert(MLInt{2, UInt32}, -(big(1) << 63 + 1))
        @test_throws InexactError convert(MLInt{2, UInt32}, big(1) << 63)
    end
end


@testcase "conversion, Bool to MLInt" begin
    @test convert(MLInt{2, UInt32}, false) == zero(MLInt{2, UInt32})
    @test convert(MLInt{2, UInt32}, true) == one(MLInt{2, UInt32})
end


@testcase "comparisons" for op in [<, >, <=, >=]
    check_function_random(MLInt{3, UInt8}, op, op, 2)
end


@testcase tags=[:exhaustive] "comparisons, exhaustive" for op in [<, >, <=, >=]
    check_function_exhaustive(MLInt{2, UInt4}, op, op, 2)
end


@testcase "utility functions" begin
    @test zero(MLInt{2, UInt64}) == MLInt{2, UInt64}((0, 0))
    @test one(MLInt{2, UInt64}) == MLInt{2, UInt64}((1, 0))
    @test oneunit(MLInt{2, UInt64}) == MLInt{2, UInt64}((1, 0))
    @test typemin(MLInt{2, UInt64}) == MLInt{2, UInt64}((zero(UInt64), unsigned(typemin(Int64))))
    @test typemax(MLInt{2, UInt64}) == MLInt{2, UInt64}((typemax(UInt64), unsigned(typemax(Int64))))

    @test zero(MLInt{2, UInt4}) == MLInt{2, UInt4}((0, 0))
    @test one(MLInt{2, UInt4}) == MLInt{2, UInt4}((1, 0))
    @test oneunit(MLInt{2, UInt4}) == MLInt{2, UInt4}((1, 0))
    @test typemin(MLInt{2, UInt4}) == MLInt{2, UInt4}((zero(UInt4), unsigned(typemin(Int4))))
    @test typemax(MLInt{2, UInt4}) == MLInt{2, UInt4}((typemax(UInt4), unsigned(typemax(Int4))))

    @test iseven(MLInt{2, UInt64}((2, 1)))
    @test !iseven(MLInt{2, UInt64}((1, 1)))
    @test !isodd(MLInt{2, UInt64}((2, 1)))
    @test isodd(MLInt{2, UInt64}((1, 1)))
    @test iszero(zero(MLInt{2, UInt64}))
    @test !iszero(one(MLInt{2, UInt64}))

    @test leading_zeros(MLInt{3, UInt64}((0, 1234, 0))) == leading_zeros(UInt64(1234)) + 64
    @test trailing_zeros(MLInt{3, UInt64}((0, 1234, 0))) == trailing_zeros(UInt64(1234)) + 64

    @test eltype(MLInt{3, UInt64}) == UInt64
    @test eltype(one(MLInt{3, UInt64})) == UInt64

    @test sizeof(MLInt{3, UInt64}) == 3 * 8
    @test bitsizeof(MLInt{3, UInt64}) == 3 * 64
    @test sizeof(MLInt{3, UInt4}) == 3
    @test bitsizeof(MLInt{3, UInt4}) == 3 * 4

    @test abs(MLInt{2, UInt64}((2, 1))) == MLInt{2, UInt64}((2, 1))
end


@testcase "type stability" begin
    tp = MLInt{4, UInt64}
    @test Base.promote_op(zero, tp) == tp
    @test Base.promote_op(one, tp) == tp
    @test Base.promote_op(+, tp, tp) == tp
    @test Base.promote_op(-, tp, tp) == tp
    @test Base.promote_op(-, tp) == tp
    @test Base.promote_op(*, tp, tp) == tp
end


end
