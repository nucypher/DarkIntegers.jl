using DarkIntegers
using DarkIntegers: _unsafe_convert


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


end
