using DarkIntegers


@testgroup "Utility functions" begin


@testcase "num_bits(signed)" begin
    @test num_bits(UInt64(0)) == 0
    @test num_bits(UInt64(1)) == 1
    @test num_bits(typemax(UInt64)) == 64
end


@testcase "num_bits(unsigned)" begin
    @test num_bits(Int64(0)) == 0
    @test num_bits(Int64(1)) == 1
    @test num_bits(Int64(-1)) == 1
    @test num_bits(typemax(Int64)) == 63
    @test num_bits(typemin(Int64)) == 64
end


@testcase "num_bits(BigInt)" begin
    @test num_bits(big(0)) == 0
    @test num_bits(big(1)) == 1
    @test num_bits(big(-1)) == 1
    @test num_bits(big(1) << 234) == 235
    @test num_bits(-(big(1) << 234)) == 235
end


end
