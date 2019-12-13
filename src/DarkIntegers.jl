module DarkIntegers

using Base: setindex
using Primes: factor, isprime


include("utils.jl")
export bitsizeof
export log_bitsizeof
export num_bits
export encompassing_type

include("uint4.jl")

include("single_limb.jl")
export addhilo
export mulhilo
export divremhilo
export divhilo
export remhilo

include("single_limb_modulo.jl")
export addmod
export submod
export mulmod
export powmod

include("multi_limb.jl")
export MLUInt

include("multi_limb_signed.jl")
export MLInt

include("residue_ring.jl")
export RRElem
export AbstractRRElem

include("montgomery_reduction.jl")

include("residue_ring_montgomery.jl")
export RRElemMontgomery

include("polynomial.jl")
export Polynomial
export shift_polynomial

include("modification.jl")
export change_modulus
export change_representation
export change_base_type
export change_length
export rescale

include("ntt.jl")

include("convolution_kernels.jl")
include("nussbaumer.jl")

end
