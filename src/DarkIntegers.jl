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

include("modulo_int.jl")
export ModUInt
export AbstractModUInt

include("montgomery_reduction.jl")

include("modulo_int_montgomery.jl")
export MgModUInt

include("polynomial.jl")
export Polynomial
export ZeroPolynomial
export shift_polynomial
export cyclic_modulus
export negacyclic_modulus

include("modulo_int_helpers.jl")
export raw_value
export raw_value_as_builtin
export value_as_builtin
export modulus
export modulus_as_builtin

export change_modulus
export change_representation
export change_base_type
export change_length
export rescale

include("ntt.jl")

include("convolution_kernels.jl")
include("nussbaumer.jl")

end
