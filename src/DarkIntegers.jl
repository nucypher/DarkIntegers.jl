module DarkIntegers

using Base: setindex
using Primes: factor, isprime


include("utils.jl")
export bitsizeof
export log_bitsizeof
export num_bits
export encompassing_type
export as_builtin

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
export mul_by_monomial
export cyclic_modulus
export negacyclic_modulus
export known_isprime
export broadcast_into_polynomial
export broadcast_into_polynomial!
export with_modulus

include("modulo_int_helpers.jl")
export value
export raw_value
export modulus

export change_length

include("ntt.jl")
export ntt
export known_generator

include("convolution_kernels.jl")
include("nussbaumer.jl")

end
