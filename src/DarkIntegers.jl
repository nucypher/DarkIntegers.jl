module DarkIntegers

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

include("mp_number.jl")
export MPNumber

include("residue_ring.jl")
export RRElem
export AbstractRRElem

include("montgomery_reduction.jl")

include("residue_ring_montgomery.jl")
export RRElemMontgomery

include("polynomial.jl")
export Polynomial
export shift_polynomial

include("misc.jl")
export encompassing_type

include("residue_ring_generic.jl")
export rr_base_type
export rr_modulus

include("ntt.jl")

end
