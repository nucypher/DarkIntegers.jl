module DarkIntegers

include("uint4.jl")
include("single_limb.jl")
include("single_limb_modulo.jl")

include("mp_number.jl")
export MPNumber

include("residue_ring.jl")
export RRElem

include("montgomery_reduction.jl")

end
