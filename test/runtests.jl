using Jute

include("TestUtils.jl")

include("single_limb.test.jl")
include("single_limb_modulo.test.jl")
include("mp_number.test.jl")
include("residue_ring.test.jl")
include("montgomery_reduction.test.jl")

exit(runtests())
