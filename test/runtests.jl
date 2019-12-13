using Jute

include("TestUtils.jl")

include("utils.test.jl")
include("single_limb.test.jl")
include("single_limb_modulo.test.jl")
include("multi_limb.test.jl")
include("multi_limb_signed.test.jl")
include("residue_ring.test.jl")
include("montgomery_reduction.test.jl")
include("residue_ring_montgomery.test.jl")
include("polynomial.test.jl")
include("modification.test.jl")
include("ntt.test.jl")

exit(runtests(options=Dict(:exclude_tags => [:performance, :exhaustive])))
