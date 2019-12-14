using Jute

include("TestUtils.jl")

include("utils.test.jl")
include("single_limb.test.jl")
include("single_limb_modulo.test.jl")
include("multi_limb.test.jl")
include("multi_limb_signed.test.jl")
include("modulo_int.test.jl")
include("montgomery_reduction.test.jl")
include("modulo_int_montgomery.test.jl")
include("polynomial.test.jl")
include("modification.test.jl")
include("ntt.test.jl")

exit(runtests(options=Dict(:exclude_tags => [:performance, :exhaustive])))
