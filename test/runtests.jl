using Jute

include("TestUtils.jl")

include("single_limb.test.jl")
include("single_limb_modulo.test.jl")
include("mp_number.test.jl")

exit(runtests())
