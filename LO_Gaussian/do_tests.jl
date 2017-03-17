include("utils.jl")
test_utils()

include("structures.jl")
include("basfun.jl")
include("overlap.jl")
include("kinetic.jl")
include("nuclear.jl")

test_pgbf()
test_cgbf()
test_kinetic()

test_fgamma()
