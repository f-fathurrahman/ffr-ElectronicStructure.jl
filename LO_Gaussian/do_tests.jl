include("utils.jl")
test_utils()

include("structures.jl")
include("basfun.jl")
include("overlap.jl")
include("kinetic.jl")
include("nuclear.jl")
#include("pyquante.jl")

test_pgbf()
test_cgbf()
test_kinetic()

test_fgamma()
test_one()
test_na2()

include("two_el.jl")
test_two_terms()
test_coul1()
test_vrr()
test_hrr()

test_geo_basis()

include("rhf.jl")
include("test_molecules.jl")

test_h2()
test_lih()
test_h2o()
