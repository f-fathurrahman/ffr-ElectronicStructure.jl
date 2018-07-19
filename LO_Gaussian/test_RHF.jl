using Printf
using LinearAlgebra
using SpecialFunctions

include("constants.jl")
include("alias.jl")
include("Atoms.jl")
include("utils.jl")
include("PGBF.jl")
include("CGBF.jl")
include("overlap.jl")
include("kinetic.jl")
include("nuclear.jl")
include("coulomb.jl")
include("sto3g.jl")
include("BasisSet.jl")
include("RHF.jl")

function test_main(filename)
    atoms = init_atoms_xyz(filename)
    energy, E, U = RHF(atoms)
end

#test_main("H2_v2.xyz")
#test_main("H2O.xyz")
#test_main("LiH.xyz")
test_main("I2.xyz")