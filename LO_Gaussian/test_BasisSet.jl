include("constants.jl")
include("alias.jl")
include("utils.jl")
include("PGBF.jl")
include("CGBF.jl")
include("overlap.jl")
include("Atoms.jl")
include("sto3g.jl")
include("BasisSet.jl")

function test_main()
    atoms = init_atoms_xyz("H2.xyz")
    println(atoms)

    bas = build_basis(atoms)
    println(bas)
end

test_main()