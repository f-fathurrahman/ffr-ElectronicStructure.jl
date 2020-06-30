include("LatticeVars.jl")
include("AtomicVars.jl")

include("r3frac.jl")
include("r3mv.jl")

function create_lattice_vars()
    LatVecs = zeros(3,3)
    LatVecs[1,:] = [5.0, 5.0, 0.0]
    LatVecs[2,:] = [5.0, 0.0, 5.0]
    LatVecs[3,:] = [0.0, 5.0, 5.0]
    lattice_vars = LatticeVars( LatVecs )
    return lattice_vars
end

function test_AtomicVars()

    lattice_vars = create_lattice_vars()

    maxatoms = 200
    maxspecies = 8
    atposl = zeros(3,maxatoms,maxspecies)
    
    nspecies = 2
    natoms = [1,1]

    # species 1, atom 1
    atposl[:,1,1] = [0.1, 0.1, 0.1]
    # species 2, atom 1
    atposl[:,1,2] = [0.25, 0.25, 0.25]

    atomic_vars = AtomicVars(nspecies, natoms, atposl, lattice_vars)
end

test_AtomicVars()
