using Printf

include("INC.jl")

function create_lattice_vars()
    LatVecs = zeros(3,3)
    A = 5.13
    LatVecs[1,:] = [A, A, 0.0]
    LatVecs[2,:] = [A, 0.0, A]
    LatVecs[3,:] = [0.0, A, A]
    lattice_vars = LatticeVars( LatVecs )
    return lattice_vars
end

function create_atomic_vars(lattice_vars)

    maxatoms = 200
    maxspecies = 8
    atposl = zeros(3,maxatoms,maxspecies)
    
    nspecies = 2
    natoms = [1,1]

    # species 1, atom 1
    atposl[:,1,1] = [0.0, 0.0, 0.0]
    # species 2, atom 1
    atposl[:,1,2] = [0.25, 0.25, 0.25]

    atomic_vars = AtomicVars(nspecies, natoms, atposl, lattice_vars)
end

function main()

    latt_vars = create_lattice_vars()
    atm_vars = create_atomic_vars(latt_vars)

    Nspecies = 2
    atsp_vars = AtomicSpeciesVars(Nspecies)
    mtr_vars = MuffinTinRadialVars(Nspecies)
    apwlo_vars = APWLOVars(Nspecies, mtr_vars.lmaxapw)

    readspecies!(1, "DATA_species/Si.in", atsp_vars, mtr_vars, apwlo_vars)
    readspecies!(2, "DATA_species/Pt.in", atsp_vars, mtr_vars, apwlo_vars)

    init_zero!( mtr_vars )

    checkmt!( latt_vars, atm_vars, atsp_vars.spsymb, mtr_vars )
    genrmesh!( atm_vars, atsp_vars, mtr_vars )
    init_packed_mtr!(mtr_vars)

    println("Before size(atsp_vars.vcln) = ", size(atsp_vars.vcln))
    init_nuclear_pot!( atsp_vars )
    println("After  size(atsp_vars.vcln) = ", size(atsp_vars.vcln))

    for is in 1:Nspecies
        println(atsp_vars.vcln[1:5,is])
    end
end

@time main()