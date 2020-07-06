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

    atsp_vars.nstspmax = maximum(atsp_vars.nstsp)

    init_zero!( mtr_vars )

    println("before nrsp = ", atsp_vars.nrsp)
    println("atsp_vars.rsp = ", size(atsp_vars.rsp))

    checkmt!( latt_vars, atm_vars, atsp_vars.spsymb, mtr_vars )
    genrmesh!( atm_vars, atsp_vars, mtr_vars )
    init_packed_mtr!(mtr_vars)

    println("after nrsp = ", atsp_vars.nrsp)
    println("atsp_vars.rsp = ", size(atsp_vars.rsp))

    nrspmax = atsp_vars.nrspmax
    nspecies = atm_vars.nspecies
    nstspmax = atsp_vars.nstspmax

    # allocate global species charge density and potential arrays
    rhosp = zeros(Float64,nrspmax,nspecies)
    vrsp = zeros(Float64,nrspmax,nspecies)
    
    xctsp = atsp_vars.xctsp
    xcgrad = false

    spzn = atsp_vars.spzn
    nstsp = atsp_vars.nstsp
    nsp = atsp_vars.nsp
    lsp = atsp_vars.lsp
    ksp = atsp_vars.ksp
    occsp = atsp_vars.occsp
    nrsp = atsp_vars.nrsp
    rsp = atsp_vars.rsp
    evalsp = atsp_vars.evalsp
    ptnucl = atsp_vars.ptnucl

    rwf = zeros(Float64,nrspmax,2,nstspmax)
    
    # speed of light in atomic units (=1/alpha) (CODATA 2018)
    sol = 137.035999084
    # scaled speed of light
    solsc = sol

    for is in 1:nspecies
        @views solve_atom!(
            solsc, ptnucl, spzn[is], nstsp[is], nsp[:,is], lsp[:,is], ksp[:,is],
            occsp[:,is], xctsp, xcgrad, nrsp[is], rsp[:,is], evalsp[:,is], rhosp[:,is],
            vrsp[:,is], rwf
        )
        for ist in 1:nstsp[is]
            @printf("%3d %18.10f\n", ist, evalsp[ist,is])
        end
    end

end

@time main()