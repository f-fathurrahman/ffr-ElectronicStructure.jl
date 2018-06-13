using Printf
using LinearAlgebra
using Random

include("../common/PWGrid_v01.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")
include("../common/gen_lattice_pwscf.jl")
include("EnergiesT.jl")
include("PotentialsT.jl")
include("gen_dr.jl")
include("init_pot_harm_3d.jl")
include("op_K.jl")
include("op_V_ps_loc.jl")
include("op_H.jl")
include("calc_rho.jl")
include("calc_grad.jl")
include("calc_Energies.jl")
include("KS_solve_Emin_sd.jl")
include("KS_solve_Emin_cg.jl")
include("Poisson_solve.jl")
include("LDA_VWN.jl")
include("Kprec.jl")

function test_main( ns1::Int64,ns2::Int64,ns3::Int64 )
    #
    Ns = [ns1,ns2,ns3]
    LatVecs = gen_lattice_sc(6.0)
    #
    pw = PWGrid( Ns, LatVecs )
    # Shortcuts
    Npoints = pw.Npoints
    Ω  = pw.Ω
    r  = pw.r
    G  = pw.G
    G2 = pw.G2
    #
    # Generate array of distances
    #
    center = sum(LatVecs,dims=2)/2
    dr = gen_dr( r, center )
    #
    # Setup potential
    #
    V_ionic = init_pot_harm_3d( pw, dr )
    #
    Nstates = 4
    Focc = 2.0*ones(Nstates)
    #
    #psi, Energies, Potentials = kssolve_Emin_sd( pw, V_ionic, Focc, Nstates, NiterMax=10 )
    #psi, Energies, Potentials = kssolve_Emin_cg( pw, V_ionic, Focc, Nstates,
    #                            NiterMax=1000, Potentials0=Potentials, psi0=psi )
    psi, Energies, Potentials = KS_solve_Emin_cg( pw, V_ionic, Focc, Nstates, NiterMax=1000 )
    #psi, Energies = kssolve_Emin_cg( pw, Vpot, psi, NiterMax=1000 )
    #
    Y = ortho_gram_schmidt(psi)
    mu = Y' * op_H( pw, Potentials, Y )
    evals, evecs = eigen(mu)
    Psi = Y*evecs
    for st = 1:Nstates
        @printf("State # %d, Energy = %18.10f\n", st, real(evals[st]))
    end
end

@time test_main( 30,30,30 )
