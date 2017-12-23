include("../common/PWGrid_v01.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

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

function test_main( ns1::Int,ns2::Int,ns3::Int )
    #
    Ns = [ns1,ns2,ns3]
    const LatVecs = diagm( [6.0, 6.0, 6.0] )
    #
    pw = PWGrid( Ns, LatVecs )
    # Shortcuts
    const Npoints = pw.Npoints
    const Ω  = pw.Ω
    const r  = pw.r
    const G  = pw.G
    const G2 = pw.G2
    #
    # Generate array of distances
    #
    center = sum(LatVecs,2)/2
    dr = gen_dr( r, center )
    #
    # Setup potential
    #
    V_ionic = init_pot_harm_3d( pw, dr )
    #
    const Nstates = 4
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
    evals, evecs = eig(mu)
    Psi = Y*evecs
    for st = 1:Nstates
        @printf("State # %d, Energy = %18.10f\n", st, real(evals[st]))
    end
end

@time test_main( 30,30,30 )
