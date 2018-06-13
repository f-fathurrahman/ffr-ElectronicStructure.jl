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

include("structure_factor.jl")
include("gen_rho.jl")
include("gen_dr.jl")
include("calc_ewald.jl")

function test_main( ns1::Int64,ns2::Int64,ns3::Int64 )

    Ns = [ns1,ns2,ns3]
    LatVecs = gen_lattice_sc(16.0)

    pw = PWGrid( Ns, LatVecs )

    Npoints = pw.Npoints
    Ω = pw.Ω
    r = pw.r
    G = pw.G
    G2 = pw.G2

    Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )

    Sf = structure_factor( Xpos, G )

    E_nn = calc_ewald( pw, Xpos, Sf )

    Vg = zeros(ComplexF64,Npoints)
    prefactor = -4*pi/Ω
    for ig=2:Npoints
        Vg[ig] = prefactor/G2[ig]
    end
    V_ionic = real( G_to_R(Ns, Vg .* Sf) ) * Npoints

    #
    Nstates = 1
    Focc = [1.0]
    psi, Energies, Potentials = KS_solve_Emin_cg( pw, V_ionic, Focc, Nstates, NiterMax=1000 )

    #
    Y = ortho_gram_schmidt(psi)
    mu = Y' * op_H( pw, Potentials, Y )
    evals, evecs = eigen(mu)
    psi = Y*evecs

    for st = 1:Nstates
        @printf("State # %d, Energy = %f\n", st, real(evals[st]))
    end

    @printf("E_nn    = %18.10f\n", E_nn)
    @printf("E total = %18.10f\n", E_nn + Energies.Total)

end

@time test_main( 64, 64, 64 )
