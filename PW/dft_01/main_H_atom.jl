include("../common/PWGrid_v01.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

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

function test_main( ns1::Int,ns2::Int,ns3::Int )

    Ns = [ns1,ns2,ns3]
    const LatVecs = 16.0*diagm( ones(3) )

    pw = PWGrid( Ns, LatVecs )

    const Npoints = pw.Npoints
    const Ω = pw.Ω
    const r = pw.r
    const G = pw.G
    const G2 = pw.G2

    Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )

    Sf = structure_factor( Xpos, G )

    E_nn = calc_ewald( pw, Xpos, Sf )

    Vg = zeros(Complex128,Npoints)
    prefactor = -4*pi/Ω
    for ig=2:Npoints
        Vg[ig] = prefactor/G2[ig]
    end
    V_ionic = real( G_to_R(Ns, Vg .* Sf) ) * Npoints

    #
    const Nstates = 1
    Focc = [1.0]
    psi, Energies, Potentials = KS_solve_Emin_cg( pw, V_ionic, Focc, Nstates, NiterMax=1000 )

    #
    Y = ortho_gram_schmidt(psi)
    mu = Y' * op_H( pw, Potentials, Y )
    evals, evecs = eig(mu)
    psi = Y*evecs

    for st = 1:Nstates
        @printf("State # %d, Energy = %f\n", st, real(evals[st]))
    end

    @printf("E_nn    = %18.10f\n", E_nn)
    @printf("E total = %18.10f\n", E_nn + Energies.Total)

end

@code_native test_main(2,2,2)
@time test_main( 64, 64, 64 )
