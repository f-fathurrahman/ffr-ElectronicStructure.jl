include("../common/PWGrid_v03.jl")
include("../common/gen_dr_center.jl")
include("../common/calc_strfact_v2.jl")
include("../common/calc_ewald_v2.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("EnergiesT.jl")
include("PotentialsT.jl")
include("op_K.jl")
include("op_V_loc.jl")
include("op_H.jl")
include("calc_rho.jl")
include("calc_grad.jl")
include("calc_Energies.jl")
include("KS_solve_Emin_cg.jl")
include("Poisson_solve.jl")
include("LDA_VWN.jl")
include("Kprec.jl")

include("diag_lobpcg.jl")
include("KS_solve_SCF.jl")
include("andersonmix.jl")
include("KS_solve_SCF_andersonmix.jl")

include("init_V_gauss_G.jl")

function test_main( ecutwfc_Ry::Float64 )
    const LatVecs = 16.0*diagm( ones(3) )

    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs )

    const Ω  = pw.Ω
    const r  = pw.r
    const G  = pw.gvec.G
    const G2 = pw.gvec.G2
    const Ns = pw.Ns
    const Npoints = prod(Ns)
    const Ngwx = pw.gvecw.Ngwx

    @printf("ecutwfc (Ha) = %f\n", ecutwfc_Ry*0.5)
    @printf("Ns = (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])
    @printf("Ngwx = %d\n", Ngwx)
    const actual = Npoints/Ngwx
    const theor = 1/(4*pi*0.25^3/3)
    @printf("Compression: actual, theor: %f , %f\n", actual, theor)

    const Nstates = 2
    Zatm = 1.0
    Focc = [2.0, 0.0]
    Nelectrons = 2
    #
    Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )
    strf = calc_strfact( Xpos, 1, [1], pw.gvec.G )
    E_nn = calc_ewald( pw, strf, Xpos, 1, [1], [Zatm] )

    Nspecies = 1
    A = zeros(Nspecies)
    alpha = zeros(Nspecies)
    #
    A_in = 10.0
    alpha_in = 0.1
    A[1] = A_in/(2.0*pi*alpha_in^2)^1.5
    alpha[1] = 0.5/alpha_in^2
    #
    V_ionic = init_V_gauss_G(strf, Ω, Ns, G2, A, alpha)
    println("Integrated V_ionic = ", sum(V_ionic)*Ω/Npoints)

    Energies, Potentials, psi, evals = KS_solve_SCF( pw, V_ionic, Focc, Nstates, E_NN=E_nn )
    for st = 1:Nstates
        @printf("State # %d, Energy = %f\n", st, real(evals[st]))
    end
    print_Energies(Energies)

end

test_main(30.0)
