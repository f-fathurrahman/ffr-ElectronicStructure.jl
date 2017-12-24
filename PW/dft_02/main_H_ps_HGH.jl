include("../common/PWGrid_v02.jl")
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
include("KS_solve_Emin_sd.jl")
include("KS_solve_Emin_cg.jl")
include("Poisson_solve.jl")
include("LDA_VWN.jl")
include("Kprec.jl")

include("calc_strfact_v1.jl")
include("gen_dr.jl")
include("calc_ewald_v1.jl")

include("diag_lobpcg.jl")
include("KS_solve_scf.jl")


function test_main( Ns )

    const LatVecs = 16.0*diagm( ones(3) )

    pw = PWGrid( Ns, LatVecs )

    const Ω  = pw.Ω
    const r  = pw.r
    const G  = pw.gvec.G
    const G2 = pw.gvec.G2
    const Npoints = prod(Ns)
    const Ngwx = pw.gvecw.Ngwx

    @printf("Ns   = (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])
    @printf("Ngwx = %d\n", Ngwx)

    const actual = Npoints/Ngwx
    const theor = 1/(4*pi*0.25^3/3)
    @printf("Compression: actual, theor: %f , %f\n", actual, theor)

    Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )

    Sf = calc_strfact( Xpos, G )

    E_nn = calc_ewald( pw, Xpos, Sf )

    # Set up potential, using HGH pseudopotential for H
    # contains only local pseudopotential
    const Zion = 1.0
    const rloc = 0.2
    const c1 = -4.180237
    const c2 = 0.725075
    Vg = zeros(Complex128,Npoints)
    pre1 = -4*pi*Zion/Ω
    pre2 = sqrt(8*pi^3)*rloc^3/Ω
    #
    for ig=2:Npoints
        Gr = sqrt(G2[ig])*rloc
        expGr2 = exp(-0.5*Gr^2)
        Vg[ig] = pre1/G2[ig]*expGr2 + pre2*expGr2 * (c1 + c2*(3-Gr^2) )
    end
    Vg[1] = 2*pi*rloc^2 + (2*pi)^1.5 * rloc^3 * (c1 + 3.0*c2) # limiting value
    V_ionic = real( G_to_R(Ns, Vg .* Sf) ) * Npoints

    #
    const Nstates = 1
    Focc = [1.0]

    psi, Energies, Potentials = KS_solve_Emin_cg( pw, V_ionic, Focc, Nstates, NiterMax=1000, E_NN = E_nn )

    Y = ortho_gram_schmidt(psi)
    mu = Y' * op_H( pw, Potentials, Y )
    evals, evecs = eig(mu)
    psi = Y*evecs

    #Energies, Potentials, psi, evals = KS_solve_scf( pw, V_ionic, Focc, Nstates )

    for st = 1:Nstates
        @printf("State # %d, Energy = %f\n", st, real(evals[st]))
    end

    print_Energies(Energies)

end

@code_native test_main( [64,64,64] )
@time test_main( [64,64,64] )
