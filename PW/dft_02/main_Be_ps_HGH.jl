include("../common/PWGrid_v02.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("EnergiesT.jl")
include("PotentialsT.jl")
include("op_K.jl")
include("op_V_loc.jl")
include("op_H.jl")
include("calc_rho.jl")
include("calc_grad_v2.jl")
include("calc_Energies.jl")
include("KS_solve_Emin_sd.jl")
include("KS_solve_Emin_cg.jl")
include("Poisson_solve.jl")
include("LDA_VWN.jl")
include("Kprec.jl")

include("../common/gen_dr_center.jl")
include("../common/calc_strfact_v2.jl")
include("../common/calc_ewald_v2.jl")

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

    # Helium atom
    Zion = 4.0
    const Nstates = 2
    Focc = [2.0, 2.0]

    Xpos = reshape( [1.0, 1.0, 1.0], (3,1) )

    Sf = calc_strfact( Xpos, 1, [1], pw.gvec.G )
    E_nn = calc_ewald( pw, Sf, Xpos, 1, [1], [Zion] )

    # Set up potential, using HGH pseudopotential for H
    # contains only local pseudopotential
    const rloc = 0.325000
    const c1 = -24.015041
    const c2 = 17.204014
    const c3 = -3.326390
    const c4 = 0.165419
    #
    Vg = zeros(Complex128,Npoints)
    pre1 = -4*pi*Zion/Ω
    pre2 = sqrt(8*pi^3)*rloc^3/Ω
    #
    for ig=2:Npoints
        Gr = sqrt(G2[ig])*rloc
        expGr2 = exp(-0.5*Gr^2)
        Vg[ig] = pre1/G2[ig]*expGr2 + pre2*expGr2 * (c1 + c2*(3-Gr^2) + 
                 c3*(15 - 10*Gr^2 + Gr^4) + c4*(105 - 105*Gr^2 + 21*Gr^4 - Gr^6) )
    end
    # limiting value
    Vg[1] = 2*pi*Zion*rloc^2 + (2*pi)^1.5 * rloc^3 * (c1 + 3.0*c2 + 15*c3 + 105*c4)
    println("Vg[1] = ", Vg[1])
    V_ionic = real( G_to_R(Ns, Vg .* Sf) ) * Npoints

    # Need to sum up over Nspecies for more than one species type
    # We simply need reshape because we only have one species type here.
    V_ionic = reshape( V_ionic, (Npoints) )

    println("sum(V_ionic) = ", sum(V_ionic))

    psi, Energies, Potentials = KS_solve_Emin_cg( pw, V_ionic, Focc, Nstates, NiterMax=1000, E_NN = E_nn )

    Y = ortho_gram_schmidt(psi)
    mu = Y' * op_H( pw, Potentials, Y )
    evals, evecs = eig(mu)
    psi = Y*evecs

    #@printf("Solution by self-consistent field method\n")
    #Energies, Potentials, psi, evals = KS_solve_scf( pw, V_ionic, Focc, Nstates, β=0.3 )

    for ist = 1:Nstates
        @printf("State # %d, Energy = %f\n", ist, real(evals[ist]))
    end

    print_Energies(Energies)

end

@time test_main( [54,54,54] )
#@time test_main( [64,64,64] )
#@time test_main( [100,100,100] )
