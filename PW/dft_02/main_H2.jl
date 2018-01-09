include("../common/PWGrid_v02.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("EnergiesT.jl")
include("PotentialsT.jl")
include("gen_dr.jl")
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

include("../common/gen_dr_center.jl")
include("../common/calc_strfact_v2.jl")
include("../common/calc_ewald_v2.jl")

include("diag_lobpcg.jl")
include("KS_solve_SCF.jl")


function test_main( Ns, xx )

    const LatVecs = 16.0*diagm( ones(3) )

    pw = PWGrid( Ns, LatVecs )

    const Ω  = pw.Ω
    const r  = pw.r
    const G  = pw.gvec.G
    const G2 = pw.gvec.G2
    const Npoints = prod(Ns)
    const Ngwx = pw.gvecw.Ngwx

    Xpos = zeros( 3, 2 )
    Xpos[:,1] = [0.0, 0.0, 0.0]
    Xpos[:,2] = [xx, 0.0, 0.0]

    println(Xpos')

    Nspecies = 1
    atmsymb = ["H"] # unique list of atomic symbols
    atm2species = [1, 1]    # mapping from atom to species
    Zv = [1.0]    # only valence ?

    Sf = calc_strfact( Xpos, Nspecies, atm2species, pw.gvec.G )

    E_nn = calc_ewald( pw, Sf, Xpos, Nspecies, atm2species, Zv )
    #@printf("E_nn = %18.10f\n", E_nn)

    Vg = zeros(Complex128,Npoints)
    prefactor = -4*pi/Ω
    for ig=2:Npoints
        Vg[ig] = prefactor/G2[ig]
    end
    V_ionic = real( G_to_R(Ns, Vg .* Sf) ) * Npoints

    # Need to sum up over Nspecies for more than one species type
    # We simply need reshape because we only have one species type here.
    V_ionic = reshape( V_ionic, (Npoints) )

    #println("sum(Sf) = ", sum(Sf))
    #println("size(V_ionic) = ", size(V_ionic))

    const Nstates = 1
    Focc = [2.0]

    psi, Energies, Potentials = KS_solve_Emin_cg( pw, V_ionic, Focc, Nstates, NiterMax=1000 )

    Y = ortho_gram_schmidt(psi)
    mu = Y' * op_H( pw, Potentials, Y )
    evals, evecs = eig(mu)
    psi = Y*evecs

    #Energies, Potentials, psi, evals = KS_solve_SCF( pw, V_ionic, Focc, Nstates )

    for st = 1:Nstates
        @printf("=== State # %d, Energy = %f ===\n", st, real(evals[st]))
    end

    #@printf("E_nn = %18.10f\n", E_nn)
    @printf("E total = %18.10f\n", E_nn + Energies.Total)

end


@time test_main( [64,64,64], 1.5 )

#for xx in [0.5, 1.0, 1.25, 1.50, 1.75, 2.0, 4.0, 6.0]
#    test_main( [64,64,64], xx )
#end
