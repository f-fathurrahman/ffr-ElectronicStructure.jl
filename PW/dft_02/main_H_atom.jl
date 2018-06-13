using Printf
using LinearAlgebra
using Random

include("../common/PWGrid_v02.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")
include("../common/gen_lattice_pwscf.jl")
include("EnergiesT.jl")
include("PotentialsT.jl")
include("op_K.jl")
include("op_V_loc.jl")
include("op_H.jl")
include("calc_rho.jl")
include("calc_grad.jl")
include("calc_Energies.jl")
include("KS_solve_Emin_sd.jl")
include("KS_solve_MGC_cg.jl")
include("KS_solve_Emin_cg.jl")
include("Poisson_solve.jl")
include("LDA_VWN.jl")
include("Kprec.jl")

include("gen_dr.jl")
include("calc_strfact_v1.jl")
include("calc_ewald_v1.jl")

include("diag_lobpcg.jl")
include("KS_solve_SCF.jl")

include("andersonmix.jl")
include("KS_solve_SCF_andersonmix.jl")

include("pulaymix.jl")
include("KS_solve_SCF_pulaymix.jl")

include("get_ub_lb_lanczos.jl")
include("KS_solve_ChebySCF.jl")
include("chebyfilt.jl")
include("norm_matrix_induced.jl")

function test_main( Ns; method="SCF" )

    LatVecs = gen_lattice_sc(16.0)

    pw = PWGrid( Ns, LatVecs )

    Ω  = pw.Ω
    r  = pw.r
    G  = pw.gvec.G
    G2 = pw.gvec.G2
    Npoints = prod(Ns)
    Ngwx = pw.gvecw.Ngwx

    @printf("Ns   = (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])
    @printf("Ngwx = %d\n", Ngwx)

    actual = Npoints/Ngwx
    theor = 1/(4*pi*0.25^3/3)
    @printf("Compression: actual, theor: %f , %f\n", actual, theor)

    Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )

    Sf = calc_strfact( Xpos, G )

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

    if method == "CG"
        psi, Energies, Potentials = KS_solve_Emin_cg( 
            pw, V_ionic, Focc, Nstates, NiterMax=1000, E_NN=E_nn )
        #
        Y = ortho_gram_schmidt(psi)
        mu = Y' * op_H( pw, Potentials, Y )
        evals, evecs = eigen(mu)
        psi = Y*evecs
    elseif method == "SCF"
        Energies, Potentials, psi, evals = KS_solve_SCF( pw, V_ionic, Focc, Nstates, E_NN=E_nn )
    elseif method == "ChebySCF"
        Energies, Potentials, psi, evals = KS_solve_ChebySCF( pw, V_ionic, Focc, Nstates, β=0.8, E_NN=E_nn )
    elseif method == "SCF_andersonmix"
        Energies, Potentials, psi, evals = KS_solve_SCF_andersonmix( pw, V_ionic, Focc, Nstates, β=0.5, E_NN=E_nn )
    elseif method == "SCF_pulaymix"
        Energies, Potentials, psi, evals = KS_solve_SCF_pulaymix( pw, V_ionic, Focc, Nstates, β=0.5, E_NN=E_nn)
    else
        println("ERROR: unknown method: ", method)
        exit()
    end

    for st = 1:Nstates
        @printf("State # %d, Energy = %18.10f\n", st, real(evals[st]))
    end

    print_Energies(Energies)

end

@time test_main( [66,66,66], method="CG" )
#@time test_main( [64,64,64], method="SCF" )
#@time test_main( [64,64,64], method="ChebySCF" )
#@time test_main( [64, 64, 64], method="SCF_andersonmix" )
#@time test_main( [64, 64, 64], method="SCF_pulaymix")
