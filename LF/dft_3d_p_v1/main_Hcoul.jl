include("../LF_common/m_LF1d.jl")
include("../LF_common/m_LF3d.jl")
include("../LF_common/m_Gvectors.jl")

include("../utils/ortho_gram_schmidt.jl")
include("../utils/orthonormalize.jl")

include("init_pot_Hcoul.jl")

include("EnergiesT.jl")
include("PotentialsT.jl")

include("op_nabla2.jl")
include("op_H.jl")

include("calc_Energies.jl")
include("calc_evals.jl")
include("calc_rho.jl")
include("calc_grad.jl")

include("KS_solve_Emin_cg.jl")
include("KS_solve_Emin_pcg.jl")
include("KS_solve_scf.jl")

include("../LF_common/sparse_LF3d.jl")
include("../LF_common/prec_mkl_ilu0.jl")
include("../LF_common/apply_prec_ilu0.jl")

include("../LF_common/solve_poisson_FFT.jl")

include("LDA_VWN.jl")
include("diag_lobpcg.jl")


function test_main( ; method = "Emin_cg" )
    # LF parameters
    NN = [37, 37, 37]
    AA = [-8.0, -8.0, -8.0]
    BB = [ 8.0,  8.0,  8.0]

    Npoints = prod(NN)

    # Initialize LF
    LF = init_LF3d_p( NN, AA, BB )
    ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h

    # Initialize G-vectors
    L = BB - AA
    Gv = GvectorsT( NN, diagm(L) )

    Vg = zeros(Complex128,Npoints)
    for ig=2:Npoints
        Vg[ig] = -4*pi/Gv.G2[ig]/Gv.Ω
    end
    V_ionic = real( G_to_R(NN,Vg) ) * Npoints

    println("sum(V_ionic): ", sum(V_ionic))

    Ncols = 1
    Focc = 1.0*ones(Ncols)

    if method == "Emin_cg_sparse"
        #
        ∇2 = get_Laplacian3d_kron(LF)
        precH = prec_mkl_ilu0( -0.5*∇2 + spdiagm(V_ionic) )
        #precH = speye(Npoints)
        Energies, evecs, Potentials = KS_solve_Emin_pcg( LF, Gv, ∇2, precH,
                                        V_ionic, Focc, Ncols, verbose=true )
        evals = calc_evals( LF, ∇2, Potentials, evecs )
        #
    elseif method == "SCF"
        ∇2 = get_Laplacian3d_kron(LF)
        precH = prec_mkl_ilu0( -0.5*∇2 + spdiagm(V_ionic) )
        #precH = speye(Npoints)
        Output = KS_solve_scf( LF, Gv, ∇2, precH, V_ionic, Focc, Ncols, verbose=true )
        exit()
    else
        #
        Energies, evecs, Potentials = KS_solve_Emin_cg( LF, Gv, V_ionic, Focc, Ncols, verbose=true )
        evals = calc_evals( LF, Potentials, evecs )
    end

    print_Energies(Energies)
    @printf("\nEigenvalues:\n")
    for i = 1:Ncols
        @printf("%8d %f\n", i, evals[i])
    end

    rho = calc_rho( Focc, evecs )
    @printf( "\nintegRho = %lf\n", sum(rho)*ΔV )

end

#@code_native test_main()
@time test_main(method="Emin_cg_sparse")
@time test_main(method="SCF")
@time test_main(method="Emin_cg")
