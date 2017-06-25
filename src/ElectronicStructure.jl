__precompile__()

module ElectronicStructure

export PWGrid

export GvectorsT
include("m_Gvectors.jl")

export LF1dGrid
export init_LF1d_p
export init_LF1d_c
export init_LF1d_sinc
include("m_LF1d.jl")

export LF3dGrid
export init_LF3d_p
export init_LF3d_c
export init_LF3d_sinc
include("m_LF3d.jl")

export ortho_gram_schmidt
include("ortho_gram_schmidt.jl")

export orthonormalize
include("orthonormalize.jl")

export print_Energies
include("EnergiesT.jl")
include("PotentialsT.jl")

include("op_nabla2.jl")
include("op_H.jl")

export calc_evals, calc_rho
include("calc_Energies.jl")
include("calc_evals.jl")
include("calc_rho.jl")
include("calc_grad.jl")

export KS_solve_scf, KS_solve_Emin_cg, KS_solve_Emin_pcg
include("KS_solve_Emin_cg.jl")
include("KS_solve_Emin_pcg.jl")
include("KS_solve_scf.jl")

export prec_mkl_ilu0, get_Laplacian3d_kron
include("sparse_LF3d.jl")
include("prec_mkl_ilu0.jl")
include("apply_prec_ilu0.jl")

include("solve_poisson_FFT.jl")

include("LDA_VWN.jl")
include("diag_lobpcg.jl")

export build_nabla2
include("build_nabla2.jl")

end
