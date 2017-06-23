__precompile__()

module ElectronicStructure

export PWGrid

export GvectorsT
include("../LF/LF_common/m_Gvectors.jl")

export LF1dGrid
export init_LF1d_p
export init_LF1d_c
export init_LF1d_sinc
include("../LF/LF_common/m_LF1d.jl")

export LF3dGrid
export init_LF3d_p
export init_LF3d_c
export init_LF3d_sinc
include("../LF/LF_common/m_LF3d.jl")

export ortho_gram_schmidt
include("../LF/utils/ortho_gram_schmidt.jl")

export orthonormalize
include("../LF/utils/orthonormalize.jl")

export init_pot_harm_3d
include("../LF/dft_3d_p_v1/init_pot_harm_3d.jl")

export print_Energies
include("../LF/dft_3d_p_v1/EnergiesT.jl")
include("../LF/dft_3d_p_v1/PotentialsT.jl")

include("../LF/dft_3d_p_v1/op_nabla2.jl")
include("../LF/dft_3d_p_v1/op_H.jl")

export calc_evals, calc_rho
include("../LF/dft_3d_p_v1/calc_Energies.jl")
include("../LF/dft_3d_p_v1/calc_evals.jl")
include("../LF/dft_3d_p_v1/calc_rho.jl")
include("../LF/dft_3d_p_v1/calc_grad.jl")

export KS_solve_scf, KS_solve_Emin_cg, KS_solve_Emin_pcg
include("../LF/dft_3d_p_v1/KS_solve_Emin_cg.jl")
include("../LF/dft_3d_p_v1/KS_solve_Emin_pcg.jl")
include("../LF/dft_3d_p_v1/KS_solve_scf.jl")

export prec_mkl_ilu0, get_Laplacian3d_kron
include("../LF/LF_common/sparse_LF3d.jl")
include("./prec_mkl_ilu0.jl")
include("./apply_prec_ilu0.jl")

include("../LF/LF_common/solve_poisson_FFT.jl")

include("../LF/dft_3d_p_v1/LDA_VWN.jl")
include("../LF/dft_3d_p_v1/diag_lobpcg.jl")

end
