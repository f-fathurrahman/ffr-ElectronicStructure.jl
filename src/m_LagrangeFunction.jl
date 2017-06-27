module LagrangeFunction

export GvectorsT
include("LagrangeFunction/Gvectors.jl")

export LF1dGrid
export init_LF1d_p
export init_LF1d_c
export init_LF1d_sinc
include("LagrangeFunction/LF1d.jl")

export LF3dGrid
export init_LF3d_p
export init_LF3d_c
export init_LF3d_sinc
include("LagrangeFunction/LF3d.jl")

export ortho_gram_schmidt
include("LagrangeFunction/ortho_gram_schmidt.jl")

export orthonormalize
include("LagrangeFunction/orthonormalize.jl")

export print_Energies
include("LagrangeFunction/EnergiesT.jl")
include("LagrangeFunction/PotentialsT.jl")

include("LagrangeFunction/op_nabla2.jl")
include("LagrangeFunction/op_H.jl")

export calc_evals, calc_rho
include("LagrangeFunction/calc_Energies.jl")
include("LagrangeFunction/calc_evals.jl")
include("LagrangeFunction/calc_rho.jl")
include("LagrangeFunction/calc_grad.jl")

export KS_solve_scf, KS_solve_Emin_cg, KS_solve_Emin_pcg
include("LagrangeFunction/KS_solve_Emin_cg.jl")
include("LagrangeFunction/KS_solve_Emin_pcg.jl")
include("LagrangeFunction/KS_solve_scf.jl")

export prec_mkl_ilu0, get_Laplacian3d_kron
include("LagrangeFunction/sparse_LF3d.jl")
include("LagrangeFunction/prec_mkl_ilu0.jl")
include("LagrangeFunction/apply_prec_ilu0.jl")

include("LagrangeFunction/solve_poisson_FFT.jl")

include("LagrangeFunction/LDA_VWN.jl")
include("LagrangeFunction/diag_lobpcg.jl")

export build_nabla2
include("LagrangeFunction/build_nabla2.jl")

end
