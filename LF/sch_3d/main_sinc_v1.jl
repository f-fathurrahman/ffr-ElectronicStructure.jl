push!(LOAD_PATH, "../LF_common/")
using m_LF3d

include("init_pot_Hcoul.jl")
include("init_pot_harm_3d.jl")
include("../utils/ortho_gram_schmidt.jl")
include("../utils/orthonormalize.jl")
include("EnergiesT.jl")
include("op_nabla2.jl")
include("op_H.jl")
include("calc_Energies.jl")
include("calc_grad.jl")
include("schsolve_Emin_cg.jl")
include("schsolve_Emin_pcg.jl")
include("calc_evals.jl")
include("calc_rho.jl")
include("../LF_common/sparse_LF3d.jl")
include("../LF_common/prec_mkl_ilu0.jl")
include("../LF_common/apply_prec_ilu0.jl")
include("diag_lobpcg.jl")

function test_main( ; method = "Emin_cg" )
  # LF parameters
  NN = [20, 20, 20]
  hh = [0.3, 0.3, 0.3]

  Npoints = prod(NN)

  # Initialize LF
  LF = init_LF3d_sinc( NN, hh )
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h

  # Parameter for potential
  center = [0.0, 0.0, 0.0]

  ω = 2.0
  Vpot = init_pot_harm_3d( LF, ω, center )
  Ncols = 4

  #Vpot = init_pot_Hcoul( LF, center )
  #Ncols = 1

  if method == "Emin_cg_sparse"
    #
    ∇2 = get_Laplacian3d_kron(LF)
    prec = prec_mkl_ilu0( -0.5*∇2 + spdiagm(Vpot) )
    Energies, evecs = schsolve_Emin_pcg( LF, ∇2, prec, Vpot, Ncols, verbose=true )
    evals = calc_evals( LF, Vpot, evecs )
    #
  elseif method == "diag_lobpcg"
    #
    ∇2 = get_Laplacian3d_kron(LF)
    prec = prec_mkl_ilu0( -0.5*∇2 + spdiagm(Vpot) )
    srand(1234)
    v = rand( Npoints, Ncols )
    evals, evecs = diag_lobpcg( LF, ∇2, prec, Vpot, v, verbose=true)
    evecs = evecs/sqrt(ΔV) # normalize evecs
    Energies = calc_Energies( LF, ∇2, Vpot, evecs )
    #
  else
    #
    Energies, evecs = schsolve_Emin_cg( LF, Vpot, Ncols, verbose=true )
    evals = calc_evals( LF, Vpot, evecs )
  end

  @printf("Etot = %18.10f\n", Energies.Total)
  @printf("\nEigenvalues:\n")
  for i = 1:Ncols
    @printf("%8d %f\n", i, evals[i])
  end

  rho = get_rho( evecs )
  @printf( "\nintegRho = %lf\n", sum(rho)*ΔV )

end

@code_native test_main()
@time test_main(method="Emin_cg_sparse")
#@time test_main(method="diag_lobpcg")
#@time test_main(method="Emin_cg")
