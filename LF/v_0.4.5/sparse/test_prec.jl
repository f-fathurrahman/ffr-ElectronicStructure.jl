push!(LOAD_PATH, "../LF_common/")
using m_LF3d
include("../LF_common/sparse_LF3d.jl")
include("../LF_common/prec_mkl_ilu0.jl")
include("../LF_common/apply_prec_ilu0.jl")
include("build_nabla2.jl")

const ⊗ = kron

function test_main( NN::Array{Int64} )

  println("NN = ", NN)

  Nx = NN[1]
  Ny = NN[2]
  Nz = NN[3]
  Npoints = Nx*Ny*Nz

  AA = [0.0, 0.0, 0.0]
  BB = [16.0, 16.0, 16.0]

  # Initialize LF
  LF = init_LF3d_c( NN, AA, BB, verbose=true )

  # Directly build Laplacian and ILU preconditioner from MKL
  @time nabla2 = build_nabla2(LF)
  @time prec_nabla2 = prec_mkl_ilu0( nabla2 )

  # Using Kronecker product and ILU0 preconditioner from MKL
  @time ∇2 = get_Laplacian3d_kron(LF)
  @time prec_∇2 = prec_mkl_ilu0( ∇2 )

  # Test the preconditioner
  srand(123)
  v = rand(Npoints)

  Pv_1 = apply_prec_ilu0( prec_nabla2, v )
  Pv_2 = apply_prec_ilu0( prec_∇2, v )

  println( sum(Pv_1 - Pv_2) )

end

test_main( [50, 50, 50])
