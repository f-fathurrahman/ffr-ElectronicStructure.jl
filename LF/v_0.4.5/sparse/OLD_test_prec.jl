push!(LOAD_PATH, "../LF_common/")
using m_LF3d
include("../LF_common/sparse_LF3d.jl")
include("../LF_common/prec_mkl_ilu0.jl")
include("../LF_common/apply_prec_ilu0.jl")
include("build_nabla2_x.jl")
include("build_nabla2_y.jl")
include("build_nabla2_z.jl")

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

  nabla2_x = build_nabla2_x(LF)
  nabla2_y = build_nabla2_y(LF)
  nabla2_z = build_nabla2_z(LF)

  # Using Kronecker product
  ∇2 = get_Laplacian3d_kron(LF)

  #
  srand(123)
  v = rand(Npoints)

  # Test multiplication operation
  #@time Lv_1 = nabla2_x*v + nabla2_y*v + nabla2_z*v
  #@time Lv_2 = ∇2*v

  #@printf("Difference: %18.10e\n", sum(abs(Lv_1 - Lv_2))/Npoints)


  # Test the preconditioner

  prec_x = prec_mkl_ilu0( nabla2_x )
  prec_y = prec_mkl_ilu0( nabla2_y )
  prec_z = prec_mkl_ilu0( nabla2_z )

  pvx = apply_prec_ilu0( prec_x, v )
  pvy = apply_prec_ilu0( prec_y, v )
  pvz = apply_prec_ilu0( prec_z, v )
  Pv_1 = pvx + pvy + pvz
  #
  prec_∇2 = prec_mkl_ilu0( ∇2 )
  Pv_2 = apply_prec_ilu0(prec_∇2, v)

  println(sum(v))

  println(sum(pvx))
  println(sum(pvy))
  println(sum(pvz))

  println(sum(Pv_1))
  println(sum(Pv_2))


end

test_main( [3, 3, 3])
