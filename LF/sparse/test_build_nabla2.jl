push!(LOAD_PATH, "../LF_common/")
using m_LF3d

include("../LF_common/sparse_LF3d.jl")
include("build_nabla2.jl")

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

  @time nabla2 = build_nabla2( LF )

  @time nabla2_v2 = get_Laplacian3d_kron( LF )

  #println(nabla2)
  #println(nabla2_v2)

  #println( full(nabla2) - full(nabla2_v2) )

  #println( nabla2 - nabla2_v2 )

  srand(123)
  x = rand(Npoints)

  # Compare multiplication result
  Lx = nabla2 * x
  Lx_v2 = nabla2_v2 * x

  println( sum(Lx-Lx_v2) )

end

test_main( [30, 30, 30] )
#test_main( [3,2,4] )
