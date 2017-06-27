using ElectronicStructure
using ElectronicStructure.LagrangeFunction

function test_main( NN::Array{Int64} )

  println("\nNN = ", NN)

  Nx = NN[1]
  Ny = NN[2]
  Nz = NN[3]
  Npoints = Nx*Ny*Nz

  AA = [0.0, 0.0, 0.0]
  BB = [16.0, 16.0, 16.0]

  # Initialize LF
  LF = init_LF3d_c( NN, AA, BB )

  @printf("build_nabla2: ")
  @time nabla2 = build_nabla2( LF )

  @printf("get_Laplacian3d_kron: ")
  @time nabla2_v2 = get_Laplacian3d_kron( LF )

  srand(123)
  x = rand(Npoints)

  # Compare multiplication result
  Lx = nabla2 * x
  Lx_v2 = nabla2_v2 * x

  @printf("Diff = ")
  println( sum(Lx-Lx_v2) )

end

test_main( [10, 10, 10] )
test_main( [20, 20, 20] )
test_main( [30, 30, 30] )
test_main( [40, 40, 40] )
#test_main( [3,2,4] )
