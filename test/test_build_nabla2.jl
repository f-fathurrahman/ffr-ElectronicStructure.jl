using ElectronicStructure

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

  srand(123)
  x = rand(Npoints)

  # Compare multiplication result
  Lx = nabla2 * x
  Lx_v2 = nabla2_v2 * x

  println( sum(Lx-Lx_v2) )

end

test_main( [10, 10, 10] )
test_main( [20, 20, 20] )
test_main( [30, 30, 30] )
test_main( [40, 40, 40] )
#test_main( [3,2,4] )
