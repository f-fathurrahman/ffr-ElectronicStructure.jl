push!(LOAD_PATH, "../LF_common/")
using m_LF3d

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

  lin2xyz = LF.lin2xyz

  nnzc = Nx + Ny + Nz - 2
  NNZ = nnzc*Npoints

  rowval = zeros( Int64, NNZ )
  nzval = zeros( Float64, NNZ )

  rowGbl_x_orig = zeros( Int64, Nx )
  rowGbl_y_orig = zeros( Int64, Ny )
  rowGbl_z_orig = zeros( Int64, Nz )

  rowGbl_x = zeros( Int64, Nx )
  rowGbl_y = zeros( Int64, Ny )
  rowGbl_z = zeros( Int64, Nz )

  for iz = 1:Nz
    rowGbl_z_orig[iz] = iz
  end

  rowGbl_y_orig[1] = 1
  for iy = 2:Ny
    rowGbl_y_orig[iy] = rowGbl_y_orig[iy-1] + Nz
  end

  rowGbl_x_orig[1] = 1
  for ix = 2:Nx
    rowGbl_x_orig[ix] = rowGbl_x_orig[ix-1] + Ny*Nz
  end

  # use @inbounds ?
  colGbl = 1
  idx_diag = 1

  for colGbl = 1:Npoints
    ix = lin2xyz[1,colGbl]
    iy = lin2xyz[2,colGbl]
    iz = lin2xyz[3,colGbl]
    println("\ncolGbl = ", colGbl, " : ix, iy, iz = ", [ix,iy,iz])
    ixx = ceil(Int,colGbl/(Ny*Nz))
    println(rowGbl_x_orig + colGbl - 1 - (ixx-1)*Nz*Ny )
    println(rowGbl_y_orig + 1 )
    println(rowGbl_z_orig + 1 )
  end


end


test_main( [3, 2, 4] )
