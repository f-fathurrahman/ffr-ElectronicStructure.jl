push!(LOAD_PATH, "../LF_common/")
using m_LF3d

include("../LF_common/sparse_LF3d.jl")

using PyPlot
const plt = PyPlot

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
  D2jl_x = LF.LFx.D2jl
  D2jl_y = LF.LFy.D2jl
  D2jl_z = LF.LFz.D2jl

  nnzc = Nx + Ny + Nz - 2
  NNZ = nnzc*Npoints

  rowval = zeros( Int64, NNZ )
  nzval = zeros( Float64, NNZ )
  colptr = zeros( Int64, Npoints + 1 )

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

  ip = 0
  # rowLoc = 1:Nx, 1:Ny, 1:Nz for each gblCol
  for colGbl = 1:Npoints

    colLoc_x = ceil(Int,colGbl/(Ny*Nz))
    #println( rowGbl_x_orig + colGbl - 1 - (colLoc_x-1)*Nz*Ny )

    izz = ceil(Int,colGbl/Nz)
    #println( rowGbl_y_orig + colGbl - 1 - (izz-1)*Nz + (colLoc_x-1)*Ny*Nz )

    #println( rowGbl_z_orig + (izz-1)*Nz )

    yy = colGbl - (colLoc_x-1)*Ny*Nz
    colLoc_y = ceil( Int, yy/Nz )

    colLoc_z = colGbl - (izz-1)*Nz

    #println("coLLoc x, y, z:", [colLoc_x, colLoc_y, colLoc_z])
    #println("")

    # diagonal element
    ip = ip + 1
    rowval[ip] = colGbl
    nzval[ip] = D2jl_x[colLoc_x,colLoc_x] + D2jl_y[colLoc_y,colLoc_y] +  D2jl_z[colLoc_z,colLoc_z]
    #println("diagonal: ", rowval[ip])

    # non
    for ir = 1:Nx
      if ir != colLoc_x
        ip = ip + 1
        rowval[ip] = rowGbl_x_orig[ir] + colGbl - 1 - (colLoc_x-1)*Nz*Ny
        nzval[ip] = D2jl_x[ir,colLoc_x]
        #println("x: ", rowval[ip])
      end
    end

    for ir = 1:Ny
      if ir != colLoc_y
        ip = ip + 1
        rowval[ip] = rowGbl_y_orig[ir] + colGbl - 1 - (izz-1)*Nz + (colLoc_x-1)*Ny*Nz
        nzval[ip] = D2jl_y[ir,colLoc_y]
        #println("y: ", rowval[ip])
      end
    end

    for ir = 1:Nz
      if ir != colLoc_z
        ip = ip + 1
        rowval[ip] = rowGbl_z_orig[ir] + (izz-1)*Nz
        nzval[ip] = D2jl_z[ir,colLoc_z]
        #println("z: ", rowval[ip])
      end
    end

  end

  # colptr
  colptr[1] = 1
  for i = 2:Npoints+1
    colptr[i] = colptr[i-1] + nnzc
  end

  #println("ip = ", ip)
  #println("NNZ = ", NNZ)

  #println(colptr)
  #println(rowval)

  nabla2 = SparseMatrixCSC( Npoints, Npoints, colptr, rowval, nzval )

  nabla2_v2 = get_Laplacian3d_kron( LF )

  #println(nabla2)
  #println(nabla2_v2)

  #println( full(nabla2) - full(nabla2_v2) )

  #println( nabla2 - nabla2_v2 )

  srand(123)
  x = rand(Npoints)

  Lx = nabla2 * x
  Lx_v2 = nabla2_v2 * x

  println( sum(Lx-Lx_v2) )

  #plt.clf()
  #plt.spy(full(nabla2))
  #plt.title("nabla")
  #plt.savefig("nabla2.pdf")

  #plt.clf()
  #plt.spy(full(nabla2_v2))
  #plt.title("nabla_v2")
  #plt.savefig("nabla2_v2.pdf")

end


test_main( [30, 30, 30] )
