push!(LOAD_PATH, "../LF_common/")
using m_LF3d
include("../LF_common/sparse_LF3d.jl")
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

  @time nabla2_x = build_nabla2_x(LF)
  @time nabla2_y = build_nabla2_y(LF)
  @time nabla2_z = build_nabla2_z(LF)

  @time combine_nabla2( NN, nabla2_x, nabla2_y, nabla2_z )

end

function combine_nabla2( NN, nabla2_x, nabla2_y, nabla2_z )

  Nx = NN[1]
  Ny = NN[2]
  Nz = NN[3]
  Npoints = Nx*Ny*Nz

  # number of non-zero per-column
  nnz = Nx + Ny + Nz - 2

  # combined results
  rowval = zeros(Int64,nnz*Npoints)
  colptr = zeros(Int64,Npoints+1)
  colptr[1] = 1
  nzval = zeros(Float64,nnz*Npoints)

  # shortcuts
  x_rowval = nabla2_x.rowval
  y_rowval = nabla2_y.rowval
  z_rowval = nabla2_z.rowval
  #
  x_nzval = nabla2_x.nzval
  y_nzval = nabla2_y.nzval
  z_nzval = nabla2_z.nzval

  c_rowval = zeros(Int64, Nx+Ny+Nz)  # temporary combined indices of row
  c_nzval  = zeros(Float64, Nx+Ny+Nz)

  idxSorted = zeros(Int64, Nx+Ny+Nz)

  idxGbl = 0

  println(size(x_rowval))

  compxyz = Array{ASCIIString}(Nx+Ny+Nz)
  compxyz[1:Nx] = "x"
  compxyz[Nx+1:Nx+Ny] = "y"
  compxyz[Nx+Ny+1:Nx+Ny+Nz] = "z"

  idxDiag = 0
  for colGbl = 1:Npoints
    #
    idxDiag = idxDiag + 1
    #println("\nidxDiag = ", idxDiag)
    #
    idx = 0
    for i in range(1,Nx) + Nx*(colGbl-1)
      idx = idx + 1
      c_rowval[idx] = x_rowval[i]
      c_nzval[idx] = x_nzval[i]
    end
    for j in range(1,Ny) + Ny*(colGbl-1)
      idx = idx + 1
      c_rowval[idx] = y_rowval[j]
      c_nzval[idx] = y_nzval[j]
    end
    for k in range(1,Nz) + Nz*(colGbl-1)
      idx = idx + 1
      c_rowval[idx] = z_rowval[k]
      c_nzval[idx] = z_nzval[k]
    end
    # Sort the indices
    idxSorted = sortperm(c_rowval)
    println(compxyz[idxSorted])
    #println(c_rowval[idxSorted])
    #
  end


end


test_main( [3, 2, 4] )
#test_main( [40, 40, 40] )
#test_main( [80, 80, 80] )
