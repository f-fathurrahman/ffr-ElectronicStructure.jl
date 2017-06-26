push!(LOAD_PATH, "../LF_common/")
using m_LF3d
include("../LF_common/sparse_LF3d.jl")

const ⊗ = kron

function test_main( NN::Array{Int64} )

  Nx = NN[1]
  Ny = NN[2]
  Nz = NN[3]

  AA = [0.0, 0.0, 0.0]
  BB = [16.0, 16.0, 16.0]

  # Initialize LF
  LF = init_LF3d_c( NN, AA, BB, verbose=true )

  @printf("Using loop")
  @time Ly_v1 = build_nabla2_y(LF)
  #println(Ly_v1)

  @printf("Using Kronecker product directly")
  @time Ly_v2 = speye(Nx) ⊗ LF.LFy.D2jl ⊗ speye(Nz)

  # Check whether there are any differences between the two matrices
  println(Ly_v1 - Ly_v2)

end

function build_nabla2_y(LF)

  D2jl = LF.LFy.D2jl

  Nx = LF.Nx
  Ny = LF.Ny
  Nz = LF.Nz

  Npoints = Nx*Ny*Nz

  # initialize with
  rowval = zeros(Int64,Ny*Npoints)
  colptr = zeros(Int64,Npoints+1)
  colptr[1] = 1
  nzval = zeros(Float64,Ny*Npoints)

  # XXX Move these arrays here ?
  rowGbl_init = zeros(Int64,Ny)
  rowGbl = zeros(Int64,Ny)

  idx = 0
  idxCol = 1

  for ix = 1:Nx

    print("\n\nix = ", ix)
    rowGbl_init[1] = 1 + (ix-1)*Ny*Nz
    for i = 2:Ny
      rowGbl_init[i] = rowGbl_init[i-1] + Nz
    end

    for colLoc = 1:Ny
      print("\n")
      colGbl_start = (colLoc-1)*Nz + 1 + (ix-1)*Ny*Nz
      colGbl_stop = colLoc*Nz + (ix-1)*Ny*Nz
      for colGbl = colGbl_start:colGbl_stop
        print("\ncolGbl: ", colGbl, " rowGbl: ")
        rowGbl[:] = rowGbl_init[:] + colGbl - 1 - (colLoc-1)*Nz - (ix-1)*Ny*Nz # shift
        for rowLoc = 1:Ny
          idx = idx + 1
          rowval[idx] = rowGbl[rowLoc]
          nzval[idx] = D2jl[colLoc,rowLoc]
          print(" ", rowGbl[rowLoc])
          print("(", rowLoc, ", ", colLoc, ")")
        end
        idxCol = idxCol + 1
        colptr[idxCol] = colptr[idxCol-1] + Ny
      end
    end

  end # ix = 1:Nx

  println("\n\nidx = ", idx)
  println("idx should be ", Ny*Npoints)

  println("\n\nidxCol = ", idxCol)
  println("idxCol should be ", Npoints+1)

  #print("\n")
  #exit()

  return SparseMatrixCSC( Npoints, Npoints, colptr, rowval, nzval )
end


test_main([6,4,5])
