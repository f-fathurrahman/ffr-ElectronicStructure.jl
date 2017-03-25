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
  #println(Lx_v1)

  #@printf("Using Kronecker product directly")
  #@time Ly_v2 = speye(Nx) ⊗ LF.LFy.D2jl ⊗ speye(Nz)

  # Check whether there are any differences between the two matrices
  #println(Ly_v1 - Ly_v2)
  #println(Lx_v2)

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

  idx = 1

  # XXX Move these arrays here ?
  rowGbl_init = zeros(Int64,Ny)
  rowGbl = zeros(Int64,Ny)

  #for ix = 1:Nx
    ix = 1
    print("\nix = ", ix)
    rowGbl_init[1] = 1 + (ix-1)*Ny*Nz
    for i = 2:Ny
      rowGbl_init[i] = rowGbl_init[i-1] + Nz
    end
  #ix = 1
  for colLoc = 1:Ny
    print("\n")
    colGbl_start = (colLoc-1)*Nz + 1
    colGbl_stop = colLoc*Nz
    for colGbl = colGbl_start:colGbl_stop
      print("\ncolGbl: ", colGbl, " rowGbl: ")
      rowGbl[:] = rowGbl_init[:] + colGbl - 1 - (colLoc-1)*Nz # shift
      for rowLoc = 1:Ny
        rowval[idx] = rowGbl[rowLoc]
        nzval[idx] = D2jl[colLoc,rowLoc]
        print(" ", rowGbl[rowLoc])
        print("(", rowLoc, ", ", colLoc, ")")
      end
    end
  end

  ix = 2
  print("\n\nix = ", ix)
  rowGbl_init[1] = 1 + (ix-1)*Ny*Nz
  for i = 2:Ny
    rowGbl_init[i] = rowGbl_init[i-1] + Nz
  end
  #print(rowGbl_init)

for colLoc = 1:Ny
  print("\n")
  colGbl_start = (colLoc-1)*Nz + 1 + (ix-1)*Ny*Nz
  colGbl_stop = colLoc*Nz + (ix-1)*Ny*Nz
  for colGbl = colGbl_start:colGbl_stop
    print("\ncolGbl: ", colGbl, " rowGbl: ")
    rowGbl[:] = rowGbl_init[:] + colGbl - 1 - (colLoc-1)*Nz - (ix-1)*Ny*Nz # shift
    for rowLoc = 1:Ny
      rowval[idx] = rowGbl[rowLoc]
      nzval[idx] = D2jl[colLoc,rowLoc]
      print(" ", rowGbl[rowLoc])
      print("(", rowLoc, ", ", colLoc, ")")
    end
  end
end

  #end # Nx

  print("\n")
  exit()

  return SparseMatrixCSC( Npoints, Npoints, colptr, rowval, nzval )
end


test_main([2,3,4])
