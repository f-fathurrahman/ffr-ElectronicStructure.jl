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

  @printf("Using Kronecker product directly")
  @time Ly_v2 = speye(Nx) ⊗ LF.LFx.D2jl ⊗ speye(Nz)

  # Check whether there are any differences between the two matrices
  println(Lx_v1 - Lx_v2)
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
  for colLoc = 1:Ny
    colGbl_start = (colLoc-1)*(Ny*Nz)+1
    colGbl_stop = colLoc*Ny*Nz
    #println("")
    #println("colGbl_start = ", colGbl_start)
    #println("colGbl_stop  = ", colGbl_stop)
    for colGbl = colGbl_start:colGbl_stop
      for rowLoc = 1:Nx
        rowGbl = (rowLoc-1)*Ny*Nz + colGbl - (colLoc-1)*Ny*Nz
        rowval[idx] = rowGbl
        nzval[idx] = D2jl[rowLoc,colLoc]
        idx = idx + 1
        #println(colLoc, " ", rowGbl)
        #append!( rowval, [rowGbl] )
        #append!( nzval, [D2jl[rowLoc,colLoc]] )
      end
      colptr[colGbl+1] = colptr[colGbl] + Nx
    end
  end # colLoc
  #@printf("Finished looping\n")


  #println(size(rowval[2:end]))
  #println(size(nzval[2:end]))
  #println(size(colptr))
  #println(colptr)
  return SparseMatrixCSC( Npoints, Npoints, colptr, rowval, nzval )
end


test_main([3,2,4])
