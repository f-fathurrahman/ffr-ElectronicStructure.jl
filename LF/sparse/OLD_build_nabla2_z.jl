
# not efficient for large Nz
function build_nabla2_z(LF)

  D2jl = LF.LFz.D2jl

  Nx = LF.Nx
  Ny = LF.Ny
  Nz = LF.Nz

  Npoints = Nx*Ny*Nz

  # initialize with
  rowval = zeros(Int64,Nz*Npoints)
  colptr = zeros(Int64,Npoints+1)
  colptr[1] = 1
  nzval = zeros(Float64,Nz*Npoints)

  idx = 0
  idxCol = 1

  #print("\n")
  for ic = 1:Nx*Ny
    colGbl_start = (ic-1)*Nz + 1
    colGbl_stop = ic*Nz
    for colGbl = colGbl_start:colGbl_stop
      colLoc = colGbl - (ic-1)*Nz
      for rowLoc = 1:Nz
        idx = idx + 1
        rowGbl = rowLoc + (ic-1)*Nz
        rowval[idx] = rowGbl
        nzval[idx] = D2jl[rowLoc,colLoc]
        #println("colGbl = ", colGbl, " rowLoc = ", rowLoc, " colLoc = ", colLoc)
      end
      idxCol = idxCol + 1
      colptr[idxCol] = colptr[idxCol-1] + Nz
    end
  end

  return SparseMatrixCSC( Npoints, Npoints, colptr, rowval, nzval )

end
