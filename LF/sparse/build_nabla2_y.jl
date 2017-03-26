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

    rowGbl_init[1] = 1 + (ix-1)*Ny*Nz
    for i = 2:Ny
      rowGbl_init[i] = rowGbl_init[i-1] + Nz
    end

    for colLoc = 1:Ny
      colGbl_start = (colLoc-1)*Nz + 1 + (ix-1)*Ny*Nz
      colGbl_stop = colLoc*Nz + (ix-1)*Ny*Nz
      for colGbl = colGbl_start:colGbl_stop
        for rowLoc = 1:Ny
          idx = idx + 1
          rowval[idx] = rowGbl_init[rowLoc] + colGbl - 1 - (colLoc-1)*Nz - (ix-1)*Ny*Nz
          nzval[idx] = D2jl[colLoc,rowLoc]
        end
        idxCol = idxCol + 1
        colptr[idxCol] = colptr[idxCol-1] + Ny
      end
    end

  end # ix = 1:Nx

  return SparseMatrixCSC( Npoints, Npoints, colptr, rowval, nzval )
end
