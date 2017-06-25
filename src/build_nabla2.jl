function build_nabla2(LF::LF3dGrid)

  Nx = LF.Nx
  Ny = LF.Ny
  Nz = LF.Nz
  Npoints = Nx*Ny*Nz

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

  ip_start = 0
  ip_stop = 0

  rowval_unsrt = zeros( Int64, nnzc )
  nzval_unsrt = zeros( Float64, nnzc )

  # rowLoc = 1:Nx, 1:Ny, 1:Nz for each gblCol
  for colGbl = 1:Npoints

    colLoc_x = ceil(Int,colGbl/(Ny*Nz))

    yy = colGbl - (colLoc_x-1)*Ny*Nz
    colLoc_y = ceil( Int, yy/Nz )

    izz = ceil(Int,colGbl/Nz)
    colLoc_z = colGbl - (izz-1)*Nz

    ip = 0  # reset ip

    # diagonal element
    ip = ip + 1
    rowval_unsrt[ip] = colGbl
    nzval_unsrt[ip] = D2jl_x[colLoc_x,colLoc_x] + D2jl_y[colLoc_y,colLoc_y] +  D2jl_z[colLoc_z,colLoc_z]

    # non diagonal element
    for ir = 1:Nx
      if ir != colLoc_x
        ip = ip + 1
        rowval_unsrt[ip] = rowGbl_x_orig[ir] + colGbl - 1 - (colLoc_x-1)*Nz*Ny
        nzval_unsrt[ip] = D2jl_x[ir,colLoc_x]
      end
    end

    for ir = 1:Ny
      if ir != colLoc_y
        ip = ip + 1
        rowval_unsrt[ip] = rowGbl_y_orig[ir] + colGbl - 1 - (izz-1)*Nz + (colLoc_x-1)*Ny*Nz
        nzval_unsrt[ip] = D2jl_y[ir,colLoc_y]
      end
    end

    for ir = 1:Nz
      if ir != colLoc_z
        ip = ip + 1
        rowval_unsrt[ip] = rowGbl_z_orig[ir] + (izz-1)*Nz
        nzval_unsrt[ip] = D2jl_z[ir,colLoc_z]
      end
    end
    ip_start = (colGbl-1)*nnzc + 1
    ip_stop = ip_start + nnzc - 1

    idx_srt = sortperm(rowval_unsrt)

    rowval[ip_start:ip_stop] = rowval_unsrt[idx_srt]
    nzval[ip_start:ip_stop] = nzval_unsrt[idx_srt]

  end

  # colptr has more predictable pattern
  colptr[1] = 1
  for i = 2:Npoints+1
    colptr[i] = colptr[i-1] + nnzc
  end

  nabla2 = SparseMatrixCSC( Npoints, Npoints, colptr, rowval, nzval )

end
