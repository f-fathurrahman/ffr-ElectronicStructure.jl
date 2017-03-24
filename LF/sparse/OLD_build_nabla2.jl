function build_nabla2_x(LF)
  D2jl = LF.LFx.D2jl
  Nx = LF.Nx
  Ny = LF.Ny
  Nz = LF.Nz
  Npoints = Nx*Ny*Nz

  # initialize with
  rowval = [0]  # start with dummy value 0
  colptr = zeros(Int64,Npoints+1)
  colptr[1] = 1
  nzval = [0.0]  # start with dummy value 0.0

  colLoc = 1
  for colGbl = 1:Ny*Nz
    for rowLoc = 1:Nx
      rowGbl = (rowLoc-1)*Ny*Nz + colGbl - 0*Ny*Nz
      #println(colLoc, " ", rowGbl)
      append!( rowval, [rowGbl] )
      append!( nzval, [D2jl[rowLoc,colLoc]] )
    end
    colptr[colGbl+1] = colptr[colGbl] + Nx
  end

  #println("")
  colLoc = 2
  colGblStart = (colLoc-1)*(Ny*Nz+1)
  colGblStop = colLoc*Ny*Nz
  for colGbl =Ny*Nz+1:2*Ny*Nz
    for rowLoc = 1:Nx
      rowGbl = (rowLoc-1)*Ny*Nz + colGbl - 1*Ny*Nz
      #println(colLoc, " ", rowGbl)
      append!( rowval, [rowGbl] )
      append!( nzval, [D2jl[rowLoc,colLoc]] )
    end
    colptr[colGbl+1] = colptr[colGbl] + Nx
  end

  #println("")
  colLoc = 3
  for colGbl = 2*Ny*Nz+1:3*Ny*Nz
    for rowLoc = 1:Nx
      rowGbl = (rowLoc-1)*Ny*Nz + colGbl - 2*Ny*Nz
      #println(colLoc, " ", rowGbl)
      append!( rowval, [rowGbl] )
      append!( nzval, [D2jl[rowLoc,colLoc]] )
    end
    colptr[colGbl+1] = colptr[colGbl] + Nx
  end

  println(size(rowval[2:end]))
  println(size(nzval[2:end]))
  println(size(colptr))
  println(colptr)
  return SparseMatrixCSC( Npoints, Npoints, colptr, rowval[2:end], nzval[2:end] )
end

function build_nabla2_x(LF)
  D2jl = LF.LFx.D2jl
  Nx = LF.Nx
  Ny = LF.Ny
  Nz = LF.Nz
  Npoints = Nx*Ny*Nz

  # initialize with
  rowval = [0]  # start with dummy value 0
  colptr = zeros(Int,Npoints)
  colptr[1] = 1
  nzval = [0.0]  # start with dummy value 0.0

  for jx=1:Nx
    j = (jx-1)*Ny*Nz + 1
    nnzc = 0
    for ix=1:Nx
      i = (ix-1)*Ny*Nz + 1
      for dd = 1:Ny*Nz
        nnzc = nnzc + 1
        append!( rowval, [i+dd-1])
        append!( nzval, [D2jl[ix,jx]] )
      end
      append!( colptr, [colptr[end] + nnzc])
    end
  end

  @printf("Converting to sparse matrix ..\n");
  println(size(colptr))
  println(rowval[2:end])
  #return sparse( rows[2:end], cols[2:end], vals[2:end] )
  #return SparseMatrixCSC( Npoints, Npoints, colptr, rowval[2:end], nzval[2:end] )
  exit()
end
