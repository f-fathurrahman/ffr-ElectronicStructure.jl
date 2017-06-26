#Version 2
function build_nabla2_x(LF)
  D2jl = LF.LFx.D2jl
  Nx = LF.Nx
  Ny = LF.Ny
  Nz = LF.Nz
  Npoints = Nx*Ny*Nz

  # initialize with
  rowval = zeros(Int64,Nx*Npoints)
  colptr = zeros(Int64,Npoints+1)
  colptr[1] = 1
  nzval = zeros(Float64,Nx*Npoints)

  idx = 1
  for colLoc = 1:Nx
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


# Version 1
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
