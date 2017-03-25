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
