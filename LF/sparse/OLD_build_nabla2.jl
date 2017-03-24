function build_nabla2_x(LF)
  D2jl = LF.LFx.D2jl
  Nx = LF.Nx
  Ny = LF.Ny
  Nz = LF.Nz
  Npoints = Nx*Ny*Nz

  # initialize with
  rowval = [0]  # start with dummy value 0
  colptr = [1]
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
