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

  #@time Lx_v1 = build_nabla2_x(LF)
  build_nabla2_x(LF)

  #@time Lx_v2 = LF.LFx.D2jl ⊗ speye(Ny) ⊗ speye(Nz)
  #println(Lx_v1 - Lx_v2)
  #println(Lx_v1)

end

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
      append!( colptr, [colptr[end] + Nx])
    end
  end

  @printf("Converting to sparse matrix ..\n");
  println(size(colptr))
  println(rowval[2:end])
  #return sparse( rows[2:end], cols[2:end], vals[2:end] )
  #return SparseMatrixCSC( Npoints, Npoints, colptr, rowval[2:end], nzval[2:end] )
  exit()
end


test_main([3,3,3])
