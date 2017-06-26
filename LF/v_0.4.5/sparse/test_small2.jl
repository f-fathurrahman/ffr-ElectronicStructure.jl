# For comparison with ffr-LFDFT

push!(LOAD_PATH, "../LF_common/")
using m_LF3d
include("../LF_common/sparse_LF3d.jl")

const ⊗ = kron

function test_main( NN::Array{Int64} )

  println("NN = ", NN)

  Nx = NN[1]
  Ny = NN[2]
  Nz = NN[3]
  Npoints = Nx*Ny*Nz

  AA = [0.0, 0.0, 0.0]
  BB = [16.0, 16.0, 16.0]

  # Initialize LF
  LF = init_LF3d_c( NN, AA, BB, verbose=true )

  ∇2 = get_Laplacian3d_kron( LF )

  nzval = ∇2.nzval
  rowval = ∇2.rowval
  for ip = 1:size(nzval)[1]
    @printf("%5d %5d %18.10f\n", ip, rowval[ip], nzval[ip])
  end

  v = 1.1*ones(Npoints)
  Av = ∇2 * v

  for ip = 1:Npoints
    @printf("%5d %18.10f %18.10f\n", ip, v[ip], Av[ip])
  end

end

test_main( [3, 2, 4])

