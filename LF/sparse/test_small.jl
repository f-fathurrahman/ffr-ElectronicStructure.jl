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

  AA = [0.0, 0.0, 0.0]
  BB = [16.0, 16.0, 16.0]

  # Initialize LF
  LF = init_LF3d_p( NN, AA, BB, verbose=true )

  #print("\n")

  print("Using Kronecker product directly: nabla2_x")
  @time Lx_v2 = LF.LFx.D2jl ⊗ speye(Ny) ⊗ speye(Nz)

  print("Using Kronecker product directly: nabla2_y")
  @time Ly_v2 = speye(Nx) ⊗ LF.LFy.D2jl ⊗ speye(Nz)

  @printf("Using Kronecker product directly: nabla2_z")
  @time Lz_v2 = speye(Nx) ⊗ speye(Ny) ⊗ LF.LFz.D2jl

  print("Forming Laplacian:")
  @time ∇2 = Lx_v2 + Ly_v2

  #println(Lx_v2)
  #println(Ly_v2)
  println( Lx_v2 + Ly_v2 )

end

test_main( [5, 3, 5])

