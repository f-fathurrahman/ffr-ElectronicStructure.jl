push!(LOAD_PATH, "../LF_common/")
using m_LF3d
include("../LF_common/sparse_LF3d.jl")
include("build_nabla2_y.jl")

const ⊗ = kron

function test_main( NN::Array{Int64} )

  Nx = NN[1]
  Ny = NN[2]
  Nz = NN[3]

  AA = [0.0, 0.0, 0.0]
  BB = [16.0, 16.0, 16.0]

  # Initialize LF
  LF = init_LF3d_c( NN, AA, BB, verbose=true )

  @printf("Using loop")
  @time Ly_v1 = build_nabla2_y(LF)
  #println(Ly_v1)

  @printf("Using Kronecker product directly")
  @time Ly_v2 = speye(Nx) ⊗ LF.LFy.D2jl ⊗ speye(Nz)

  # Check whether there are any differences between the two matrices
  println(Ly_v1 - Ly_v2)

end


test_main([80,80,80])
