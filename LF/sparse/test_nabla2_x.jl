push!(LOAD_PATH, "../LF_common/")
using m_LF3d
include("../LF_common/sparse_LF3d.jl")
include("build_nabla2_x.jl")

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
  @time Lx_v1 = build_nabla2_x(LF)
  #println(Lx_v1)

  @printf("Using Kronecker product directly")
  @time Lx_v2 = LF.LFx.D2jl ⊗ speye(Ny) ⊗ speye(Nz)

  # Check whether there are any differences between the two matrices
  println(Lx_v1 - Lx_v2)
  #println(Lx_v2)

end

test_main([65,65,65])
