push!(LOAD_PATH, "../LF_common/")
using m_LF3d
include("../LF_common/sparse_LF3d.jl")
include("build_nabla2_x.jl")
include("build_nabla2_y.jl")
include("build_nabla2_z.jl")

const ⊗ = kron

function test_main( NN::Array{Int64} )

  println("NN = ", NN)

  Nx = NN[1]
  Ny = NN[2]
  Nz = NN[3]

  AA = [0.0, 0.0, 0.0]
  BB = [16.0, 16.0, 16.0]

  # Initialize LF
  LF = init_LF3d_c( NN, AA, BB, verbose=true )

  print("Building nabla2_x")
  @time Lx_v1 = build_nabla2_x(LF)

  print("Building nabla2_y")
  @time Ly_v1 = build_nabla2_y(LF)

  print("Building nabla2_z")
  @time Lz_v1 = build_nabla2_z(LF)

  #print("Forming Laplacian:")
  #@time ∇2 = Lx_v1 + Ly_v1

  #print("\n")

  #print("Using Kronecker product directly: nabla2_x")
  #@time Lx_v2 = LF.LFx.D2jl ⊗ speye(Ny) ⊗ speye(Nz)

  #print("Using Kronecker product directly: nabla2_y")
  #@time Ly_v2 = speye(Nx) ⊗ LF.LFy.D2jl ⊗ speye(Nz)

  #@printf("Using Kronecker product directly: nabla2_z")
  #@time Lz_v2 = speye(Nx) ⊗ speye(Ny) ⊗ LF.LFz.D2jl

  #print("Forming Laplacian:")
  #@time ∇2 = Lx_v2 + Ly_v2

end

@code_native test_main([3,3,3])

test_main( [80, 80, 80])
