push!(LOAD_PATH, "../LF_common/")
using m_LF3d

include("../utils/ortho_gram_schmidt.jl")
include("orthonormalize.jl")
include("test_ortho.jl")

function test_main()
  NN = [25; 25; 25]
  AA = [0.0; 0.0; 0.0]
  BB = [6.0; 6.0; 6.0]

  Npoints = prod(NN)

  # Initialize LF
  LF = init_LF3d_c( NN, AA, BB )
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h

  srand(1234)
  v = rand(Npoints,4)
  orthonormalize!(LF,v)
  test_ortho( ΔV, 4, v )
end

test_main()
