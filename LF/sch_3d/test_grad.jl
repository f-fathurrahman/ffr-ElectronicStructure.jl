push!(LOAD_PATH, "../LF_common/")
using m_LF3d

include("init_pot_Hcoul.jl")
include("init_pot_harm_3d.jl")
include("../utils/ortho_gram_schmidt.jl")
include("../utils/orthonormalize.jl")
include("EnergiesT.jl")
include("op_nabla2.jl")
include("op_H.jl")
include("calc_Energies.jl")
include("gradE.jl")

function test_grad()
  # LF parameters
  NN = [25, 25, 25]
  AA = [0.0, 0.0, 0.0]
  BB = [6.0, 6.0, 6.0]

  Npoints = prod(NN)

  # Initialize LF
  LF = init_LF3d_p( NN, AA, BB )
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h

  # Parameter for potential
  center = AA + 0.5*(BB-AA)
  # Potential
  ω = 2.0
  Vpot = init_pot_harm_3d( LF, ω, center )

  Ncols = 4

  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h

  println("sum(Vpot) = ", sum(Vpot))
  println("ΔV = ", ΔV)

  v = zeros(Npoints,Ncols)
  for i = 1:Ncols
    v[i,i] = 1.0/sqrt(ΔV)
  end

  grad = gradE( LF, Vpot, v )
  println(sum(grad))

  Energies = calc_Energies( LF, Vpot, v)
  println(Energies)
end

test_grad()
