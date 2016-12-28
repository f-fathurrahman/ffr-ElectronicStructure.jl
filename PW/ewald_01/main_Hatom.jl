include("../common/PWGrid_v01.jl")
include("../common/wrappers_fft.jl")

include("structure_factor.jl")
include("gen_dr.jl")
include("gen_rho.jl")
include("calc_ewald.jl")

function test_main()
  #
  const Ns = [64, 64, 64]
  const LatVecs = diagm( [16.0, 16.0, 16.0] )
  #
  pw_grid = PWGrid( Ns, LatVecs )
  #
  const Npoints = pw_grid.Npoints
  const Ω = pw_grid.Ω
  const R = pw_grid.R
  const G = pw_grid.G
  const G2 = pw_grid.G2
  const Ns = pw_grid.Ns
  # Atomic positions and nuclear charge
  Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )
  Z = 1.0
  # Calculate structure factor
  Sf = structure_factor( Xpos, G )

  E_nn = calc_ewald( pw_grid, Xpos, Sf )
  @printf("E_nn = %18.10f\n", E_nn)
end

@time test_main()
