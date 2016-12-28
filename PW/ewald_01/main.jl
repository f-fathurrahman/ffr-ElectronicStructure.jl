include("../common/PWGrid_v01.jl")
include("../common/wrappers_fft.jl")

include("structure_factor.jl")
include("gen_dr.jl")
include("gen_rho.jl")
include("calc_ewald.jl")

function test_main()
  Ns = [64,64,64]
  LatVecs = diagm( [16.0,16.0,16.0] )
  pw_grid = PWGrid(Ns,LatVecs)
  # Atomic positions and nuclear charge
  Xpos = reshape( [0.0, 0.0, 0.0,
                   1.5, 0.0, 0.0], (3,2) )
  Z = 1.0
  # Calculate structure factor
  Sf = structure_factor( Xpos, pw_grid.G )

  E_nn = calc_ewald( pw_grid, Xpos, Sf )
  @printf("E_nn = %18.10f\n", E_nn)
end

@time test_main()
