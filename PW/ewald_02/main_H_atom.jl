include("../common/PWGrid_v01.jl")
include("../common/wrappers_fft.jl")
include("../common/gen_dr_center.jl")
include("../common/calc_strfact_v1.jl")
include("../common/calc_ewald_v1.jl")

function test_main()

  Ns = [64,64,64]
  LatVecs = diagm( [16.0,16.0,16.0] )
  pw = PWGrid(Ns,LatVecs)

  # Atomic positions and nuclear charge
  Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )

  Nspecies = 1
  atmsymb = ["H"] # unique list of atomic symbols
  atm2species = [1]  # mapping from atom to species
  Zv = [1.0]  # only valence ?

  # Calculate structure factor
  Sf = calc_strfact( Xpos, Nspecies, atm2species, pw.G )

  print("sum(Sf)="); println(sum(Sf))

  E_nn = calc_ewald( pw, Sf, Xpos, Nspecies, atm2species, Zv )
  @printf("E_nn = %18.10f\n", E_nn)
end

@time test_main()
