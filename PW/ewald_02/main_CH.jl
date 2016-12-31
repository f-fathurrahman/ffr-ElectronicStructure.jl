include("../common/PWGrid_v01.jl")
include("../common/wrappers_fft.jl")
include("gen_dr_center.jl")
include("calc_strfact.jl")
include("calc_ewald.jl")

function test_main()

  Ns = [64,64,64]
  LatVecs = diagm( [16.0,16.0,16.0] )
  pw = PWGrid(Ns,LatVecs)

  # Atomic positions and nuclear charge
  Xpos = reshape( [0.0, 0.0, 0.0,
                   1.5, 0.0, 0.0], (3,2) )

  Nspecies = 2
  atmsymb = ["H", "C"] # unique list of atomic symbols
  atm2species = [1, 2]  # mapping from atom to species
  Zv = [1.0, 4.0]  # only valence ?

  # Calculate structure factor
  Sf = calc_strfact( Xpos, Nspecies, atm2species, pw.G )

  print("sum(Sf)="); println(sum(Sf))

  E_nn = calc_ewald( pw, Sf, Xpos, Nspecies, atm2species, Zv, sigma=[1.0,1.2] )
  @printf("E_nn = %18.10f\n", E_nn)

  E_nn = calc_ewald( pw, Sf, Xpos, Nspecies, atm2species, Zv, sigma=[0.3,1.2] )
  @printf("E_nn = %18.10f\n", E_nn)

end

@time test_main()
