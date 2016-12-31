include("../common/PWGrid_v01.jl")
include("../common/wrappers_fft.jl")

include("structure_factor.jl")
include("gen_dr.jl")
include("gen_rho.jl")
include("calc_ewald.jl")

function test_main()
  #
  #const Ns = [64, 64, 64]
  const Ns = [80, 80, 80]
  const LatVecs = diagm( [16.0, 16.0, 16.0] )
  #
  pw = PWGrid( Ns, LatVecs )
  #
  const Npoints = pw.Npoints
  const Ω = pw.Ω
  const R = pw.R
  const G = pw.G
  const G2 = pw.G2
  const Ns = pw.Ns

  @printf("dr, dVol = %10.5e %10.5e\n", (Ω/Npoints)^(1./3.), Ω/Npoints)

  # Atomic positions and nuclear charge
  Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )

  println("Xpos =")
  println(Xpos)

  Z = 1.0
  # Calculate structure factor
  Sf = structure_factor( Xpos, G )

  print("sum(Sf)="); println(sum(Sf))

  E_nn = calc_ewald( pw, Xpos, Sf, sigma=0.3 )
  @printf("E_nn = %18.10f\n", E_nn)
end

@time test_main()
