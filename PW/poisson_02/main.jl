include("../pwgrid_03/PWGrid_v02.jl")
include("../common/wrappers_fft_v01.jl")

include("gen_dr.jl")
include("gen_rho.jl")
include("solve_poisson.jl")

function test_main()

  const ecutwfc = 35 * 0.5 # in Ha
  const LatVecs = diagm( [16.0, 16.0, 16.0] )
  #
  pw = PWGrid( ecutwfc, LatVecs )

  #
  Ns = pw.Ns
  Npoints = prod(Ns)
  Ω = pw.Ω
  R = pw.R

  @printf("Number of G-vectors: %d\n", pw.gvectors.Ng)
  @printf("Number of real space sampling points: %d\n", Npoints)
  @printf("Percentage G/R = %f\n", 100.0*pw.gvectors.Ng/Npoints)

  #
  # Generate array of distances
  #
  center = sum(LatVecs,2)/2
  dr = gen_dr( R, center )
  #
  # Generate charge density
  #
  const σ1 = 0.75
  const σ2 = 0.50
  rho = gen_rho( dr, σ1, σ2 )
  #
  # Solve Poisson equation and calculate Hartree energy
  #
  phi = solve_poisson( pw, rho )
  Ehartree = 0.5*dot( phi, rho ) * Ω/Npoints
  #
  Uanal = ( (1/σ1 + 1/σ2)/2 - sqrt(2) / sqrt( σ1^2 + σ2^2 ) ) / sqrt(pi)
  @printf("Num, ana, diff = %18.10f %18.10f %18.10e\n", Ehartree, Uanal, abs(Ehartree-Uanal))
end

#@code_native test_main()
test_main()
