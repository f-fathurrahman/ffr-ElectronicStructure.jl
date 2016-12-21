include("PWGrid.jl")
include("R_to_G.jl")
include("G_to_R.jl")
include("gen_dr.jl")
include("gen_rho.jl")
include("solve_poisson.jl")

function test_main()
  #
  const Ns = [64, 64, 64]
  const LatVecs = diagm( [16.0, 16.0, 16.0] )
  #
  pw = PWGrid( Ns, LatVecs )
  #
  const Npoints = pw.Npoints
  const Ω = pw.Ω
  const R = pw.R
  const Ns = pw.Ns
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
