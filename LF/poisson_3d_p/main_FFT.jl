include("../LF_common/m_LF1d.jl")
include("../LF_common/m_LF3d.jl")
include("../LF_common/m_Gvectors.jl")

include("../LF_common/apply_Laplacian.jl")
include("../LF_common/solve_poisson_cg.jl")
include("linsolve_cg_v2.jl")
include("linsolve_cr.jl")
include("solve_poisson_FFT.jl")

function test_main( NN::Array{Int64} )

  AA = [0.0, 0.0, 0.0]
  BB = [16.0, 16.0, 16.0]

  # Initialize LF
  LF = init_LF3d_p( NN, AA, BB, verbose=true )

  # Box dimensions
  Lx = BB[1] - AA[1]
  Ly = BB[2] - AA[2]
  Lz = BB[3] - AA[3]

  # Center of the box
  x0 = Lx/2.0
  y0 = Ly/2.0
  z0 = Lz/2.0

  # Parameters for two gaussian functions
  sigma1 = 0.75
  sigma2 = 0.50
  Npoints = LF.Nx * LF.Ny * LF.Nz

  rho = zeros( Npoints )
  phi = zeros( Npoints )

  # Initialization of charge density
  for ip = 1 : Npoints
    r = norm( LF.lingrid[:,ip] - [ x0, y0, z0 ] )
    rho[ip] = exp( -r^2 / (2.0*sigma2^2) ) / (2.0*pi*sigma2^2)^1.5 -
              exp( -r^2 / (2.0*sigma1^2) ) / (2.0*pi*sigma1^2)^1.5
  end

  deltaV = LF.LFx.h * LF.LFy.h * LF.LFz.h

  @printf("#Test norm charge: %f\n", sum(rho)*deltaV)
  @printf("#Solving Poisson equation:\n")
  #phi = solve_poisson_cg( LF, -4.0*pi*rho, 1000, verbose=true, TOL=5e-9 )

  phi = zeros(Npoints)
  linsolve_cg_v2!( LF, phi, -4.0*pi*rho, verbose=true, TOL=5e-9 )
  #linsolve_cr!( LF, phi, -4.0*pi*rho, verbose=true, TOL=5e-9 )

  #phi = linsolve_bicgstab( LF, -4.0*pi*rho, 200, verbose=true )

  Gvec    = GvectorsT( NN, diagm( [Lx,Ly,Lz] ) )
  phi_FFT = solve_poisson_FFT( Gvec, rho )

  @printf("LF ΔV = %18.10f\n", deltaV)
  @printf("Gv ΔV = %18.10f\n", Gvec.Ω/prod(NN))

  # Calculation of Hartree energy
  Unum = 0.5*sum( rho .* phi ) * deltaV
  Unum_FFT = 0.5*sum( rho .* phi_FFT ) * deltaV
  Uana = ( ( 1.0/sigma1 + 1.0/sigma2 ) / 2.0 - sqrt(2.0)/sqrt(sigma1^2 + sigma2^2) ) / sqrt(pi)
  @printf("#CG, analytic, diff: %20.10f %20.10f %20.10f\n", Unum, Uana, abs(Unum-Uana))
  @printf("#PW, analytic, diff: %20.10f %20.10f %20.10f\n", Unum_FFT, Uana, abs(Unum_FFT-Uana))
end

#@code_native test_main([3,3,3])
@time test_main([25,25,25])
