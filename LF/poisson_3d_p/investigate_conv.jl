push!(LOAD_PATH, "../LF_common/")
using m_LF3d

include("../LF_common/apply_Laplacian.jl")
include("../LF_common/solve_poisson_cg.jl")

function do_calc( NN::Array{Int64} )

  AA = [0.0, 0.0, 0.0]
  BB = [16.0, 16.0, 16.0]

  # Initialize LF
  LF = init_LF3d_c( NN, AA, BB )

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
  phi = solve_poisson_cg( LF, -4.0*pi*rho, 1000, verbose=true )

  # Calculation of Hartree energy
  Unum = 0.5*sum( rho .* phi ) * deltaV
  Uana = ( ( 1.0/sigma1 + 1.0/sigma2 ) / 2.0 - sqrt(2.0)/sqrt(sigma1^2 + sigma2^2) ) / sqrt(pi)
  @printf("#Numeric, analytic, diff: %20.10f %20.10f %20.10f\n", Unum, Uana, abs(Unum-Uana))
end

import PyPlot
const plt = PyPlot
using PyCall
@pyimport numpy
const np = numpy

function do_plot( N::Int64 )
  data = np.loadtxt( "log_"*string(N) )
  plt.clf()
  plt.plot( data[:,1], log10(data[:,2]), marker="o" )
  plt.savefig( "log_"*string(N)*".png", dpi=300 )
end

function do_plot( )
  plt.clf()
  for ii in [19, 25, 29, 35, 39, 45]
    data = np.loadtxt( "log_"*string(ii) )
    plt.plot( data[:,1], log10(data[:,2]), marker="o", label="N-"*string(ii) )
  end
  plt.legend()
  plt.title("Convergence of CG Iteration")
  plt.xlabel("Iteration")
  plt.ylabel("log(error)")
  plt.grid()
  plt.savefig( "log_all.png", dpi=300 )
  plt.savefig( "log_all.pdf" )
end

#do_calc( 45*ones(Int64,3) )
do_plot()
