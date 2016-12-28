"""
An implementation of simple method to calculate Ewald energy
"""
function calc_ewald( pw_grid::PWGrid, Xpos, Sf; sigma=0.25 )
  #
  const Npoints = pw_grid.Npoints
  const Ω  = pw_grid.Ω
  const R  = pw_grid.R
  const Ns = pw_grid.Ns
  const G2 = pw_grid.G2
  #
  # Generate array of distances
  center = sum(pw_grid.LatVecs,2)/2
  dr = gen_dr( R, center )
  #
  # Generate charge density
  rho = gen_rho( Ns, dr, sigma, Sf )
  #intrho = sum(rho)*Ω/Npoints
  #@printf("sigma, int_rho: %4.1f %20.16f\n", sigma, intrho)
  #
  # Solve Poisson equation and calculate Hartree energy
  ctmp = 4.0*pi*R_to_G( Ns, rho )
  ctmp[1] = 0.0
  for ip = 2:Npoints
    ctmp[ip] = ctmp[ip] / G2[ip]
  end
  phi = real( G_to_R( Ns, ctmp ) )
  Ehartree = 0.5*dot( phi, rho ) * Ω/Npoints
  #
  Eself = 1.0/(2*sqrt(pi))*(1.0/sigma)*size(Xpos,2)
  E_nn = Ehartree - Eself
  #@printf("Ehartree, Eself, E_nn = %20.16f %20.16f %20.16f\n", Ehartree, Eself, E_nn)
  return E_nn
end
