function solve_poisson( pw::PWGrid, rho )
  #
  Ω  = pw.Ω
  G2 = pw.gvectors.G2
  Ns = pw.Ns
  Ng = pw.gvectors.Ng
  idx_g2r = pw.gvectors.idx_g2r
  #
  ctmp = 4.0*pi*R_to_G( Ns, rho )
  for ig = 2:Ng
    ip = idx_g2r[ig]
    ctmp[ip] = ctmp[ip] / G2[ig]  # BEWARE: ig vs ip
  end
  ctmp[1] = 0.0
  phi = real( G_to_R( Ns, ctmp ) )
  return phi
end


# The same as `solve_poisson`, but using FFT that calls to C directly
function solve_poisson_new( pw_grid::PWGrid, rho )
  #
  Ω  = pw_grid.Ω
  G2 = pw_grid.G2
  Ns = pw_grid.Ns
  Npoints = pw_grid.Npoints
  #
  ctmp = 4.0*pi*c_R_to_G( Ns, convert(Array{Complex128},rho) )
  for ip = 2:Npoints
    ctmp[ip] = ctmp[ip] / G2[ip]
  end
  ctmp[1] = 0.0
  phi = real( c_G_to_R( Ns, ctmp ) )
  return phi
end
