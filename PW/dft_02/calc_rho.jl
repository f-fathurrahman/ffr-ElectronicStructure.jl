function calc_rho( pw::PWGrid, Focc, psi::Array{Complex128,2} )
  Ω  = pw.Ω
  Ns = pw.Ns
  Npoints = prod(Ns)
  Nstates = size(psi)[2]
  #
  rho = zeros(Complex128,Npoints)

  # Transform to real space
  cpsi = zeros( Complex128, Npoints, Nstates )
  cpsi[pw.gvecw.idx_gw2r,:] = psi[:,:]
  psiR = G_to_R(Ns, cpsi)

  # orthonormalization in real space
  ortho_gram_schmidt!( Nstates, psiR )
  scale!( sqrt(Npoints/Ω), psiR )

  for is = 1:Nstates
    for ip = 1:Npoints
      rho[ip] = rho[ip] + conj(psiR[ip,is])*psiR[ip,is]*Focc[is]
    end
  end
  return real(rho)
end
