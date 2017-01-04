function calc_rho( pw::PWGrid, psi::Array{Complex128,2} )
  Ω = pw.Ω
  Ns = pw.Ns
  Npoints = prod(pw.Ns)
  Nstates = size(psi)[2]

  #
  ρ = zeros(Complex128,Npoints) # FIXME: or use Float64 instead?

  # Transform to real space
  cpsi = zeros( Complex128, Npoints, Nstates )
  cpsi[pw.gvecw.idx_gw2r,:] = psi[:,:]
  psiR = G_to_R(Ns, cpsi)

  # orthonormalization in real space
  ortho_gram_schmidt!(Nstates,psiR); scale!(sqrt(Npoints/Ω),psiR)
  for is = 1:Nstates
    for ip = 1:Npoints
      ρ[ip] = ρ[ip] + conj(psiR[ip,is])*psiR[ip,is]
    end
  end
  return real(ρ)
end
