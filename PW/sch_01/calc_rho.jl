function calc_rho( pw::PWGrid, psi::Array{Complex128,2} )
  Ω = pw.Ω
  Ns = pw.Ns
  Npoints = pw.Npoints
  Nstates = size(psi)[2]
  #
  ρ = zeros(Complex128,Npoints)
  # Transform to real space
  psiR = G_to_R(Ns,psi)
  # orthonormalization in real space
  ortho_gram_schmidt!(Nstates,psiR); scale!(sqrt(Npoints/Ω),psiR)
  for is = 1:Nstates
    for ip = 1:Npoints
      ρ[ip] = ρ[ip] + conj(psiR[ip,is])*psiR[ip,is]
    end
  end
  return real(ρ)
end
