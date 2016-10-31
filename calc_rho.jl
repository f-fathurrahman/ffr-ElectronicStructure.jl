#
# Calculate electronic density in real space
#
function calc_rho( PW::PWGrid, Focc::Array{Float64,1}, psi::Array{Complex128,2} )
  Ω = PW.Ω
  Ns = PW.Ns
  Npoints = PW.Npoints
  Nstates = size(psi)[2]
  #
  ρ = zeros(Complex128,Npoints)
  psiR = G_to_R(Ns,psi)
  #
  ortho_gram_schmidt!(Nstates,psiR); scale!(sqrt(Npoints/Ω),psiR)
  #
  for is = 1:Nstates
    for ip = 1:Npoints
      ρ[ip] = ρ[ip] + conj(psiR[ip,is])*psiR[ip,is]*Focc[is]
    end
  end
  return real(ρ)
end

