# Quadratic line minimization as suggested in original Beigi-Arias paper

function quad_lin_min( pw, Vpot, α_t, psi, d, E, g )
  # compute directional derivative
  Δ = 2.0*real( sum( conj(d).*g ) )
  #
  psi_t = ortho_gram_schmidt( psi + α_t*d )
  E_trial = calc_Etot( pw, Vpot, psi_t )
  curvature = ( E_trial - ( E + α_t*Δ ) ) /α_t^2
  α = -Δ/(2*curvature)

  #println("Etot = ", E)
  #println("E_trial = ", E_trial)
  #println("curvature = ", curvature)
  #println("Δ = ", Δ)
  #println("α = ", α)

  return abs(α)   # need to use absolute value of α !!!
end
