# Quadratic line minimization as suggested in original Beigi-Arias paper

function linmin1( pw, Vpot, α_t, psi, d, E, g )
    # compute directional derivative
    Δ = 2.0*real( sum( conj(d).*g ) )
    #
    psi_t = ortho_gram_schmidt( psi + α_t*d )
    E_trial = calc_Etot( pw, Vpot, psi_t )
    curvature = ( E_trial - ( E + α_t*Δ ) ) /α_t^2
    α = -Δ/(2*curvature)

    if abs(curvature) < eps()
        α = 0.0
    end
  
    return abs(α)   # need to use absolute value of α !!!
end
