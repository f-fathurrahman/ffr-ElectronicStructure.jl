"""
Input psi is assumed to be already orthonormalized
"""
function calc_grad( PW::PWGrid, Potentials, Focc, psi::Array{Complex128,2} )
  Npoints = size(psi)[1]
  Nstates = size(psi)[2]
  Ω = PW.Ω
  Ns = PW.Ns
  #
  grad = zeros( Complex128, Npoints, Nstates )
  #
  rho = calc_rho( PW, Focc, psi )
  V_Hartree = real( G_to_R( Ns, solve_poisson( PW, rho ) ) )
  V_xc = excVWN( rho ) + rho .* excpVWN( rho )
  V_loc = Potentials.Ionic + V_Hartree + V_xc

  H_psi = apply_H( PW, V_loc, psi )
  for i = 1:Nstates
    grad[:,i] = H_psi[:,i]
    for j = 1:Nstates
      grad[:,i] = grad[:,i] - dot( psi[:,j], H_psi[:,i] ) * psi[:,j]
    end
    grad[:,i] = Focc[i]*grad[:,i]
  end
  return grad
end

