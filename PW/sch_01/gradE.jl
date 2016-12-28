function gradE( pw::PWGrid, Vpot, psi::Array{Complex128,2} )

  Npoints = size(psi)[1]
  Nstates = size(psi)[2]
  Ω = pw.Ω
  Ns = pw.Ns
  #
  grad = zeros( Complex128, Npoints, Nstates )

  H_psi = apply_H( pw, Vpot, psi )
  for i = 1:Nstates
    grad[:,i] = H_psi[:,i]
    for j = 1:Nstates
      grad[:,i] = grad[:,i] - dot( psi[:,j], H_psi[:,i] ) * psi[:,j]
    end
  end
  return grad

end
