function apply_H( pw, Vpot, psi )
  return apply_K( pw, psi ) + apply_Vpot( pw, Vpot, psi )
end
