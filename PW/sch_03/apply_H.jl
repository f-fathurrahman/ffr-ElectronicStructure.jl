function apply_H( pw, Vpot, psi, ik )
  return apply_K( pw, psi, ik ) + apply_Vpot( pw, Vpot, psi, ik )
end
