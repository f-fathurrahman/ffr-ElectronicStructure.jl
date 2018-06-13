function op_H( pw, Vpot, psi, ik )
  return op_K( pw, psi, ik ) + op_Vpot( pw, Vpot, psi, ik )
end
