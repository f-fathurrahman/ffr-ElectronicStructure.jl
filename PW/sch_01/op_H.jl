function op_H( pw_grid, Vpot, psi )
  return op_K( pw_grid, psi ) + op_Vpot( pw_grid, Vpot, psi )
end
