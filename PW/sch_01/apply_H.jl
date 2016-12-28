function apply_H( pw_grid, Vpot, psi )
  return apply_K( pw_grid, psi ) + apply_Vpot( pw_grid, Vpot, psi )
end
