function apply_H( pw, Potentials, psi )
  V_loc = Potentials.Ionic + Potentials.Hartree + Potentials.XC
  return apply_K( pw, psi ) + apply_V_loc( pw, V_loc, psi )
end
