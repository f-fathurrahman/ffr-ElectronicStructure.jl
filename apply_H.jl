function apply_H( PW::PWGrid, V_loc, psi )
  return apply_K( PW, psi ) + apply_V_loc( PW, V_loc, psi )
end

"""
XXX: Recalculate V_loc ?
XXX: I think apply_H should not modify the Hamiltonian itself, i.e.
XXX: the potentials should not be recalculated here. Even if the potentials
XXX: need to be modified, it should be done outside this function.
"""
function apply_H( PW::PWGrid, Potentials::PotentialsT, Focc, psi::Array{Complex128,2} )
  rho = calc_rho( PW, Focc, psi )
  V_Hartree = real( G_to_R( PW.Ns, solve_poisson( PW, rho ) ) )
  V_xc = excVWN( rho ) + rho .* excpVWN( rho )
  #
  V_loc = Potentials.Ionic + V_Hartree + V_xc
  return apply_K( PW, psi ) + apply_V_loc( PW, V_loc, psi )
end
