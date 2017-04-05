# psi is assumed to be already orthonormalized elsewhere
# Potentials are not updated
function calc_Energies( PW::PWGrid, Potentials, Focc::Array{Float64},
                        psi::Array{Complex128,2} )
  Ω = PW.Ω
  Ns = PW.Ns
  Npoints = size(psi)[1]
  Nstates = size(psi)[2]
  #
  Kpsi = op_K( PW, psi )
  E_kin = 0.0
  for is = 1:Nstates
    E_kin = E_kin + Focc[is] * real( dot( psi[:,is], Kpsi[:,is] ) )
  end
  #
  # Compute rho
  rho = calc_rho( PW, Focc, psi )
  V_Hartree = real(G_to_R( Ns, Poisson_solve( PW, rho ) ))
  E_Hartree = 0.5*dot( V_Hartree, rho ) * Ω/Npoints
  #
  E_xc = dot( excVWN(rho), rho ) * Ω/Npoints
  #
  E_ionic = dot( Potentials.Ionic, rho ) * Ω/Npoints
  #
  Energies = EnergiesT(0.0, 0.0, 0.0, 0.0, 0.0)
  Energies.Kinetic = E_kin
  Energies.Ionic   = E_ionic
  Energies.Hartree = E_Hartree
  Energies.XC      = E_xc
  Energies.Total   = E_kin + E_ionic + E_Hartree + E_xc
  #
  return Energies

end
