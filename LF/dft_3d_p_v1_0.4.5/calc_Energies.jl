"""
Calculate total energy components.
Potentials are NOT updated.
psi is assumed to be already orthogonalized properly
rho is calculated outside this function.
"""
function calc_Energies( LF, ∇2, Potentials::PotentialsT, rho, Focc, psi::Array{Float64,2} )
  #
  Nstates = size(psi)[2]
  ΔV  = LF.LFx.h * LF.LFy.h * LF.LFz.h
  #
  E_Hartree = 0.5*sum( rho .* Potentials.Hartree ) * ΔV
  #
  E_xc = sum(excVWN( rho ) .* rho) * ΔV
  #
  E_ionic = sum( rho .* Potentials.Ionic ) * ΔV

  #
  E_kin = 0.0
  for ist = 1:Nstates
    nabla2_v = ∇2 * psi[:,ist]
    E_kin = E_kin + -0.5*dot( psi[:,ist], nabla2_v ) * ΔV * Focc[ist]
  end

  Energies = EnergiesT(0.0, 0.0, 0.0, 0.0, 0.0)
  Energies.Kinetic = E_kin
  Energies.Ionic   = E_ionic
  Energies.Hartree = E_Hartree
  Energies.XC      = E_xc
  Energies.Total   = E_kin + E_ionic + E_Hartree + E_xc
  #
  return Energies
end


"""
Calculate total energy components AND update Potentials.
psi is assumed to be already orthogonalized properly
This function modifies `Energies` and `Potentials`.
"""
function calc_Energies!( LF, Gv, Energies, Potentials, Focc, psi::Array{Float64,2} )

  Npoints = size(psi)[1]
  Nstates = size(psi)[2]
  ΔV  = LF.LFx.h * LF.LFy.h * LF.LFz.h

  # Electron density
  rho = calc_rho( Focc, psi )
  #
  # solve Poisson equation to get Hartree potential
  V_Hartree = solve_poisson_FFT( Gv, rho )
  # Calculate Hartree energy
  E_Hartree = 0.5*sum( rho .* V_Hartree ) * ΔV
  #
  V_xc = excVWN( rho ) + rho .* excpVWN( rho )
  E_xc = sum(excVWN( rho ) .* rho) * ΔV
  #
  E_ionic = sum( rho .* Potentials.Ionic ) * ΔV

  #
  E_kin = 0.0
  for ist = 1:Nstates
    nabla2_v = op_nabla2( LF, psi[:,ist] )
    E_kin = E_kin + -0.5*dot( psi[:,ist], nabla2_v ) * ΔV * Focc[ist]
  end

  #
  Energies.Kinetic = E_kin
  Energies.Ionic   = E_ionic
  Energies.Hartree = E_Hartree
  Energies.XC      = E_xc
  Energies.Total   = E_kin + E_ionic + E_Hartree + E_xc
  #
  Potentials.Hartree = V_Hartree[:]
  Potentials.XC      = V_xc[:]

  return
end


function calc_Energies!( LF, Gv, ∇2, Energies, Potentials, Focc, psi::Array{Float64,2} )

  Npoints = size(psi)[1]
  Nstates = size(psi)[2]
  ΔV  = LF.LFx.h * LF.LFy.h * LF.LFz.h

  # Electron density
  rho = calc_rho( Focc, psi )
  #
  # solve Poisson equation to get Hartree potential
  V_Hartree = solve_poisson_FFT( Gv, rho )
  # Calculate Hartree energy
  E_Hartree = 0.5*sum( rho .* V_Hartree ) * ΔV
  #
  V_xc = excVWN( rho ) + rho .* excpVWN( rho )
  E_xc = sum(excVWN( rho ) .* rho) * ΔV
  #
  E_ionic = sum( rho .* Potentials.Ionic ) * ΔV

  #
  E_kin = 0.0
  for ist = 1:Nstates
    nabla2_v = ∇2 * psi[:,ist]
    E_kin = E_kin + -0.5*dot( psi[:,ist], nabla2_v ) * ΔV * Focc[ist]
  end

  #
  Energies.Kinetic = E_kin
  Energies.Ionic   = E_ionic
  Energies.Hartree = E_Hartree
  Energies.XC      = E_xc
  Energies.Total   = E_kin + E_ionic + E_Hartree + E_xc
  #
  Potentials.Hartree = V_Hartree[:]
  Potentials.XC      = V_xc[:]

  return
end
