# psi is assumed to be already orthonormalized elsewhere
function calc_Energies!( PW::PWGrid, Energies, Potentials, Focc::Array{Float64}, 
                         psi::Array{Complex128,2} )
  Ω = PW.Ω
  Ns = PW.Ns
  Npoints = size(psi)[1]
  Nstates = size(psi)[2]
  #
  Kpsi = apply_K( PW, psi )
  E_kin = 0.0
  for is = 1:Nstates
    E_kin = E_kin + real( Focc[is]*dot( psi[:,is], Kpsi[:,is] ) )
  end
  #println("sum Kpsi = $(sum(Kpsi))")
  #println("E_kin = $(E_kin)")
  #
  # Compute rho
  rho = calc_rho( PW, Focc, psi )
  #println("sum(rho)*Ω/Npoints = $(sum(rho)*Ω/Npoints)");
  V_Hartree = real(G_to_R( Ns, solve_poisson( PW, rho ) )) # 4*pi factor is handled in solve_poisson
  E_Hartree = 0.5*dot( V_Hartree, rho ) * Ω/Npoints
  #println("Hartree energy: $(E_Hartree)")
  #
  E_xc = dot( excVWN(rho), rho ) * Ω/Npoints
  #println("XC energy: $(E_xc)")
  #
  E_ionic = dot( Potentials.Ionic, rho ) * Ω/Npoints
  #println("Ionic energy: $(E_ionic)")
  #
  Energies.Kinetic = E_kin
  Energies.Hartree = E_Hartree
  Energies.XC = E_xc
  Energies.Ionic = E_ionic
  Energies.Total = E_kin + E_Hartree + E_xc + E_ionic

  return Energies.Total
end

