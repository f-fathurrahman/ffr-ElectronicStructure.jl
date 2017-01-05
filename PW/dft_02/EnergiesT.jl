type EnergiesT
  Total::Float64
  Kinetic::Float64
  Ionic::Float64
  Hartree::Float64
  XC::Float64
end

function print_Energies( Energies::EnergiesT )
  @printf("Total   energy: %18.10f\n", Energies.Total )
  @printf("Kinetic energy: %18.10f\n", Energies.Kinetic )
  @printf("Ionic   energy: %18.10f\n", Energies.Ionic )
  @printf("Hartree energy: %18.10f\n", Energies.Hartree )
  @printf("XC      energy: %18.10f\n", Energies.XC )
end
