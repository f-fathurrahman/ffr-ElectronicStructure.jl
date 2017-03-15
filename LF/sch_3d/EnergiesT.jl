type EnergiesT
  Total::Float64
  Kinetic::Float64
  Ionic::Float64
end

function print_Energies( Energies::EnergiesT )
  @printf("Total   energy: %18.10f\n", Energies.Total )
  @printf("Kinetic energy: %18.10f\n", Energies.Kinetic )
  @printf("Ionic   energy: %18.10f\n", Energies.Ionic )
end
