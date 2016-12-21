function E_min_sd!( PW::PWGrid, Energies::EnergiesT, Potentials::PotentialsT, Focc,
                    psi0::Array{Complex128,2} ;
                    α_t=3e-5, NiterMax=500, verbose=false )
  #
  psi = copy(psi0)
  Ω = PW.Ω
  #
  Etot_old = 0.0
  Etot = 0.0

  for iter = 1:NiterMax

    ortho_gram_schmidt!( psi )

    psi = psi - ( α_t * calc_grad( PW, Potentials, Focc, psi ) )

    ortho_gram_schmidt!( psi )
    Etot = calc_Energies!( PW, Energies, Potentials, Focc, psi )

    conv = abs(Etot-Etot_old)
    if verbose
      @printf("E min SD: %8d %18.10f %18.10e\n",
              iter, Etot, conv )
    end
    if conv < 1e-6
      print("Convergence achieved\n")
      break
    end
    Etot_old = Etot
  end
  return Etot, psi
end
