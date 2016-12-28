function schsolve_Emin_sd( pw::PWGrid, Vpot, psi::Array{Complex128,2};
                           NiterMax=1000 )
  α = 3e-5
  Etot_old = 0.0
  Etot = 0.0

  for iter = 1:NiterMax
    psi = psi - α*gradE( pw, Vpot, psi )

    psi  = ortho_gram_schmidt(psi)
    Etot = calc_Etot( pw, Vpot, psi )

    conv = abs(Etot-Etot_old)
    @printf("Steepest descent: (s, step, conv) %8d %18.10f %18.10e\n",
            iter, Etot, conv )
    if conv < 1e-6
      print("Convergence achieved\n")
      break
    end
    Etot_old = Etot
  end
  return psi, Etot
end
