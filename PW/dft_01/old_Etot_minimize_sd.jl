function Etot_minimize_sd( PW::PWGrid, psi::Array{Complex128,2},
                   Vpot::Array{Complex128,1}, Niter )
  α = 3e-5;
  Etot_old = 0.0
  Etot = 0.0

  for iter = 1:Niter
    psi = psi - (α*gradE( pw_grid, psi, Vpot ))
    Etot = get_Etot( pw_grid, psi, Vpot )
    conv = abs(Etot-Etot_old)
    @printf("Steepest descent: (s, step, conv) %8d %20.16f %20.16e\n",
            iter, Etot, conv )
    if conv < 1e-6
      print("Convergence achieved\n")
      break
    end
    Etot_old = Etot
  end
  return psi, Etot
end
