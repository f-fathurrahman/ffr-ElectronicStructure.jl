function kssolve_Emin_cg( pw::PWGrid, V_ionic, Focc, Nstates::Int;
                         psi0=nothing, Potentials0 = nothing,
                         α_t = 3e-5, NiterMax=1000, verbose=false )

  Ns = pw.Ns
  Npoints = prod(Ns)

  if psi0 == nothing
    srand(2222)
    psi = randn(Npoints,Nstates) + im*randn(Npoints,Nstates)
    psi = ortho_gram_schmidt(psi)
  else
    psi = copy(psi0)
  end

  if Potentials0 == nothing
    Potentials = PotentialsT( V_ionic, zeros(Npoints), zeros(Npoints) )
    rho = calc_rho( pw, Focc, psi )
    Potentials.Hartree = real( G_to_R( Ns, solve_poisson(pw, rho) ) )
    Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )
  else
    Potentials = PotentialsT( Potentials0.Ionic,
                              Potentials0.Hartree,
                              Potentials0.XC )
  end

  d = zeros(Complex128, Npoints, Nstates)
  g_old = zeros(Complex128, Npoints, Nstates)
  d_old = zeros(Complex128, Npoints, Nstates)

  try_Potentials = PotentialsT( V_ionic, zeros(Npoints), zeros(Npoints) )

  beta     = 0.0
  Etot_old = 0.0
  Etot     = 0.0
  Energies = EnergiesT( 0.0, 0.0, 0.0, 0.0, 0.0 )

  for iter = 1:NiterMax

    g = gradE( pw, Potentials, Focc, psi)
    nrm = 0.0
    for is = 1:Nstates
      nrm = nrm + real( dot( g[:,is], g[:,is] ) )
    end
    if iter != 1
      beta = real(trace(g'*Kprec(pw,g)))/real(trace(g_old'*Kprec(pw,g_old)))
      #beta = real(trace((g-g_old)'*Kprec(pw,g)))/real(trace(g_old'*Kprec(pw,g_old)))
      #beta = real(trace((g-g_old)'*Kprec(pw,g)))/real(trace((g-g_old)'*d))
      #@printf("\nbeta = %f\n", beta)
    end

    d = -Kprec(pw, g) + beta * d_old

    psic = ortho_gram_schmidt(psi + α_t*d)
    rho = calc_rho( pw, Focc, psic )
    try_Potentials.Hartree = real( G_to_R( Ns, solve_poisson(pw, rho) ) )
    try_Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )
    gt = gradE( pw, try_Potentials, Focc, psic )

    if real(trace((g-gt)'*d)) != 0.0
      α = abs( α_t*real(trace(g'*d))/real(trace((g-gt)'*d)))
    else
      α = 0.0
    end

    # Update wavefunction
    psi = psi[:,:] + α*d[:,:]

    # Update potentials
    psi = ortho_gram_schmidt(psi)
    rho = calc_rho( pw, Focc, psi )
    Potentials.Hartree = real( G_to_R( Ns, solve_poisson(pw, rho) ) )
    Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )

    Energies = calc_Energies( pw, Potentials, Focc, psi )
    Etot = Energies.Total

    diff = abs(Etot-Etot_old)
    @printf("E step %8d = %18.10f %18.10f %18.10f\n", iter, Etot, diff, nrm/Nstates)
    if diff < 1e-6
      @printf("CONVERGENCE ACHIEVED\n")
      break
    end

    g_old = copy(g)
    d_old = copy(d)
    Etot_old = Etot
  end
  return psi, Energies, Potentials
  #
end
