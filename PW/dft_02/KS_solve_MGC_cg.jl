function KS_solve_MGC_cg( pw::PWGrid, V_ionic, Focc, Nstates::Int;
                           psi0=nothing, Potentials0 = nothing,
                           α_t = 3e-5, NiterMax=1000, verbose=false )

  Ns = pw.Ns
  Npoints = prod(Ns)
  Ngwx = pw.gvecw.Ngwx

  if psi0 == nothing
    srand(2222)
    psi = randn(Ngwx,Nstates) + im*randn(Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)
  else
    psi = copy(psi0)
  end

  if Potentials0 == nothing
    Potentials = PotentialsT( V_ionic, zeros(Npoints), zeros(Npoints) )
    rho = calc_rho( pw, Focc, psi )
    Potentials.Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
    Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )
  else
    Potentials = PotentialsT( Potentials0.Ionic,
                              Potentials0.Hartree,
                              Potentials0.XC )
  end

  d = zeros(Complex128, Ngwx, Nstates)
  g_old = zeros(Complex128, Ngwx, Nstates)
  d_old = zeros(Complex128, Ngwx, Nstates)
  Kg = zeros(Complex128, Ngwx, Nstates)
  Kg_old = zeros(Complex128, Ngwx, Nstates)

  β        = 0.0
  Etot_old = 0.0
  Etot     = 0.0
  Energies = EnergiesT( 0.0, 0.0, 0.0, 0.0, 0.0 )

  for iter = 1:NiterMax

    g = calc_grad_MGC( pw, Potentials, Focc, psi)
    Kg = Kprec(pw,g)

    if iter != 1
      #β = real(sum(conj(g).*Kg))/real(sum(conj(g_old).*Kg_old))
      β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g_old).*Kg_old))
      #β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g-g_old).*d))
      #β = real(sum(conj(g).*Kg))/real(sum((g-g_old).*conj(d_old)))
    end

    d = -Kprec(pw, g) + β * d_old

    psic = ortho_gram_schmidt(psi + α_t*d)
    #psic = psi + α_t*d
    rho = calc_rho( pw, Focc, psic )
    Potentials.Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
    Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )
    gt = calc_grad_MGC( pw, Potentials, Focc, psic )

    denum = real(sum(conj(g-gt).*d))
    if denum != 0.0
      α = abs( α_t*real(sum(conj(g).*d))/denum )
    else
      α = 0.0
    end

    # Update wavefunction
    psi = psi[:,:] + α*d[:,:]

    # Update potentials
    psi = ortho_gram_schmidt(psi)
    rho = calc_rho( pw, Focc, psi )

    Potentials.Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
    Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )

    Energies = calc_Energies( pw, Potentials, Focc, psi )
    Etot = Energies.Total

    diff = abs(Etot-Etot_old)
    @printf("CG step %8d = %18.10f %18.10f\n", iter, Etot, diff)
    if diff < 1e-6
      @printf("CONVERGENCE ACHIEVED\n")
      break
    end

    g_old = copy(g)
    d_old = copy(d)
    Kg_old = copy(Kg)
    Etot_old = Etot
  end
  return psi, Energies, Potentials
  #
end


function calc_grad_MGC( pw::PWGrid, Potentials, Focc, psi::Array{Complex128,2} )

  # Focc are assumed to be 2s
  Ngwx    = size(psi)[1]
  Nstates = size(psi)[2]

  g = zeros(Complex128,Ngwx,Nstates)
  H = zeros(Complex128,Nstates,Nstates)
  S = zeros(Complex128,Nstates,Nstates)
  Q = zeros(Complex128,Nstates,Nstates)

  Hpsi = op_H( pw, Potentials, psi )

  for m = 1:Nstates
    for n = 1:Nstates
      H[m,n] = sum( conj(psi[:,m]).*Hpsi[:,n] )
      S[m,n] = sum( conj(psi[:,m]).*psi[:,n] )
      if m == n
        Q[m,n] = 2 - S[m,m]
      end
    end
  end

  η = 1.0  # ????

  for n = 1:Nstates
    for m = 1:Nstates
      g[:,n] = Hpsi[:,m]*Q[m,n] - psi[:,m]*Q[m,n]*η - psi[:,m]*(H[m,n] - η*S[m,n])
    end
    g[:,n] = 2*g[:,n]
  end

  return g

end
