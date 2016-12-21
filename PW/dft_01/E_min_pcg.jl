function E_min_pcg!( PW::PWGrid, Energies::EnergiesT, Potentials::PotentialsT, Focc,
                     psi0::Array{Complex128,2} ;
                     α_t=3e-5, NiterMax=500, verbose=false )
  #
  psi = copy(psi0)
  #
  Npoints = size(psi)[1]
  Nstates = size(psi)[2]
  Ω = PW.Ω
  # orthogonalize
  # TODO add option to avoid initial orthogonalization ?
  ortho_gram_schmidt!( psi )
  #
  Etot = calc_Energies!( PW, Energies, Potentials, Focc, psi )
  #
  Etot_old = Etot
  #
  g_old = zeros(Complex128,Npoints,Nstates)
  d     = zeros(Complex128,Npoints,Nstates)
  d_old = zeros(Complex128,Npoints,Nstates)
  #
  β = 0.0
  #@printf("Etot initial %18.10f\n", Etot)
  for iter = 1:NiterMax
    g = calc_grad( PW, Potentials, Focc, psi )
    #@printf("\nNew iter sum(g) = %18.10f\n", sum(g))
    #
    #norm_grad = 0.0
    #for ic = 1:Ncol
    #  norm_grad = norm_grad + norm( grad[:,ic] )
    #end
    #norm_grad = norm_grad/Ncol
    #
    if iter != 1
      β = trace( g' * Kprec(PW,g) )/trace( g_old' * Kprec(PW,g_old) )
      #β = trace( (g - g_old)' * g ) / trace( g_old' * g_old )
      #β = trace( (g-g_old)' * g ) / trace( (g-g_old)' * d )
    end
    d[:,:] = -Kprec(PW,g[:,:]) + β * d_old[:,:]
    #
    # compute gradient at trial step
    psi_t = psi + α_t * d
    ortho_gram_schmidt!( psi_t )
    g_t = calc_grad( PW, Potentials, Focc, psi_t )
    #
    # compute estimate of best step and update current trial vectorsd
    denum = trace( (g - g_t)' * d )
    if denum != 0.0
      α = abs( α_t * trace( g' * d )/denum )
    else
      α = 0.0
    end
    #
    psi[:,:] = psi[:,:] + α * d[:,:]
    #
    ortho_gram_schmidt!( psi )
    Etot = calc_Energies!( PW, Energies, Potentials, Focc, psi )
    #
    if verbose
      @printf("E min PCG: %8d %18.10f %18.10e\n", iter, Etot, abs(Etot-Etot_old))
    end
    if abs(Etot-Etot_old) < 1.e-7
      if verbose
        @printf("Emin CG converges in %8d iterations\n", iter)
      end
      break
    end
    Etot_old = Etot
    # deepcopy or reference?
    g_old[:,:] = g[:,:]
    d_old[:,:] = d[:,:]
  end
  return Etot, psi
end


function Kprec( PW::PWGrid, psi::Array{Complex128,2} )
  Ncol = size(psi)[2]
  G2 = PW.G2
  Npoints = PW.Npoints
  Kpsi = zeros( Complex128, size(psi) )
  for ic = 1:Ncol
    for ip = 1:Npoints
      Kpsi[ip,ic] = psi[ip,ic] / ( 1.0 + G2[ip] )
    end
  end
  return Kpsi
end
