function schsolve_Emin_cg( LF::LF3dGrid, Vpot, Ncol::Int64;
                           v0 = nothing,
                           α_t=3e-5, Niter=1000,
                           verbose=false )
  #
  Npoints = LF.Nx * LF.Ny * LF.Nz
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
  #
  # Setup trial solution
  #
  if v0 == nothing
    # random trial solution
    srand(1234)
    v = rand( Npoints, Ncol )
    v = orthonormalize( LF, v )
  else
    v = copy(v0)
    # v0 is assumed to be ortonormalized
  end
  #
  Energies = calc_Energies( LF, Vpot, v )
  #
  Etot_old = Energies.Total
  #
  g     = zeros(Float64,Npoints,Ncol)
  g_old = zeros(Float64,Npoints,Ncol)
  d     = zeros(Float64,Npoints,Ncol)
  d_old = zeros(Float64,Npoints,Ncol)
  #
  β = 0
  α = 0
  for iter = 1:Niter
    g = calc_grad( LF, Vpot, v )
    #
    if iter != 1
      β = trace( g' * g)/trace( g_old'*g_old )
      #β = trace( (g-g_old)'*g ) / trace( g_old'*g_old )
      #β = trace( (g-g_old)'*g )/ trace( (g-g_old)'*d )
    end
    d = -g + β*d_old
    #
    # compute gradient at trial step
    #
    v2 = orthonormalize(LF, v + α_t*d)
    g_t = calc_grad( LF, Vpot, v2 )
    #
    # compute estimate of best step and update current trial vectors
    #
    denum = trace( (g - g_t)' * d )
    if denum != 0.0
      α = abs( α_t * trace( g'*d )/denum )
    else
      α = 0.0
    end
    v = v + α*d
    #
    v = orthonormalize( LF, v )
    Energies = calc_Energies( LF, Vpot, v )
    Etot = Energies.Total
    #
    if verbose
      @printf("%8d %16.9e %16.9e\n", iter, Etot, abs(Etot-Etot_old))
    end
    if abs(Etot-Etot_old) < 1.e-7
      if verbose
        @printf("Emin CG converges in %8d iterations\n", iter)
      end
      break
    end
    Etot_old = Etot
    #
    g_old = g[:,:]
    d_old = d[:,:]
  end
  return Energies, v
end
