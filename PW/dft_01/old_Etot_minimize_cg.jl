function Etot_minimize_cg(  pw_grid::PWGrid, psi::Array{Complex128,2},
                   Vpot::Array{Complex128,1}, Niter )
  #
  d = zeros(Complex128,size(psi))
  g_old = zeros(Complex128,size(psi))
  d_old = zeros(Complex128,size(psi))
  Nstates = size(psi,2)
  #
  alphat = 1.e-5
  beta = 0.0
  Etot_old = 0.0
  Etot = 0.0
  #
  for iter = 1:Niter
    g = gradE( pw_grid, psi, Vpot )
    nrm = 0.0
    for is = 1:Nstates
      nrm = nrm + real( dot( g[:,is], g[:,is] ) )
    end
    if iter != 1
      beta = real(trace(g'*Kprec(pw_grid,g)))/real(trace(g_old'*Kprec(pw_grid,g_old)))
      #beta = real(trace((g-g_old)'*Kprec(pw_grid,g)))/real(trace(g_old'*Kprec(pw_grid,g_old)))
      #beta = real(trace((g-g_old)'*Kprec(pw_grid,g)))/real(trace((g-g_old)'*d))
      #@printf("\nbeta = %f\n", beta)
    end
    d = -Kprec(pw_grid, g) + beta * d_old
    gt = gradE( pw_grid, psi + alphat*d, Vpot )
    if real(trace((g-gt)'*d)) != 0.0
      alpha = abs(alphat*real(trace(g'*d))/real(trace((g-gt)'*d)))
      #@printf("alpha = %f\n",alpha)
    else
      #@printf("alpha is set to zero\n")
      alpha = 0.0
    end
    psi = psi + alpha*d
    Etot = get_Etot( pw_grid, psi, Vpot )
    diff = abs(Etot-Etot_old)
    @printf("E step %8d = %20.16f %20.16f %20.16f\n", iter, Etot, diff, nrm/Nstates)
    if diff < 1e-6
      @printf("CONVERGENCE ACHIEVED\n")
      break
    end
    g_old = g
    d_old = d
    Etot_old = Etot
  end
  return psi, Etot
  #
end


function Kprec( pw_grid::PWGrid, psi::Array{Complex128,2} )
  Ncol = size(psi)[2]
  G2 = pw_grid.G2
  Npoints = pw_grid.Npoints
  Kpsi = zeros( Complex128, size(psi) )
  for ic = 1:Ncol
    for ip = 1:Npoints
      Kpsi[ip,ic] = psi[ip,ic] / ( 1.0 + G2[ip] )
    end
  end
  return Kpsi
end

function Kprec(psi)
  return psi
end
