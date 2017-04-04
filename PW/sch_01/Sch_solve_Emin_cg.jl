function Sch_solve_Emin_cg( pw_grid::PWGrid, Vpot, psi::Array{Complex128,2};
                            NiterMax=1000 )
  #
  Npoints = size(psi)[1]
  Nstates = size(psi)[2]
  d = zeros(Complex128, Npoints, Nstates)
  g_old = zeros(Complex128, Npoints, Nstates)
  d_old = zeros(Complex128, Npoints, Nstates)
  #
  alphat = 1.e-5
  beta = 0.0
  Etot_old = 0.0
  Etot = 0.0
  #
  for iter = 1:NiterMax
    g = calc_grad( pw_grid, Vpot,  psi)
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
    psic = ortho_gram_schmidt(psi + alphat*d)
    gt = calc_grad( pw_grid, Vpot, psic )
    if real(trace((g-gt)'*d)) != 0.0
      alpha = abs(alphat*real(trace(g'*d))/real(trace((g-gt)'*d)))
    else
      alpha = 0.0
    end
    # Update wavefunction
    psi = psi[:,:] + alpha*d[:,:]

    psi = ortho_gram_schmidt(psi)
    Etot = calc_Etot( pw_grid, Vpot, psi )

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
  return psi, Etot
  #
end
