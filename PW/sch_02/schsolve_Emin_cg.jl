function schsolve_Emin_cg( pw::PWGrid, Vpot, psi::Array{ComplexF64,2};
                           NiterMax=1000 )
  #
  Npoints = size(psi)[1]
  Nstates = size(psi)[2]
  d = zeros(ComplexF64, Npoints, Nstates)
  g_old = zeros(ComplexF64, Npoints, Nstates)
  d_old = zeros(ComplexF64, Npoints, Nstates)
  #
  alphat = 1.e-5
  beta = 0.0
  Etot_old = 0.0
  Etot = 0.0
  #
  for iter = 1:NiterMax
    g = calc_grad( pw, Vpot,  psi)
    #println("sum(g) = ", sum(g))
    #println("sum(Kg) = ", sum(Kprec(pw,g)))
    nrm = 0.0
    for is = 1:Nstates
      nrm = nrm + real( dot( g[:,is], g[:,is] ) )
    end
    if iter != 1
      #beta = real(tr(g'*Kprec(pw,g)))/real(tr(g_old'*Kprec(pw,g_old)))
      beta = real(tr((g-g_old)'*Kprec(pw,g)))/real(tr(g_old'*Kprec(pw,g_old)))
      #beta = real(tr((g-g_old)'*Kprec(pw,g)))/real(tr((g-g_old)'*d))
    end
    d = -Kprec(pw, g) + beta * d_old
    #println("sum(d) = ", sum(d))

    psic = ortho_gram_schmidt(psi + alphat*d)
    #println("sum(tv) = ", sum(psi+alphat*d))
    #println("sum(psic) = ", sum(psic))

    gt = calc_grad( pw, Vpot, psic )

    denum = real(tr((g-gt)'*d))
    #println("denum = ", denum)
    if denum != 0.0
      alpha = abs(alphat*real(tr(g'*d))/denum )
    else
      alpha = 0.0
    end
    #@printf("\nalpha, beta = %f %f\n", alpha, beta)
    # Update wavefunction
    psi = psi[:,:] + alpha*d[:,:]

    psi = ortho_gram_schmidt(psi)
    Etot = calc_Etot( pw, Vpot, psi )

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
