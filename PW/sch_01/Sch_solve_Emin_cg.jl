function Sch_solve_Emin_cg( pw::PWGrid, Vpot, psi::Array{Complex128,2};
                            NiterMax=1000 )
  #
  Npoints = size(psi)[1]
  Nstates = size(psi)[2]
  d = zeros(Complex128, Npoints, Nstates)
  g_old  = zeros(Complex128, Npoints, Nstates)
  d_old  = zeros(Complex128, Npoints, Nstates)
  Kg     = zeros(Complex128, Npoints, Nstates)
  Kg_old = zeros(Complex128, Npoints, Nstates)
  #
  α_t = 1.e-5
  β = 0.0
  Etot_old = 0.0
  Etot = 0.0
  #
  for iter = 1:NiterMax
    g = calc_grad( pw, Vpot,  psi)
    nrm = 0.0
    for is = 1:Nstates
      nrm = nrm + real( dot( g[:,is], g[:,is] ) )
    end
    Kg = Kprec(pw,g)
    if iter != 1
      #β = real(sum(conj(g).*Kg))/real(sum(conj(g_old).*Kg_old))
      β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g_old).*Kg_old))
      #β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g-g_old).*d_old))
      #β = real(sum(conj(g).*Kg))/real(sum((g-g_old).*conj(d_old)))
    end
    d = -Kg + β * d_old
    psic = ortho_gram_schmidt(psi + α_t*d)
    gt = calc_grad( pw, Vpot, psic )
    denum = real(sum(conj(g-gt).*d))
    if denum != 0.0
      α = abs(α_t*real(sum(conj(g).*d))/denum)
    else
      α = 0.0
    end
    # Update wavefunction
    psi = psi[:,:] + α*d[:,:]

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
    Kg_old = copy(Kg)
    Etot_old = Etot
  end
  return psi, Etot
  #
end
