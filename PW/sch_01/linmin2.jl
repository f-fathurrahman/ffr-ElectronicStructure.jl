function linmin2( pw, Vpot, α_t, psi, d, g )

  psic = ortho_gram_schmidt(psi + α_t*d)
  gt = calc_grad( pw, Vpot, psic )
  denum = real(sum(conj(g-gt).*d))
  if denum != 0.0
    α = abs(α_t*real(sum(conj(g).*d))/denum)
  else
    α = 0.0
  end

  return α

end
