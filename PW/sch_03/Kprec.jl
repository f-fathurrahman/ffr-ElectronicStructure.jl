function Kprec( pw::PWGrid, psi::Array{Complex128,2}, ik::Int )

  Ngwx  = size(psi)[1]
  Ncols = size(psi)[2]

  Ngwk = pw.gkvec.Ngw[ik]
  idx_gkk = pw.gkvec.idx_gk[ik]
  G = pw.gvec.G[:,idx_gkk]

  Kpsi  = zeros( Complex128, Ngwk, Ncols )

  for ic = 1:Ncols
    for ig = 1:Ngwk
      Gk = G[:,ig] + pw.gkvec.kpts[:,ik]
      Kpsi[ig,ic] = psi[ig,ic] / ( 1.0 + Gk[1]^2 + Gk[2]^2 + Gk[3]^2 )
    end
  end
  return Kpsi
end

function Kprec(psi)
  return psi
end
