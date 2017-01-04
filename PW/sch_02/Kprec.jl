function Kprec( pw::PWGrid, psi::Array{Complex128,2} )

  Ngwx  = size(psi)[1]
  Ncols = size(psi)[2]
  G2    = pw.gvec.G2[pw.gvecw.idx_gw2r]
  Kpsi  = zeros( Complex128, Ngwx, Ncols )

  for ic = 1:Ncols
    for ig = 1:Ngwx
      Kpsi[ig,ic] = psi[ig,ic] / ( 1.0 + G2[ig] )
    end
  end
  return Kpsi
end

function Kprec(psi)
  return psi
end
