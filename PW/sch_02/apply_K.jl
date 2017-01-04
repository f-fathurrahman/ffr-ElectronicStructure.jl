# Apply kinetic operator to wave function in reciprocal space

function apply_K( pw::PWGrid, psi::Array{Complex128,2} )
  #
  out = zeros(Complex128,size(psi))
  Ncol = size(psi)[2]

  Ngwx = pw.gvecw.Ngwx
  G2   = pw.gvec.G2[pw.gvecw.idx_gw2r]

  for is = 1:Ncol
    for ig = 1:Ngwx
      out[ig,is] = psi[ig,is]*G2[ig]
    end
  end
  return 0.5*out # two minus signs
end



function apply_K( pw::PWGrid, psi::Array{Complex128,1} )

  G2   = pw.G2[pw.gvecw.idx_gw2r]
  Ngwx = pw.gvecw.Npoints
  out  = zeros( Complex128, size(psi) )

  for ig = 1:Ngwx
    out[ig] = psi[ig]*G2[ig]
  end

  return 0.5*out
end
