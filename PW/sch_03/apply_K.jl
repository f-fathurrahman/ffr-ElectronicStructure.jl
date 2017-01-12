# Apply kinetic operator to wave function in reciprocal space

function apply_K( pw::PWGrid, psi::Array{Complex128,2}, ik::Int )
  #
  out = zeros(Complex128,size(psi))
  Ncol = size(psi)[2]

  Ngwk = pw.gkvec.Ngw[ik]
  idx_gkk = pw.gkvec.idx_gk[ik]
  G = pw.gvec.G[:,idx_gkk]

  for is = 1:Ncol
    for ig = 1:Ngwk
      Gk = G[:,ig] + pw.gkvec.kpts[:,ik]
      out[ig,is] = psi[ig,is]*( Gk[1]^2 + Gk[2]^2 + Gk[3]^2 )
    end
  end
  return 0.5*out # two minus signs
end
