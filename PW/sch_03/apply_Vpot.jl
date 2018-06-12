function apply_Vpot( pw::PWGrid, Vpot, psi::Array{ComplexF64,2}, ik::Int )
  #
  Ns = pw.Ns
  Ω  = pw.Ω
  Npoints = prod(Ns)
  Ncols   = size(psi)[2]

  ctmp = zeros(ComplexF64, Npoints, Ncols)
  idx = pw.gkvec.idx_gk[ik]
  for ic = 1:Ncols
    ctmp[idx,ic] = psi[:,ic]
  end

  # get values of psi in real space grid via forward transform
  ctmp = G_to_R( Ns, ctmp )

  cVpsi = R_to_G( Ns, Diagprod(Vpot, ctmp) )
  return cVpsi[idx,:]
end

# B is usually consists of more than one-column
function Diagprod( a,B )
  Ncol    = size(B)[2]
  Npoints = size(B)[1]
  out = zeros( ComplexF64, size(B) )
  for ic = 1:Ncol
    for ip = 1:Npoints
      out[ip,ic] = a[ip]*B[ip,ic]
    end
  end
  return out
end
