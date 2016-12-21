function apply_V_loc( PW::PWGrid, V_loc::Array{Float64,1},
                      psi::Array{Complex128,2} )
  #
  Ns = PW.Ns
  Ω = PW.Ω
  Npoints = prod(Ns)
  # get values of psi in real space grid via forward transform
  ctmp = G_to_R( Ns, psi )
  return R_to_G( Ns, Diagprod(V_loc, ctmp) )
end

# B is usually consists of more than one-column
function Diagprod( a,B )
  Ncol    = size(B)[2]
  Npoints = size(B)[1]
  out = zeros( Complex128, size(B) )  # in general the output will be complex numbers
  for ic = 1:Ncol
    for ip = 1:Npoints
      out[ip,ic] = a[ip]*B[ip,ic]
    end
  end
  return out
end
