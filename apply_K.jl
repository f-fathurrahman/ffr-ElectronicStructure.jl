function apply_K( PW::PWGrid, psi::Array{Complex128,2} )
  out = zeros(Complex128,size(psi))
  Ncol = size(psi,2)
  Ω = PW.Ω
  G2 = PW.G2
  Npoints = PW.Npoints
  for is = 1:Ncol
    for ip = 1:Npoints
      out[ip,is] = psi[ip,is]*G2[ip]
    end
  end
  return scale(0.5,out) # two minus signs
end

function apply_K( PW::PWGrid, psi::Array{Complex128,1} )
  G2 = PW.G2
  Npoints = PW.Npoints
  out = zeros( size(psi) )
  for ip = 1:Npoints
    out[ip] = psi[ip]*G2[ip]
  end
  return scale(-0.5,out)
end

