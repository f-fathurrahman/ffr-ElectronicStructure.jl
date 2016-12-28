function apply_K( pw_grid::PWGrid, psi::Array{Complex128,2} )
  out = zeros(Complex128,size(psi))
  Ncol = size(psi,2)
  Ω = pw_grid.Ω
  G2 = pw_grid.G2
  Npoints = pw_grid.Npoints
  for is = 1:Ncol
    for ip = 1:Npoints
      out[ip,is] = psi[ip,is]*G2[ip]
    end
  end
  return 0.5*out # two minus signs
end

function apply_K( pw_grid::PWGrid, psi::Array{Complex128,1} )
  G2 = pw_grid.G2
  Npoints = pw_grid.Npoints
  out = zeros( size(psi) )
  for ip = 1:Npoints
    out[ip] = psi[ip]*G2[ip]
  end
  return 0.5*out
end
