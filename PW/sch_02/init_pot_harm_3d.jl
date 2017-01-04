function init_pot_harm_3d( pw::PWGrid, dr::Array{Float64,1} )
  Npoints = size(dr)[1]
  V = zeros( Npoints )
  for ip = 1:Npoints
    V[ip] = 2.0*(dr[ip]^2)
  end
  return V
end
