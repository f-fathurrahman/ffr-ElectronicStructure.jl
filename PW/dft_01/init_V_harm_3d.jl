# Initialize 3d harmonic potential
# $$
# V(r) = \frac{1}{2} r^2
# $$

function init_V_harm_3d( PW::PWGrid, dr::Array{Float64,1} )
  Npoints = size(dr)[1]
  Ns = PW.Ns
  Ω = PW.Ω
  V = zeros( Npoints )
  for ip = 1:Npoints
    V[ip] = 2.0*(dr[ip]^2)
  end
  return V
end

