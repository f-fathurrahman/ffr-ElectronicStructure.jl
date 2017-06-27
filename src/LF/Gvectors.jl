struct GvectorsT
  Ns::Array{Int64,1}
  LatVecs::Array{Float64,2}
  RecVecs::Array{Float64,2}
  Ω::Float64
  G::Array{Float64,2}
  G2::Array{Float64}
end

function GvectorsT( Ns::Array{Int,1},LatVecs::Array{Float64,2} )
  RecVecs = 2*pi*inv(LatVecs')
  Ω = det(LatVecs)
  G,G2 = init_grids( Ns, LatVecs, RecVecs )
  return GvectorsT( Ns, LatVecs, RecVecs, Ω, G, G2 )
end

function mm_to_nn(mm::Int,S::Int)
  if mm > S/2
    return mm - S
  else
    return mm
  end
end

function init_grids( Ns, LatVecs, RecVecs )
  #
  Npoints = prod(Ns)
  #
  G  = Array{Float64}(3,Npoints)
  G2 = Array{Float64}(Npoints)
  ip    = 0
  for k in 0:Ns[3]-1
  for j in 0:Ns[2]-1
  for i in 0:Ns[1]-1
    gi = mm_to_nn( i, Ns[1] )
    gj = mm_to_nn( j, Ns[2] )
    gk = mm_to_nn( k, Ns[3] )
    ip = ip + 1
    G[1,ip] = RecVecs[1,1]*gi + RecVecs[1,2]*gj + RecVecs[1,3]*gk
    G[2,ip] = RecVecs[2,1]*gi + RecVecs[2,2]*gj + RecVecs[2,3]*gk
    G[3,ip] = RecVecs[3,1]*gi + RecVecs[3,2]*gj + RecVecs[3,3]*gk
    G2[ip] = G[1,ip]^2 + G[2,ip]^2 + G[3,ip]^2
  end
  end
  end
  return G,G2
end
