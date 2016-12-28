type PWGrid
  Ns::Array{Int64,1}
  LatVecs::Array{Float64,2}
  RecVecs::Array{Float64,2}
  Npoints::Int
  Ω::Float64
  R::Array{Float64,2}
  G::Array{Float64,2}
  G2::Array{Float64}
end

function PWGrid( Ns::Array{Int,1},LatVecs::Array{Float64,2} )
  Npoints = prod(Ns)
  RecVecs = 2*pi*inv(LatVecs')
  Ω = det(LatVecs)
  R,G,G2 = init_grids( Ns, LatVecs, RecVecs )
  return PWGrid( Ns, LatVecs, RecVecs, Npoints, Ω, R, G, G2 )
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
  R = Array(Float64,3,Npoints)
  ip = 0
  for k in 0:Ns[3]-1
  for j in 0:Ns[2]-1
  for i in 0:Ns[1]-1
    ip = ip + 1
    R[1,ip] = LatVecs[1,1]*i/Ns[1] + LatVecs[1,2]*j/Ns[2] + LatVecs[1,3]*k/Ns[3]
    R[2,ip] = LatVecs[2,1]*i/Ns[1] + LatVecs[2,2]*j/Ns[2] + LatVecs[2,3]*k/Ns[3]
    R[3,ip] = LatVecs[3,1]*i/Ns[1] + LatVecs[3,2]*j/Ns[2] + LatVecs[3,3]*k/Ns[3]
  end
  end
  end
  #
  G  = Array(Float64,3,Npoints)
  G2 = Array(Float64,Npoints)
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
  return R,G,G2
end
