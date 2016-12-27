
type GVectors
  Ng::Int
  G::Array{Float64,2}
  G2::Array{Float64}
  idx_g2r::Array{Int}
end

type PWGrid
  ecutwfc::Float64
  ecutrho::Float64
  Ns::Array{Int64,1}
  LatVecs::Array{Float64,2}
  RecVecs::Array{Float64,2}
  Ω::Float64
  R::Array{Float64,2}
  gvectors::GVectors
end

function PWGrid( ecutwfc::Float64, LatVecs::Array{Float64,2} )
  ecutrho = 4.0*ecutwfc
  #
  RecVecs = 2*pi*inv(LatVecs')
  Ω = det(LatVecs)
  #
  LatVecsLen = zeros(3)
  LatVecsLen[1] = norm(LatVecs[1,:])
  LatVecsLen[2] = norm(LatVecs[2,:])
  LatVecsLen[3] = norm(LatVecs[3,:])

  Ns = zeros(Int64,3)
  Ns[1] = 2*round( Int, sqrt(ecutrho/2)*LatVecsLen[1]/pi ) + 2
  Ns[2] = 2*round( Int, sqrt(ecutrho/2)*LatVecsLen[2]/pi ) + 2
  Ns[3] = 2*round( Int, sqrt(ecutrho/2)*LatVecsLen[3]/pi ) + 2

  @printf("Sampling points: (%8d,%8d,%8d)\n", Ns[1], Ns[2], Ns[3])

  Npoints = prod(Ns)
  R = init_grids_R( Ns, LatVecs )

  gvectors = init_gvectors( ecutrho, Ns, RecVecs )

  return PWGrid( ecutwfc, ecutrho, Ns, LatVecs, RecVecs, Ω, R, gvectors )
end

function mm_to_nn(mm::Int,S::Int)
  if mm > S/2
    return mm - S
  else
    return mm
  end
end


function init_gvectors( ecutrho, Ns, RecVecs )

  ig = 0
  G_temp = zeros(3)

  #
  for k in 0:Ns[3]-1
  for j in 0:Ns[2]-1
  for i in 0:Ns[1]-1
    gi = mm_to_nn( i, Ns[1] )
    gj = mm_to_nn( j, Ns[2] )
    gk = mm_to_nn( k, Ns[3] )
    G_temp[1] = RecVecs[1,1]*gi + RecVecs[2,1]*gj + RecVecs[3,1]*gk
    G_temp[2] = RecVecs[1,2]*gi + RecVecs[2,2]*gj + RecVecs[3,2]*gk
    G_temp[3] = RecVecs[1,3]*gi + RecVecs[2,3]*gj + RecVecs[3,3]*gk
    Glen = norm(G_temp)^2
    if Glen <= ecutrho
      ig = ig + 1
    end
  end
  end
  end
  Ngvec = ig

  #
  G  = Array(Float64,3,Ngvec)
  G2 = Array(Float64,Ngvec)
  #
  ig = 0
  ir = 0
  idx_g2r = zeros(Int,Ngvec)
  for k in 0:Ns[3]-1
  for j in 0:Ns[2]-1
  for i in 0:Ns[1]-1
    ir = ir + 1
    gi = mm_to_nn( i, Ns[1] )
    gj = mm_to_nn( j, Ns[2] )
    gk = mm_to_nn( k, Ns[3] )
    G_temp[1] = RecVecs[1,1]*gi + RecVecs[2,1]*gj + RecVecs[3,1]*gk
    G_temp[2] = RecVecs[1,2]*gi + RecVecs[2,2]*gj + RecVecs[3,2]*gk
    G_temp[3] = RecVecs[1,3]*gi + RecVecs[2,3]*gj + RecVecs[3,3]*gk
    Glen = norm(G_temp)^2
    if Glen <= ecutrho
      ig = ig + 1
      idx_g2r[ig] = ir
      G[:,ig] = G_temp[:]
      G2[ig]  = Glen
    end
  end
  end
  end

  return GVectors( Ngvec, G, G2, idx_g2r )
end

function init_grids_R( Ns, LatVecs )
  #
  Npoints = prod(Ns)
  #
  R = Array(Float64,3,Npoints)
  ip = 0
  for k in 0:Ns[3]-1
  for j in 0:Ns[2]-1
  for i in 0:Ns[1]-1
    ip = ip + 1
    R[1,ip] = LatVecs[1,1]*i/Ns[1] + LatVecs[2,1]*j/Ns[2] + LatVecs[3,1]*k/Ns[3]
    R[2,ip] = LatVecs[1,2]*i/Ns[1] + LatVecs[2,2]*j/Ns[2] + LatVecs[3,2]*k/Ns[3]
    R[3,ip] = LatVecs[1,3]*i/Ns[1] + LatVecs[2,3]*j/Ns[2] + LatVecs[3,3]*k/Ns[3]
  end
  end
  end
  #
  return R
end
