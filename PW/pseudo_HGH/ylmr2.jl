include("spharm.jl")

function ylmr2( lmax::Int, Gvec::Array{Float64,2} )

  Ng = size(Gvec)[2]
  lmax2 = (lmax+1)^2

  const SMALL = 1e-9

  cost = zeros(Float64, Ng)
  sint = zeros(Float64, Ng)
  phi  = zeros(Float64, Ng)

  for ig = 1:Ng
    gmod = sqrt( Gvec[1,ig]^2 + Gvec[2,ig]^2 + Gvec[3,ig]^2 )
    if gmod < SMALL
       cost[ig] = 0.0
    else
       cost[ig] = Gvec[3,ig]/gmod
    end
    #  beware the arc tan, it is defined modulo pi
    if Gvec[1,ig] > SMALL
      phi[ig] = atan( Gvec[2,ig]/Gvec[1,ig] )
    elseif Gvec[1,ig] < -SMALL
      phi[ig] = atan( Gvec[2,ig]/g[1,ig] ) + pi
    else
      phi[ig] = if Gvec[2,ig] >= 0 pi/2 else -pi/2 end
    end
    sint[ig] = sqrt( max( 0.0, 1-cost[ig]^2 ) )
  end

  ylm = zeros( Float64, Ng, lmax2 )

  ilm = 0
  for l = 0:lmax
    # m = 0
    ilm = ilm + 1
    ylm[:,ilm] = spharm( l, 0, Gvec, cost, sint, phi )
    # m > 0 and m < 0
    for m = 1:l
      ilm = ilm + 1
      ylm[:,ilm] = spharm( l, m, Gvec, cost, sint, phi )*(-1.0)^m
      #
      ilm = ilm + 1
      ylm[:,ilm] = spharm( l, -m, Gvec, cost, sint, phi )*(-1.0)^m
    end
  end

  return ylm

end
