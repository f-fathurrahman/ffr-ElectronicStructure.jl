function spharm( l, m, Gvec )

  Ng = size(Gvec)[2]

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
    #
    #  beware the arc tan, it is defined modulo pi
    #
    if Gvec[1,ig] > SMALL
      phi[ig] = atan( Gvec[2,ig]/Gvec[1,ig] )
    elseif Gvec[1,ig] < -SMALL
      phi[ig] = atan( Gvec[2,ig]/g[1,ig] ) + pi
    else
      phi[ig] = if Gvec[2,ig] >= 0 pi/2 else -pi/2 end
    end
    sint[ig] = sqrt( max( 0.0, 1-cost[ig]^2 ) )
  end

  ylm = zeros( Float64, Ng )

  if l == 0
    for ig = 1:Ng
      ylm[ig] = 0.5*sqrt(1.0/pi)
    end
    return ylm

  elseif l == 1
    # py
    if m == -1
      for ig = 1:Ng
        ylm[ig] = 0.5*sqrt(3/pi)*sint[ig]*sin(phi[ig])
      end
      return ylm
    # pz
    elseif m == 0
      for ig = 1:Ng
        ylm[ig] = 0.5*sqrt(3/pi)*cost[ig]
      end
      return ylm
    # px
    elseif m == 1
      for ig = 1:Ng
        ylm[ig] = 0.5*sqrt(3/pi)*sint[ig]*cos(phi[ig])
      end
      return ylm
    end

  elseif l == 2
    # dxy
    if m == -2
      for ig = 1:Ng
        #ylm[ig] = sqrt(oppi*15/16)*sin(phi[ig])^2*sin(2*phi[ig])
        ylm[ig] = sqrt(15/16/pi)*sint[ig]^2*sin(2*phi[ig])
      end
      return ylm
    # dyz
    elseif m == -1
      for ig = 1:Ng
        #ylm[ig] = sqrt(oppi*15/4)*cost[ig]*sint[ig]*sin(phi[ig])
        ylm[ig] = sqrt(15/(4*pi))*cost[ig]*sint[ig]*sin(phi[ig])
      end
      return ylm
    # dz2
    elseif m == 0
      for ig = 1:Ng
        #ylm[ig] = sqrt(oppi*5/8)*( 2*cost[ig]^2 - sin(phi[ig])^2 )
        ylm[ig] = 0.25*sqrt(5/pi)*( 3*cost[ig]^2 - 1 )
      end
      return ylm
    # dxz
    elseif m == 1
      for ig = 1:Ng
        #ylm[ig] = sqrt(oppi*15/4)*cost[ig]*sint[ig]*cos(phi[ig])
        ylm[ig] = sqrt(15/4/pi)*cost[ig]*sint[ig]*cos(phi[ig])
      end
      return ylm
    # dx2-y2
    elseif m == 2
      for ig = 1:Ng
        #ylm[ig] = sqrt(oppi*15/16)*sin(phi[ig])^2*cos(2*phi[ig])
        ylm[ig] = 0.5*sqrt(15/4/pi)*sint[ig]^2*cos(2*phi[ig])
      end
      return ylm

    end

  end

end
