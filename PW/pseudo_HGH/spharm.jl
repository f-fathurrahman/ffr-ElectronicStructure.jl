function spharm( l, m, Gvec, cost, sint, phi )

  Ng = size(Gvec)[2]

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
        ylm[ig] = sqrt(15/16/pi)*sint[ig]^2*sin(2*phi[ig])
      end
      return ylm
    # dyz
    elseif m == -1
      for ig = 1:Ng
        ylm[ig] = sqrt(15/(4*pi))*cost[ig]*sint[ig]*sin(phi[ig])
      end
      return ylm
    # dz2
    elseif m == 0
      for ig = 1:Ng
        ylm[ig] = 0.25*sqrt(5/pi)*( 3*cost[ig]^2 - 1 )
      end
      return ylm
    # dxz
    elseif m == 1
      for ig = 1:Ng
        ylm[ig] = sqrt(15/4/pi)*cost[ig]*sint[ig]*cos(phi[ig])
      end
      return ylm
    # dx2-y2
    elseif m == 2
      for ig = 1:Ng
        ylm[ig] = 0.5*sqrt(15/4/pi)*sint[ig]^2*cos(2*phi[ig])
      end
      return ylm

    end

  elseif l == 3

    if m == -3
      for ig = 1:Ng
        ylm[ig] = 0.25*sqrt(35/2/pi)*sint[ig]^3*sin(3*phi[ig])
      end
      return ylm

    elseif m == -2
      for ig = 1:Ng
        ylm[ig] = 0.25*sqrt(105/pi)*sint[ig]^2*cost[ig]*sin(2*phi[ig])
      end
      return ylm

    elseif m == -1
      for ig = 1:Ng
        ylm[ig] = 0.25*sqrt(21/2/pi)*sint[ig]*( 5*cost[ig]^2 - 1 )*sin(phi[ig])
      end
      return ylm

    elseif m == 0
      for ig = 1:Ng
        ylm[ig] = 0.25*sqrt(7/pi)*( 5*cost[ig]^3 - 3*cost[ig] )
      end
      return ylm

    elseif m == 1
      for ig = 1:Ng
        ylm[ig] = 0.25*sqrt(21/2/pi)*sint[ig]*( 5*cost[ig]^2 - 1 )*cos(phi[ig])
      end
      return ylm

    elseif m == 2
      for ig = 1:Ng
        ylm[ig] = 0.25*sqrt(105/pi)*sint[ig]^2*cost[ig]*cos(2*phi[ig])
      end
      return ylm

    elseif m == 3
      for ig = 1:Ng
        ylm[ig] = 0.25*sqrt(35/2/pi)*sint[ig]^3*cos(3*phi[ig])
      end
      return ylm

    end

  end

end
