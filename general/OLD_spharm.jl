function spharm( l, m, Gvec )

  Ng = size(Gvec)[2]

  r = zeros( Float64, Ng )
  for ig = 1:Ng
    r[ig] = sqrt( Gvec[1,ig]^2 + Gvec[2,ig]^2 + Gvec[3,ig]^2 )
  end
  # This will fail for very small r[ig]

  ylm = zeros( Float64, Ng )

  if l == 0
    for ig = 1:Ng
      ylm[ig] = 0.5*sqrt(1.0/pi)
    end
    return ylm

  elseif l == 1
    #
    if m == -1
      for ig = 1:Ng
        ylm[ig] = 0.5*sqrt(3/pi)*Gvec[2,ig]/r[ig]
      end
      return ylm
    #
    elseif m == 0
      for ig = 1:Ng
        ylm[ig] = 0.5*sqrt(3/pi)*Gvec[3,ig]/r[ig]
      end
      return ylm
    #
    elseif m == 1
      for ig = 1:Ng
        ylm[ig] = 0.5*sqrt(3/pi)*Gvec[1,ig]/r[ig]
      end
      return ylm
    end

  elseif l == 2
    #
    if m == -2
      for ig = 1:Ng
        ylm[ig] = 0.5*sqrt(15/pi)*Gvec[1,ig]*Gvec[2,ig]/r[ig]^2
      end
      return ylm
    #
    elseif m == -1
      for ig = 1:Ng
        ylm[ig] = 0.5*sqrt(15/pi)*Gvec[2,ig]*Gvec[3,ig]/r[ig]^2
      end
      return ylm
    #
    elseif m == 0
      for ig = 1:Ng
        ylm[ig] = 0.25*sqrt(5/pi)*( -Gvec[1,ig]^2 - Gvec[2,ig]^2 + 2*Gvec[3,ig]^2 )/r[ig]^2
      end
      return ylm
    #
    elseif m == 1
      for ig = 1:Ng
        ylm[ig] = 0.5*sqrt(15/pi)*Gvec[3,ig]*Gvec[1,ig]/r[ig]^2
      end
      return ylm
    #
    elseif m == 2
      for ig = 1:Ng
        ylm[ig] = 0.25*sqrt(15/pi)*(Gvec[1,ig]^2 - Gvec[2,ig]^2)/r[ig]^2
      end
      return ylm

    end

  end

end
