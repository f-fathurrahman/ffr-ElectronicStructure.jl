# Gvec is of shape (3,Ng)
# the argument lmax is different from the corresponding
# Fortran subroutine
function ylmr2( lmax::Int, Gvec::Array{Float64,2} )

  Ng = size(Gvec)[2]
  lmax2 = (lmax+1)^2

  ylm = zeros( Float64, Ng, lmax2 )

  gg = zeros(Float64,Ng)

  for ig = 1:Ng
    gg[ig] = Gvec[1,ig]^2 + Gvec[2,ig]^2 + Gvec[3,ig]^2
  end

  const SMALL = 1.e-9
  const FPI = 4.0*pi

  if (Ng < 1) || (lmax2 < 1)
    return ylm
  end

  if lmax == 0
    ylm[:,1] = sqrt(1.0/FPI)
    return ylm
  end

  cost = zeros(Float64,Ng)
  sent = zeros(Float64,Ng)
  phi  = zeros(Float64,Ng)
  Q    = zeros(Float64,Ng,lmax+1,lmax+1)  # originally 0:lmax

  for ig = 1:Ng

    gmod = sqrt( gg[ig] )

    if gmod < SMALL
      cost[ig] = 0.d0
    else
      cost[ig] = Gvec[3,ig]/gmod
    end

    # beware the arc tan, it is defined modulo pi
    # NOTE: verify this in Julia
    if Gvec[1,ig] > SMALL
      phi[ig] = atan( Gvec[2,ig]/Gvec[1,ig] )
      #
    elseif Gvec[1,ig] < -SMALL
      phi[ig] = atan( Gvec[2,ig]/Gvec[1,ig] ) + pi
      #
    else
      phi(ig) = sign( pi/2., Gvec[2,ig] )  # XXX CHECK THIS
    end

    sent[ig] = sqrt( max( 0.0, 1.0 - cost[ig]^2) )

  end # ig

  #
  # Q(:,l,m) are defined as sqrt ((l-m)!/(l+m)!) * P(:,l,m) where
  # P(:,l,m) are the Legendre Polynomials (0 <= m <= l)
  #
  lm = 0
  for l = 0:lmax
    @printf("l = %d\n", l)
    c = sqrt( (2*l+1) / FPI )

    if l == 0
      for ig = 1:Ng
        Q[ig,1,1] = 1.0  # l=0 -> 1
      end

    elseif l == 1
      for ig = 1:Ng
        Q[ig,2,1] =  cost[ig]
        Q[ig,2,2] = -sent[ig]/sqrt(2.0)
      end

    else
      # recursion on l for Q(:,l,m)
      for m = 0:(l-2)
        println("Enter this ...")
        for ig = 1:Ng
          Q[ig,l+1,m+1] = cost[ig]*(2*l-1)/sqrt(l*l-m*m) * Q[ig,l,m+1]
                         - sqrt((l-1)*(l-1)-m*m)/sqrt(l*l-m*m) * Q[ig,l-1,m+1]
        end
        println("sum(Q) = ", sum(Q[:,l+1,m+1]))
      end
      for ig = 1:Ng
        Q[ig,l+1,l] = cost[ig] * sqrt(2*l-1) * Q[ig,l,l]
      end
      for ig = 1:Ng
        Q[ig,l+1,l+1] = - sqrt(2*l-1)/sqrt(2*l)*sent[ig]*Q[ig,l,l]
      end

    end

    #
    # Y_lm, m = 0
    #
    lm = lm + 1
    @printf("m = 0: lm = %d\n", lm)
    for ig = 1:Ng
      ylm[ig, lm] = c * Q[ig,l+1,1]
    end

    for m = 1:l
      #
      # Y_lm, m > 0
      #
      lm = lm + 1
      @printf("m > 0: lm = %d\n", lm)
      for ig = 1:Ng
        ylm[ig, lm] = c * sqrt(2.0) * Q[ig,l+1,m+1] * cos( m*phi[ig] )
      end
      #
      # Y_lm, m < 0
      #
      lm = lm + 1
      @printf("m < 0: lm = %d\n", lm)
      for ig = 1:Ng
        ylm[ig,lm] = c * sqrt(2.0) * Q[ig,l+1,m+1] * sin(m*phi[ig])
      end
    end

  end # lm

  return ylm

end
