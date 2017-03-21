"""
Two electron integrals
"""

function coulomb(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,
         aI::Int64,aJ::Int64,aK::Int64,
         bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
         bI::Int64,bJ::Int64,bK::Int64,
         cexpn::Float64,cx::Float64,cy::Float64,cz::Float64,
         cI::Int64,cJ::Int64,cK::Int64,
         dexpn::Float64,dx::Float64,dy::Float64,dz::Float64,
         dI::Int64,dJ::Int64,dK::Int64)
  # This is the slow method of computing integrals from Huzinaga et al.
  # Use the HRR/VRR scheme from Head-Gordon & Pople instead

  rab2 = dist2(ax-bx,ay-by,az-bz)
  rcd2 = dist2(cx-dx,cy-dy,cz-dz)

  px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
  qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,dx,dy,dz)
  rpq2 = dist2(px-qx,py-qy,pz-qz)

  g1 = aexpn+bexpn
  g2 = cexpn+dexpn
  delta = 0.25*(1/g1+1/g2)

  Bx = Barray(aI,bI,cI,dI,px,ax,bx,qx,cx,dx,g1,g2,delta)
  By = Barray(aJ,bJ,cJ,dJ,py,ay,by,qy,cy,dy,g1,g2,delta)
  Bz = Barray(aK,bK,cK,dK,pz,az,bz,qz,cz,dz,g1,g2,delta)

  s = 0
  #println("$(aI+bI+cI+dI),$(aJ+bJ+cJ+dJ),$(aK+bK+cK+dK)")
  for I in 0:(aI+bI+cI+dI)
    for J in 0:(aJ+bJ+cJ+dJ)
      for K in 0:(aK+bK+cK+dK)
        #println("coul: $I,$J,$K,$(Bx[I+1]),$(By[J+1]),$(Bz[K+1])")
        s += Bx[I+1]*By[J+1]*Bz[K+1]*Fgamma(I+J+K,0.25*rpq2/delta)
      end
    end
  end
  return 2pi^(2.5)/(g1*g2*sqrt(g1+g2))*exp(-aexpn*bexpn*rab2/g1)*exp(-cexpn*dexpn*rcd2/g2)*s
end

function coulomb(a::PGBF,b::PGBF,c::PGBF,d::PGBF)
  return a.norm*b.norm*c.norm*d.norm*coulomb(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
    b.expn,b.x,b.y,b.z,b.I,b.J,b.K,
    c.expn,c.x,c.y,c.z,c.I,c.J,c.K,
    d.expn,d.x,d.y,d.z,d.I,d.J,d.K)
end

fB(i::Int64,l1::Int64,l2::Int64,p::Float64,a::Float64,b::Float64,r::Int64,g::Float64) = binomial_prefactor(i,l1,l2,p-a,p-b)*B0(i,r,g)
B0(i::Int64,r::Int64,g::Float64) = fact_ratio2(i,r)*(4g)^(r-i)
fact_ratio2(a::Int64,b::Int64) = factorial(a,b)/factorial(a-2b)

function Bterm(i1::Int64,i2::Int64,r1::Int64,r2::Int64,u::Int64,
         l1::Int64,l2::Int64,l3::Int64,l4::Int64,
         Px::Float64,Ax::Float64,Bx::Float64,Qx::Float64,Cx::Float64,Dx::Float64,
         gamma1::Float64,gamma2::Float64,delta::Float64)
  # THO eq. 2.22
  #print("Bterm($i1,$i2,$r1,$r2,$u,$l1,$l2,$l3,$l4,$Px,$Ax,$Bx,$Qx,$Cx,$Dx,$gamma1,$gamma2,$delta)=")
  val = (-1)^(i2+u)*fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1)*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2)*(
      fact_ratio2(i1+i2-2*(r1+r2),u)*(Qx-Px)^(i1+i2-2*(r1+r2)-2*u)/delta^(i1+i2-2*(r1+r2)-u))
  #println("$val")
  return val
end

function Barray(l1::Int64,l2::Int64,l3::Int64,l4::Int64,p::Float64,a::Float64,b::Float64,
        q::Float64,c::Float64,d::Float64,g1::Float64,g2::Float64,delta::Float64)
  Imax = l1+l2+l3+l4+1
  B = zeros(Float64,Imax)
  for i1 in 0:(l1+l2)
    for i2 in 0:(l3+l4)
      for r1 in 0:div(i1,2)
        for r2 in 0:div(i2,2)
          for u in 0:(div(i1+i2,2)-r1-r2)
            I = i1+i2-2*(r1+r2)-u
            B[I+1] += Bterm(i1,i2,r1,r2,u,l1,l2,l3,l4,p,a,b,q,c,d,g1,g2,delta)
          end
        end
      end
    end
  end
  return B
end

coulomb(a::CGBF,b::CGBF,c::CGBF,d::CGBF) = contract(coulomb,a,b,c,d)


function test_two_terms()
  @printf("Calling test_two_terms\n")
  @assert fB(0,0,0,0.0,0.0,0.0,0,2.0) == 1
  @assert fB(0,0,0,1.0,1.0,1.0,0,2.0) == 1
  @assert fB(0,0,0,0.0,0.0,0.0,0,2.0 ) == 1
  @assert fB(1,0,1,0.0,0.0,0.0,0,2.0 ) == 0.125
  @assert B0(0,0,2.0) == 1
  @assert fact_ratio2(0,0) == 1
  @assert Bterm(0,0,0,0,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,2.0,2.0,0.25)==1
  @assert Bterm(0,1,0,0,0,0,0,0,1,0.0,0.0,0.0,0.0,0.0,0.0,2.0,2.0,0.25)==0

  @printf("Pass test_two_terms\n")
end


function test_coul1()
  @printf("Calling test_coul1\n")
  s = pgbf(1.0)
  px = pgbf(1.0,0,0,0,1,0,0)
  @assert coulomb(s,s,s,px)==0 # 0
  @assert isapprox(coulomb(s,s,px,px), 0.9403159725793305 )
  @printf("Pass test_coul1\n")
end

function coulomb_hgp(a::PGBF,b::PGBF,c::PGBF,d::PGBF)
  return a.norm*b.norm*c.norm*d.norm*hrr(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
    b.expn,b.x,b.y,b.z,b.I,b.J,b.K,
    c.expn,c.x,c.y,c.z,c.I,c.J,c.K,
    d.expn,d.x,d.y,d.z,d.I,d.J,d.K)
end
coulomb_hgp(a::CGBF,b::CGBF,c::CGBF,d::CGBF) = contract(coulomb_hgp,a,b,c,d)

function hrr(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,aI::Int64,aJ::Int64,aK::Int64,
  bexpn::Float64,bx::Float64,by::Float64,bz::Float64,bI::Int64,bJ::Int64,bK::Int64,
  cexpn::Float64,cx::Float64,cy::Float64,cz::Float64,cI::Int64,cJ::Int64,cK::Int64,
  dexpn::Float64,dx::Float64,dy::Float64,dz::Float64,dI::Int64,dJ::Int64,dK::Int64)
  if bI > 0
    return hrr(aexpn,ax,ay,az,aI+1,aJ,aK,bexpn,bx,by,bz,bI-1,bJ,bK,
      cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK) +
    (ax-bx)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI-1,bJ,bK,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
  elseif bJ > 0
    return hrr(aexpn,ax,ay,az,aI,aJ+1,aK,bexpn,bx,by,bz,bI,bJ-1,bK,
      cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK) +
      (ay-by)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ-1,bK,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
  elseif bK > 0
    return hrr(aexpn,ax,ay,az,aI,aJ,aK+1,bexpn,bx,by,bz,bI,bJ,bK-1,
      cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK) +
      (az-bz)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK-1,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
  elseif dI > 0
    return hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
      cexpn,cx,cy,cz,cI+1,cJ,cK,dexpn,dx,dy,dz,dI-1,dJ,dK) +
      (cx-dx)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI-1,dJ,dK)
  elseif dJ > 0
    return hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
      cexpn,cx,cy,cz,cI,cJ+1,cK,dexpn,dx,dy,dz,dI,dJ-1,dK) +
      (cy-dy)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ-1,dK)
  elseif dK > 0
    return hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
      cexpn,cx,cy,cz,cI,cJ,cK+1,dexpn,dx,dy,dz,dI,dJ,dK-1) +
      (cz-dz)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK-1)
  end
  return vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,0)
end


function vrr(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,aI::Int64,aJ::Int64,aK::Int64,
    bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
    cexpn::Float64,cx::Float64,cy::Float64,cz::Float64,cI::Int64,cJ::Int64,cK::Int64,
    dexpn::Float64,dx::Float64,dy::Float64,dz::Float64,m::Int64)
  px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
  qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,dx,dy,dz)
  zeta,eta = aexpn+bexpn,cexpn+dexpn
  wx,wy,wz = gaussian_product_center(zeta,px,py,pz,eta,qx,qy,qz)
  #println("P: $px,$py,$pz, Q: $qx,$qy,$qz, W: $wx,$wy,$wz, $zeta,$eta")

  val = 0
  if cK>0
    val = (qz-cz)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
      cexpn,cx,cy,cz,cI,cJ,cK-1,dexpn,dx,dy,dz,m) +
      (wz-qz)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK-1,dexpn,dx,dy,dz,m+1)
    #println("val1=$val")
    if cK>1
      val += 0.5*(cK-1)/eta*(
        vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
          cexpn,cx,cy,cz,cI,cJ,cK-2,dexpn,dx,dy,dz,m) -
      zeta/(zeta+eta)*
        vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
          cexpn,cx,cy,cz,cI,cJ,cK-2,dexpn,dx,dy,dz,m+1) )
    #println("val2=$val")
    end
    if aK>0
      val += 0.5*aK/(zeta+eta)*
        vrr(aexpn,ax,ay,az,aI,aJ,aK-1,bexpn,bx,by,bz,
          cexpn,cx,cy,cz,cI,cJ,cK-1,dexpn,dx,dy,dz,m+1)
    end
    #println("val3=$val")
    return val
  elseif cJ>0
    val = (qy-cy)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
      cexpn,cx,cy,cz,cI,cJ-1,cK,dexpn,dx,dy,dz,m) +
    (wy-qy)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ-1,cK-1,dexpn,dx,dy,dz,m+1)
    #println("val4=$val")
    if cJ>1
      val += 0.5*(cJ-1)/eta*(
      vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ-2,cK,dexpn,dx,dy,dz,m) -
      zeta/(zeta+eta)*
      vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ-2,cK,dexpn,dx,dy,dz,m+1)
      )
    #println("val5=$val")
    end
    if aJ>0
      val += 0.5*aJ/(zeta+eta)*
      vrr(aexpn,ax,ay,az,aI,aJ-1,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ-1,cK,dexpn,dx,dy,dz,m+1)
    end
    #println("val6=$val")
    return val
  elseif cI>0
    val = (qx-cx)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
      cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m) +
    (wx-qx)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI-1,cJ,cK-1,dexpn,dx,dy,dz,m+1)
    #println("val7=$val")
    if cI>1
      val += 0.5*(cI-1)/eta*(
      vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI-2,cJ,cK,dexpn,dx,dy,dz,m) -
      zeta/(zeta+eta)*
      vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI-2,cJ,cK,dexpn,dx,dy,dz,m+1)
      )
    end
    #println("val8=$val")
    if aI>0
      val += 0.5*aI/(zeta+eta)*
      vrr(aexpn,ax,ay,az,aI-1,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m+1)
    end
    #println("val9=$val")
    return val
  elseif aK>0
    val = (pz-az)*vrr(aexpn,ax,ay,az,aI,aJ,aK-1,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) +
    (wz-pz)*vrr(aexpn,ax,ay,az,aI,aJ,aK-1,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
    #println("val10=$val")
    if aK>1
      val += 0.5*(aK-1)/zeta*(
      vrr(aexpn,ax,ay,az,aI,aJ,aK-2,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m) -
      eta/(zeta+eta)*
      vrr(aexpn,ax,ay,az,aI,aJ,aK-2,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m+1)
      )
    end
    #println("val11=$val")
    return val
  elseif aJ>0
    val = (py-ay)*vrr(aexpn,ax,ay,az,aI,aJ-1,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m)+
    (wy-py)*vrr(aexpn,ax,ay,az,aI,aJ-1,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
    #println("val12=$val")
    if aJ>1
      val += 0.5*(aJ-1)/zeta*(
      vrr(aexpn,ax,ay,az,aI,aJ-2,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) -
      eta/(zeta+eta)*
      vrr(aexpn,ax,ay,az,aI,aJ-2,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
      )
    end
    #println("val13=$val")
    return val
  elseif aI>0
    val = (px-ax)*vrr(aexpn,ax,ay,az,aI-1,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) +
    (wx-px)*vrr(aexpn,ax,ay,az,aI-1,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
    #println("val14=$val")
    if aI>1
      val += 0.5*(aI-1)/zeta*(
      vrr(aexpn,ax,ay,az,aI-2,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) -
      eta/(zeta+eta)*
      vrr(aexpn,ax,ay,az,aI-2,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
      )
    end
    #println("val15=$val")
    return val
  end

  rab2 = dist2(ax-bx,ay-by,az-bz)
  rcd2 = dist2(cx-dx,cy-dy,cz-dz)
  rpq2 = dist2(px-qx,py-qy,pz-qz)
  T = zeta*eta/(zeta+eta)*rpq2
  Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
  Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
  #println("rab2=$rab2,rcd2=$rcd2,rpq2=$rpq2,T=$T,Kab=$Kab,Kcd=$Kcd")
  return Kab*Kcd/sqrt(zeta+eta)*Fgamma(m,T)
end


function vrr_iter(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,aI::Int64,aJ::Int64,aK::Int64,
    bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
    cexpn::Float64,cx::Float64,cy::Float64,cz::Float64,cI::Int64,cJ::Int64,cK::Int64,
    dexpn::Float64,dx::Float64,dy::Float64,dz::Float64,M::Int64)
  px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
  qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,dx,dy,dz)
  zeta,eta = aexpn+bexpn,cexpn+dexpn
  wx,wy,wz = gaussian_product_center(zeta,px,py,pz,eta,qx,qy,qz)
  rab2 = dist2(ax-bx,ay-by,az-bz)
  rcd2 = dist2(cx-dx,cy-dy,cz-dz)
  rpq2 = dist2(px-qx,py-qy,pz-qz)
  T = zeta*eta/(zeta+eta)*rpq2
  Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
  Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
  mtot = aI+aJ+aK+cI+cJ+cK+M

  vrr_terms = zeros(Float64,(aI+1,aJ+1,aK+1,cI+1,cJ+1,cK+1,mtot+1))

  for m in 0:mtot
    vrr_terms[1,1,1, 1,1,1, m+1] = Fgamma(m,T)*Kab*Kcd/sqrt(zeta+eta)
  end

  for i in 0:(aI-1)
    for m in 0:(mtot-i-1)
      vrr_terms[i+2,1,1, 1,1,1, m+1] = (
         (px-ax)*vrr_terms[i+1,1,1, 1,1,1, m+1] +
         (wx-px)*vrr_terms[i+1,1,1, 1,1,1, m+2])

      if i>0
        vrr_terms[i+2,1,1, 1,1,1, m+1] += i/2/zeta*(
          vrr_terms[i,1,1, 1,1,1, m+1] -
          eta/(zeta+eta)*vrr_terms[i,1,1, 1,1,1, m+2])
      end
    end
  end
#=
  for j in 0:(aJ-1)
    for i in 0:aI
      for m in 0:(mtot-i-j-1)
        println(("b",i,j,m))
        vrr_terms[i+1,j+2,1, 1,1,1, m+1] = (
        (py-ay)*vrr_terms[i+1,j+1,1, 1,1,1, m+1] +
        (wy-py)*vrr_terms[i+1,j+1,1, 1,1,1, m+2])
        if j>0
          vrr_terms[i+1,j+2,1, 1,1,1, m+1] += j/2/zeta*(
            vrr_terms[i+1,j,1, 1,1,1, m+1] -
            eta/(zeta+eta)*vrr_terms[i+1,j,1, 1,1,1, m+2])
        end
      end
    end
  end

  for k in 0:(aK-1)
    for j in 0:aJ
      for i in 0:aI
        for m in 0:(mtot-i-j-k-1)
          println(("c",i,j,k,m))
          vrr_terms[i+1,j+1,k+2, 1,1,1, m+1] = (
          (pz-az)*vrr_terms[i+1,j+1,k+1, 1,1,1, m+1] +
          (wz-pz)*vrr_terms[i+1,j+1,k+1, 1,1,1, m+2])
          if k>0
            vrr_terms[i+1,j+1,k+2, 1,1,1, m+1] += k/2/zeta*(
            vrr_terms[i+1,j+1,k, 1,1,1, m+1] -
            eta/(zeta+eta)*vrr_terms[i+1,j+1,k, 1,1,1, m+2])
          end
        end
      end
    end
  end

  for q in 0:(cI-1)
    for k in 0:aK
      for j in 0:aJ
        for i in 0:aI
          for m in 0:(mtot-i-j-k-q-1)
            println(("d",i,j,k,q,m))
            vrr_terms[i+1,j+1,k+1, q+2,1,1, m+1] = (
            (qx-cx)*vrr_terms[i+1,j+1,k+1, q+1,1,1, m+1] +
            (wx-qx)*vrr_terms[i+1,j+1,k+1, q+1,1,1, m+2])
            if q>0
              vrr_terms[i+1,j+1,k+1, q+2,1,1, m+1] += q/2/eta*(
              vrr_terms[i+1,j+1,k+1, q,1,1, m+1] -
              eta/(zeta+eta)*vrr_terms[i+1,j+1,k+1, q,1,1, m+2])
            end
            if i>0
              vrr_terms[i+1,j+1,k+1, q+2,1,1, m+1] += (
              i/2/(zeta+eta)*vrr_terms[i,j+1,j+1, q+1,1,1, m+2])
            end
          end
        end
      end
    end
  end

  for r in 0:(cJ-1)
    for q in 0:cI
      for k in 0:aK
        for j in 0:aJ
          for i in 0:aI
            for m in 0:(mtot-i-j-k-q-r-1)
              println(("e",i,j,k,q,r,m))
              vrr_terms[i+1,j+1,k+1, q+1,r+2,1, m+1] = (
              (qy-cy)*vrr_terms[i+1,j+1,k+1, q+1,r+1,1, m+1] +
              (wy-qy)*vrr_terms[i+1,j+1,k+1, q+1,r+1,1, m+2])
              if r>0
                vrr_terms[i+1,j+1,k+1, q+1,r+2,1, m+1] += r/2/eta*(
                vrr_terms[i+1,j+1,k+1, q+1,r+1,1, m+1] -
                zeta/(zeta+eta)*vrr_terms[i+1,j+1,k+1, q+1,r,1, m+2])
              end
              if j>0
                vrr_terms[i+1,j+1,k+1, q+1,r+2,1, m+1] += (
                j/2/(zeta+eta)*vrr_terms[i+1,j,k+1, q+1,r+1,1, m+2])
              end
            end
          end
        end
      end
    end
  end

  for s in 0:(cK-1)
    for r in 0:cJ
      for q in 0:cI
        for k in 0:aK
          for j in 0:aJ
            for i in 0:aI
              for m in 0:(mtot-i-j-k-q-r-s-1)
                println(("f",i,j,k,q,r,s,m))
                vrr_terms[i+1,j+1,k+1, q+1,r+1,s+2, m+1] = (
                (qz-cz)*vrr_terms[i+1,j+1,k+1, q+1,r+1,s+1, m+1]+
                (wz-qz)*vrr_terms[i+1,j+1,k+1, q+1,r+1,s+1, m+2])
                if s>0
                  vrr_terms[i+1,j+1,k+1, q+1,r+1,s+2,m+1] += s/2/eta*(
                  vrr_terms[i+1,j+1,k+1,q+1,r+1,s, m+1] -
                  zeta/(zeta+eta)*vrr_terms[i+1,j+1,k+1,q+1,r+1,s,m+2])
                end
                if k>0
                  vrr_terms[i+1,j+1,k+1, q+1,r+1,s+2,m+1] += (
                  k/2/(zeta+eta)*vrr_terms[i+1,j+1,k, q+1,r+1,s+1,m+2])
                end
              end
            end
          end
        end
      end
    end
  end
  println("before return")
  =#
  @show vrr_terms[aI+1,aJ+1,aK+1,cI+1,cJ+1,cK+1,M+1]
  vrr_terms[aI+1,aJ+1,aK+1,cI+1,cJ+1,cK+1,M+1]
end

function test_vrr()
  ax=ay=az=bx=by=bz=cx=cy=cz=dx=dy=dz=0.0
  aexpn=bexpn=cexpn=dexpn=1.0
  aI=aJ=aK=0
  cI=cJ=cK=0
  M=0

  for (ax,ay,az, aI,aJ,aK, cI,cJ,cK, result) in Any[
      (0.,0.,0., 0,0,0, 0,0,0, 4.37335456733),
      (0.,0.,0., 1,0,0, 1,0,0, 0.182223107579),
      (0.,0.,0., 0,1,0, 0,1,0, 0.182223107579),
      (0.,0.,0., 0,0,1, 0,0,1, 0.182223107579),

      (0.,0.,0., 2,0,0, 2,0,0,  0.223223306785),
      (0.,0.,0., 0,2,0, 0,2,0,  0.223223306785),
      (0.,0.,0., 0,0,2, 0,0,2,  0.223223306785),

      (1.,2.,3., 1,0,0, 1,0,0, -5.63387712455e-06),
      (1.,2.,3., 0,1,0, 0,1,0, -0.000116463120359),
      (1.,2.,3., 0,0,1, 0,0,1, -0.000301178525749),

      (1.,2.,3., 2,0,0, 2,0,0, 0.000225033081978),
      (1.,2.,3., 0,2,0, 0,2,0, 0.000610247078796),
      (1.,2.,3., 0,0,2, 0,0,2, 0.00134278307956),

      (0.,0.,0., 1,1,0, 1,1,0, 0.0136667330685),
      (0.,0.,0., 0,1,1, 0,1,1, 0.0136667330685),
      (0.,0.,0., 1,0,1, 1,0,1, 0.0136667330685),

      (3.,2.,1., 1,1,0, 1,1,0, 5.97677147819e-05),
      (3.,2.,1., 0,1,1, 0,1,1, 1.57429039496e-06),
      (3.,2.,1., 1,0,1, 1,0,1, 4.00292836291e-06)
    ]

    val1 = vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
      cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,M)
    val2 = vrr(cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,
      aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,M)
    @assert isapprox(val1,val2)
    #@assert isapprox(val1,result)
    println(val1, " ", result)
    val3 = vrr_iter(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
      cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,M)
    val4 = vrr_iter(cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,
      aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,M)
    @show (val1,val2,val3,val4)
  end
end


function test_hrr()

  @printf("test_hrr()\n")

  ax=ay=az=bx=by=bz=cx=cy=cz=dx=dy=dz=0.0
  aexpn=bexpn=cexpn=dexpn=1.0
  aI=aJ=aK=0
  bI,bJ,bK = 1,0,1
  cI=cJ=cK=0
  dI,dJ,dK = 1,0,1


  for (ax,ay,az, aI,aJ,aK, cI,cJ,cK, result) in Any[
      (0.,0.,0., 0,0,0, 0,0,0, 0.0136667330685),
      (0.,0.,0., 1,0,0, 1,0,0, 0.00821630976139),
      (0.,0.,0., 0,1,0, 0,1,0, 0.00122024402397),
      (0.,0.,0., 0,0,1, 0,0,1, 0.00821630976139),

      (0.,0.,0., 2,0,0, 2,0,0,   0.0039759617781),
      (0.,0.,0., 0,2,0, 0,2,0,   0.000599953311785),
      (0.,0.,0., 0,0,2, 0,0,2,  0.0039759617781),

      (1.,2.,3., 1,0,0, 1,0,0, -1.18513079462e-06),
      (1.,2.,3., 0,1,0, 0,1,0,  -4.66999128258e-06),
      (1.,2.,3., 0,0,1, 0,0,1, -3.47437366868e-05),

      (1.,2.,3., 2,0,0, 2,0,0, 2.81002247462e-06),
      (1.,2.,3., 0,2,0, 0,2,0, 7.09856891538e-06),
      (1.,2.,3., 0,0,2, 0,0,2, 3.62153023224e-05),

      (0.,0.,0., 1,1,0, 1,1,0, 0.000599953311785),
      (0.,0.,0., 0,1,1, 0,1,1, 0.000599953311785),
      (0.,0.,0., 1,0,1, 1,0,1, 0.0116431617287),

      (3.,2.,1., 1,1,0, 1,1,0, 7.37307761485e-06),
      (3.,2.,1., 0,1,1, 0,1,1, 2.53332441858e-07),
      (3.,2.,1., 1,0,1, 1,0,1, 2.4521155336e-06)
    ]
    #println("hrr($aexpn,$ax,$ay,$az,$aI,$aJ,$aK,$bexpn,$bx,$by,$bz,$bI,$bJ,$bK,")
    #println("  $cexpn,$cx,$cy,$cz,$cI,$cJ,$cK,$dexpn,$dx,$dy,$dz,$dI,$dJ,$dK)")
    val1 = hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
      cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
    val2 = hrr(cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK,
      aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)
    #@show val1,val2,result
    @assert isapprox(val1,val2)
    @assert isapprox(val1,result)
  end

  @printf("Pass test_hrr\n")
end
