"""
### Nuclear attraction term
"""
function Aterm(i::Int64,r::Int64,u::Int64,l1::Int64,l2::Int64,ax::Float64,bx::Float64,cx::Float64,gamma::Float64)
  term1 = (-1)^(i+u)*binomial_prefactor(i,l1,l2,ax,bx)
  term2 = factorial(i)*cx^(i-2r-2u)
  term3 = (1/4/gamma)^(r+u)/factorial(r)/factorial(u)/factorial(i-2r-2u)
  return term1*term2*term3
end

function Aarray(l1::Int64,l2::Int64,a::Float64,b::Float64,c::Float64,g::Float64)
  Imax = l1+l2+1
  A = zeros(Float64,Imax)
  for i in 0:(Imax-1)
    for r in 0:div(i,2)
      for u in 0:div(i-2r,2)
        I = i-2r-u+1
        A[I] += Aterm(i,r,u,l1,l2,a,b,c,g)
      end
    end
  end
  return A
end

function nuclear_attraction(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,
              aI::Int64,aJ::Int64,aK::Int64,
              bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
              bI::Int64,bJ::Int64,bK::Int64,
              cx::Float64,cy::Float64,cz::Float64)
  #print("na($aexpn,$ax,$ay,$az,$aI,$aJ,$aK,$bexpn,$bx,$by,$bz,$bI,$bJ,$bK,$cx,$cy,$cz)=")
  px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
  gamma = aexpn+bexpn
  rab2 = dist2(ax-bx,ay-by,az-bz)
  rcp2 = dist2(cx-px,cy-py,cz-pz)
  Ax = Aarray(aI,bI,px-ax,px-bx,px-cx,gamma)
  Ay = Aarray(aJ,bJ,py-ay,py-by,py-cy,gamma)
  Az = Aarray(aK,bK,pz-az,pz-bz,pz-cz,gamma)
  total = 0
  for I in 0:(aI+bI)
    for J in 0:(aJ+bJ)
      for K in 0:(aK+bK)
        total += Ax[I+1]*Ay[J+1]*Az[K+1]*Fgamma(I+J+K,rcp2*gamma)
      end
    end
  end
  val=-2pi*exp(-aexpn*bexpn*rab2/gamma)*total/gamma
  #println(val)
  #println((Ax,Ay,Az,rcp2*gamma,Fgamma(0,rcp2*gamma)))
  return val
end

function nuclear_attraction(a::PGBF,b::PGBF,cx::Float64,cy::Float64,cz::Float64)
  return a.norm*b.norm*nuclear_attraction(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
                      b.expn,b.x,b.y,b.z,b.I,b.J,b.K,cx,cy,cz)
end
nuclear_attraction(a::PGBF,b::PGBF,c::Atom) = c.atno*nuclear_attraction(a,b,c.x,c.y,c.z)
nuclear_attraction(a::PGBF,b::PGBF,m::Molecule) = sum([nuclear_attraction(a,b,c) for c in m.atomlist])

function Fgamma(m::Int64,x::Float64,SMALL::Float64=1e-12)
  #println("Fgamma($m,$x)")
  x = max(x,SMALL) # Evidently needs underflow protection
  return 0.5*x^(-m-0.5)*gammainc(m+0.5,x)
end

function gammainc(a::Float64,x::Float64)
  # This is the series version of gamma from pyquante. For reasons I don't get, it
  # doesn't work around a=1. This works alright, but is only a stopgap solution
  # until Julia gets an incomplete gamma function programmed
  if abs(a-1) < 1e-3
    println("Warning: gammainc_series is known to have problems for a ~ 1")
  end
  if x < (a+1.0)
    #Use the series representation
    gam,gln = gser(a,x)
  else
    #Use continued fractions
    gamc,gln = gcf(a,x)
    gam = 1-gamc
  end
  return exp(gln)*gam
end

function gser(a::Float64,x::Float64,ITMAX::Int64=100,EPS::Float64=3e-9)
  # Series representation of Gamma. NumRec sect 6.1.
  gln=lgamma(a)
  if x == 0
    return 0,gln
  end
  ap = a
  delt = s = 1/a
  for i in 1:ITMAX
    ap += 1
    delt *= (x/ap)
    s += delt
    if abs(delt) < abs(s)*EPS
      break
    end
  end
  return s*exp(-x+a*log(x)-gln),gln
end

function gcf(a::Float64,x::Float64,ITMAX::Int64=200,EPS::Float64=3e-9,FPMIN::Float64=1e-30)
  #Continued fraction representation of Gamma. NumRec sect 6.1"
  gln=lgamma(a)
  b=x+1.-a
  c=1./FPMIN
  d=1./b
  h=d
  for i in 1:ITMAX
    an=-i*(i-a)
    b=b+2.
    d=an*d+b
    if abs(d) < FPMIN
      d=FPMIN
    end
    c=b+an/c
    if abs(c) < FPMIN
      c=FPMIN
    end
    d=1./d
    delt=d*c
    h=h*delt
    if abs(delt-1.) < EPS
      break
    end
  end
  gammcf = exp(-x+a*log(x)-gln)*h
  return gammcf,gln
end


# Need a nested scope to squeeze this into the contract function
function nuclear_attraction(a::CGBF,b::CGBF,cx::Float64,cy::Float64,cz::Float64)
  na(a,b) = nuclear_attraction(a,b,cx,cy,cz)
  contract(na,a,b)
end
function nuclear_attraction(a::CGBF,b::CGBF,c::Atom)
  na(a,b) = nuclear_attraction(a,b,c)
  contract(na,a,b)
end
function nuclear_attraction(a::CGBF,b::CGBF,m::Molecule)
  na(a,b) = nuclear_attraction(a,b,m)
  contract(na,a,b)
end
