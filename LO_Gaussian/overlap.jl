"""
## One-electron integrals
### Overlap matrix elements
"""

function overlap( a::PGBF, b::PGBF )
    s = overlap( a.expn, a.center, a.power, b.expn, b.center, b.power)
    return a.NORM*b.NORM*s
end

# overlap with itself
function overlap(aexpn::Float64, acenter::Tuple3F64, apower::Tuple3I64)
    return overlap(aexpn,acenter,apower,aexpn,acenter,apower)
end

# Low level method
function overlap(aexpn::Float64, acenter::Tuple3F64, apower::Tuple3I64,
    bexpn::Float64, bcenter::Tuple3F64, bpower::Tuple3I64)

    ax, ay, az = acenter[1],acenter[2],acenter[3]
    bx, by, bz = bcenter[1],bcenter[2],bcenter[3]

    aI,aJ,aK = apower[1],apower[2],apower[3]
    bI,bJ,bK = bpower[1],bpower[2],bpower[3]

    gamma = aexpn + bexpn
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    rab2 = dist2(ax-bx,ay-by,az-bz)
    pre = (pi/gamma)^1.5*exp(-aexpn*bexpn*rab2/gamma)
    wx = overlap1d(aI,bI,px-ax,px-bx,gamma)
    wy = overlap1d(aJ,bJ,py-ay,py-by,gamma)
    wz = overlap1d(aK,bK,pz-az,pz-bz,gamma)
    return pre*wx*wy*wz
end

function gaussian_product_center(a::PGBF,b::PGBF)
  return (a.expn*[a.x,a.y,a.z]+b.expn*[b.x,b.y,b.z])/(a.expn+b.expn)
end

function gaussian_product_center(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,
                  bexpn::Float64,bx::Float64,by::Float64,bz::Float64)
  return (aexpn*[ax,ay,az]+bexpn*[bx,by,bz])/(aexpn+bexpn)
end

function overlap1d(la::Int64,lb::Int64,ax::Float64,bx::Float64,gamma::Float64)
  total = 0
  for i in 0:div(la+lb,2)
    total += binomial_prefactor(2i,la,lb,ax,bx)*factorial2(2i-1)/(2gamma)^i
  end
  return total
end

function binomial_prefactor(s::Int64,ia::Int64,ib::Int64,xpa::Float64,xpb::Float64)
  total = 0
  for t in 0:s
    if (s-ia) <= t <= ib
      total += binomial(ia,s-t)*binomial(ib,t)*xpa^(ia-s+t)*xpb^(ib-t)
    end
  end
  return total
end

# where should we put this ??
overlap(a::CGBF,b::CGBF) = contract(overlap,a,b)
