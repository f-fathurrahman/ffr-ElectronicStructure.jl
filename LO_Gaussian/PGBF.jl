type PGBF
  expn::Float64
  x::Float64
  y::Float64
  z::Float64
  I::Int64
  J::Int64
  K::Int64
  norm::Float64
end

# Enable keyword-based
function init_PGBF(expn; x=0,y=0,z=0,I=0,J=0,K=0,norm=1)
  p = PGBF(expn,x,y,z,I,J,K,norm)
  normalize!(p)
  return p
end

function init_PGBF(expn, x=0,y=0,z=0, I=0, J=0, K=0, norm=1 )
  p = PGBF(expn,x,y,z,I,J,K,norm)
  normalize!(p)
  return p
end

function amplitude( bf::PGBF, x, y, z)
  #
  dx, dy, dz = x-bf.x, y-bf.y, z-bf.z
  r2 = dist2(dx, dy, dz)
  #
  return bf.norm*(dx^bf.I)*(dy^bf.J)*(dz^bf.K)*exp(-bf.expn*r2)
end

function normalize!( pbf::PGBF )
  pbf.norm /= sqrt( overlap(pbf,pbf) )
end
