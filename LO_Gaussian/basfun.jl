# ## Basis Set Data
# Note use of curly braces here. Julia assumes that if you have square braces, you want
# things flattened as much as possible (to be as fast as possible, I guess). Curlys
# preserve the list structure the way I would expect from Python

include("sto3g.jl")

# ## Basis function definitions

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

function pgbf(expn,x=0,y=0,z=0,I=0,J=0,K=0,norm=1)
  p = PGBF(expn,x,y,z,I,J,K,norm)
  normalize!(p)
  return p
end

function amplitude(bf::PGBF,x,y,z)
  dx,dy,dz = x-bf.x,y-bf.y,z-bf.z
  r2 = dist2(dx,dy,dz)
  return bf.norm*(dx^bf.I)*(dy^bf.J)*(dz^bf.K)*exp(-bf.expn*r2)
end

function normalize!(pbf::PGBF)
  pbf.norm /= sqrt(overlap(pbf,pbf))
end


type CGBF
  x::Float64
  y::Float64
  z::Float64
  I::Int64
  J::Int64
  K::Int64
  norm::Float64
  pgbfs::Array{PGBF,1}
  coefs::Array{Float64,1}
end

cgbf(x=0,y=0,z=0,I=0,J=0,K=0) = CGBF(x,y,z,I,J,K,1.0,PGBF[],Float64[])

function amplitude(bf::CGBF,x,y,z)
  s = 0
  for (c,pbf) in primitives(bf)
    s += c*amplitude(pbf,x,y,z)
  end
  return bf.norm*s
end

function normalize!(bf::CGBF)
  bf.norm /= sqrt(overlap(bf,bf))
end

primitives(a::CGBF) = zip(a.coefs,a.pgbfs)

function push!(cbf::CGBF,expn,coef)
  Base.push!(cbf.pgbfs,pgbf(expn,cbf.x,cbf.y,cbf.z,cbf.I,cbf.J,cbf.K))
  Base.push!(cbf.coefs,coef)
  normalize!(cbf)
end

function contract(f,a::CGBF,b::CGBF)
  s = 0
  for (ca,abf) in primitives(a)
    for (cb,bbf) in primitives(b)
      s += ca*cb*f(abf,bbf)
    end
  end
  return a.norm*b.norm*s
end

function contract(f,a::CGBF,b::CGBF,c::CGBF,d::CGBF)
  s = 0
  for (ca,abf) in primitives(a)
    for (cb,bbf) in primitives(b)
      for (cc,cbf) in primitives(c)
        for (cd,dbf) in primitives(d)
          s += ca*cb*cc*cd*f(abf,bbf,cbf,dbf)
        end
      end
    end
  end
  return a.norm*b.norm*c.norm*d.norm*s
end


function test_pgbf()
  s = pgbf(1.0)
  px = pgbf(1.0,0,0,0,1,0,0)
  @assert isapprox(amplitude(s,0,0,0),0.71270547)
  @assert isapprox(amplitude(px,0,0,0),0)
  @printf "test_pgbf is passed\n"
end

function test_cgbf()
  c = cgbf(0.0,0.0,0.0)
  push!(c,1,1)
  @assert isapprox(amplitude(c,0,0,0),0.71270547)
  c2 = cgbf(0,0,0)
  push!(c2,1,0.2)
  push!(c2,0.5,0.2)
  @assert isapprox(overlap(c2,c2),1)
  @printf "test_cgbf is passed\n"
end

type BasisSet # list of CGBFs
  bfs::Array{CGBF,1}
end

basisset() = BasisSet(CGBF[])

function push!(basis::BasisSet,cbf::CGBF)
  Base.push!(basis.bfs,cbf)
end

function build_basis(mol::Molecule,name="sto3g")
  data = basis_set_data[name]
  basis_set = basisset()
  for atom in mol.atomlist
    for btuple in data[atom.atno]
      sym,primlist = btuple
      for (I,J,K) in sym2power[sym]
        cbf = cgbf(atom.x,atom.y,atom.z,I,J,K)
        push!(basis_set,cbf)
        for (expn,coef) in primlist
          push!(cbf,expn,coef)
        end
      end
    end
  end
  return basis_set
end

const sym2power = Dict{Any,Any}(
  'S' => [(0,0,0)],
  'P' => [(1,0,0),(0,1,0),(0,0,1)],
  'D' => [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]
  )
