# ## Basis Set Data
# Note use of curly braces here. Julia assumes that if you have square braces, you want
# things flattened as much as possible (to be as fast as possible, I guess). Curlys
# preserve the list structure the way I would expect from Python

const sto3g = Any[
  # H
  [('S',
    [(3.4252509099999999, 0.15432897000000001),
     (0.62391373000000006, 0.53532813999999995),
     (0.16885539999999999, 0.44463454000000002)])],
  # He
  [('S',
    [(6.3624213899999997, 0.15432897000000001),
     (1.1589229999999999, 0.53532813999999995),
     (0.31364978999999998, 0.44463454000000002)])],
  # Li
  [('S',
    [(16.119575000000001, 0.15432897000000001),
     (2.9362007000000001, 0.53532813999999995),
     (0.79465050000000004, 0.44463454000000002)]),
   ('S',
    [(0.63628969999999996, -0.099967230000000004),
     (0.14786009999999999, 0.39951282999999999),
     (0.048088699999999998, 0.70011546999999996)]),
   ('P',
    [(0.63628969999999996, 0.15591627),
     (0.14786009999999999, 0.60768372000000004),
     (0.048088699999999998, 0.39195739000000002)])],
  # Be
  [('S',
    [(30.167871000000002, 0.15432897000000001),
     (5.4951153000000001, 0.53532813999999995),
     (1.4871927, 0.44463454000000002)]),
   ('S',
    [(1.3148331, -0.099967230000000004),
     (0.3055389, 0.39951282999999999),
     (0.099370700000000006, 0.70011546999999996)]),
   ('P',
    [(1.3148331, 0.15591627),
     (0.3055389, 0.60768372000000004),
     (0.099370700000000006, 0.39195739000000002)])],
  # B
  [('S',
    [(48.791113000000003, 0.15432897000000001),
     (8.8873622000000001, 0.53532813999999995),
     (2.4052669999999998, 0.44463454000000002)]),
   ('S',
    [(2.2369561, -0.099967230000000004),
     (0.51982050000000002, 0.39951282999999999),
     (0.16906180000000001, 0.70011546999999996)]),
   ('P',
    [(2.2369561, 0.15591627),
     (0.51982050000000002, 0.60768372000000004),
     (0.16906180000000001, 0.39195739000000002)])],
  # C
  [('S',
    [(71.616837000000004, 0.15432897000000001),
     (13.045095999999999, 0.53532813999999995),
     (3.5305122, 0.44463454000000002)]),
   ('S',
    [(2.9412493999999998, -0.099967230000000004),
     (0.68348310000000001, 0.39951282999999999),
     (0.22228990000000001, 0.70011546999999996)]),
   ('P',
    [(2.9412493999999998, 0.15591627),
     (0.68348310000000001, 0.60768372000000004),
     (0.22228990000000001, 0.39195739000000002)])],
  # N
  [('S',
    [(99.106168999999994, 0.15432897000000001),
     (18.052312000000001, 0.53532813999999995),
     (4.8856602000000002, 0.44463454000000002)]),
   ('S',
    [(3.7804559000000002, -0.099967230000000004),
     (0.87849659999999996, 0.39951282999999999),
     (0.28571439999999998, 0.70011546999999996)]),
   ('P',
    [(3.7804559000000002, 0.15591627),
     (0.87849659999999996, 0.60768372000000004),
     (0.28571439999999998, 0.39195739000000002)])],
  # O
  [('S',
    [(130.70931999999999, 0.15432897000000001),
     (23.808861, 0.53532813999999995),
     (6.4436083000000002, 0.44463454000000002)]),
   ('S',
    [(5.0331513000000001, -0.099967230000000004),
     (1.1695960999999999, 0.39951282999999999),
     (0.38038899999999998, 0.70011546999999996)]),
   ('P',
    [(5.0331513000000001, 0.15591627),
     (1.1695960999999999, 0.60768372000000004),
     (0.38038899999999998, 0.39195739000000002)])],
  # F
  [('S',
    [(166.67912999999999, 0.15432897000000001),
     (30.360811999999999, 0.53532813999999995),
     (8.2168206999999995, 0.44463454000000002)]),
   ('S',
    [(6.4648032000000004, -0.099967230000000004),
     (1.5022812000000001, 0.39951282999999999),
     (0.48858849999999998, 0.70011546999999996)]),
   ('P',
    [(6.4648032000000004, 0.15591627),
     (1.5022812000000001, 0.60768372000000004),
     (0.48858849999999998, 0.39195739000000002)])],
  # Ne
  [('S',
     [(207.01561000000001, 0.15432897000000001),
    (37.708151000000001, 0.53532813999999995),
    (10.205297, 0.44463454000000002)]),
    ('S',
     [(8.2463151000000003, -0.099967230000000004),
    (1.9162661999999999, 0.39951282999999999),
    (0.62322929999999999, 0.70011546999999996)]),
    ('P',
     [(8.2463151000000003, 0.15591627),
    (1.9162661999999999, 0.60768372000000004),
      (0.62322929999999999, 0.39195739000000002)])]
]
const basis_set_data = Dict{Any,Any}("sto3g" => sto3g);


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
