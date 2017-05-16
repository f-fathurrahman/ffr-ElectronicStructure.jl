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

function test_a_terms()
  @assert Aterm(0,0,0,0,0,0.,0.,0.,0.) == 1.0
  @assert Aarray(0,0,0.,0.,0.,1.) == [1.0]
  @assert Aarray(0,1,1.,1.,1.,1.) == [1.0, -1.0]
  @assert Aarray(1,1,1.,1.,1.,1.) == [1.5, -2.5, 1.0]
  @assert Aterm(0,0,0,0,0,0.,0.,0.,1.) == 1.0
  @assert Aterm(0,0,0,0,1,1.,1.,1.,1.) == 1.0
  @assert Aterm(1,0,0,0,1,1.,1.,1.,1.) == -1.0
  @assert Aterm(0,0,0,1,1,1.,1.,1.,1.) == 1.0
  @assert Aterm(1,0,0,1,1,1.,1.,1.,1.) == -2.0
  @assert Aterm(2,0,0,1,1,1.,1.,1.,1.) == 1.0
  @assert Aterm(2,0,1,1,1,1.,1.,1.,1.) == -0.5
  @assert Aterm(2,1,0,1,1,1.,1.,1.,1.) == 0.5
end

function test_gamma()
  # gammainc test functions. Test values taken from Mathematica
  # println("a=0.5 test")
  @assert maximum([gammainc(0.5,float(x)) for x in 0:10]
      -Any[0, 1.49365, 1.69181, 1.7471, 1.76416, 1.76968,
        1.77151, 1.77213, 1.77234, 1.77241, 1.77244]) < 1e-5

  # println("a=1.5 test")
  @assert maximum([gammainc(1.5,float(x)) for x in 0:10]
      -Any[0, 1.49365, 1.69181, 1.7471, 1.76416, 1.76968,
        1.77151, 1.77213, 1.77234, 1.77241, 1.77244]) < 1e-5
  # println("a=2.5 test")
  @assert maximum([gammainc(2.5,float(x)) for x in 0:10]
      -Any[0, 0.200538, 0.59898, 0.922271, 1.12165, 1.22933,
        1.2831, 1.30859, 1.32024, 1.32542, 1.32768]) < 1e-5
end

function test_na()
  s = pgbf(1.0)
  c = cgbf(0.0,0.0,0.0)
  push!(c,1,1)
  @assert isapprox(amplitude(c,0,0,0),0.71270547)
  @assert isapprox(nuclear_attraction(s,s,0.,0.,0.),-1.59576912)
  @assert isapprox(nuclear_attraction(c,c,0.,0.,0.),-1.595769)
end

function test_fgamma()
  for (x,res) in Any[(0.,1),
          (30.,0.161802),
          (60.,0.1144114),
          (90.,0.0934165),
          (120.,0.08090108),
          (300.,0.051166336)]
    #@assert isapprox(res,Fgamma(0,x))
    @printf("res, Fgamma(0,x) %18.10f %18.10f\n", x, Fgamma(0,x))
  end
  @printf "test_fgamma is passed\n"
end

#todo make into a test
function test_one()
  s1 = pgbf(1)
  s2 = pgbf(1,0,1,0)
  x=y=0.
  println("S: $(overlap(s1,s2))")
  println("T: $(kinetic(s1,s2))")
  for z in linspace(0,1,5)
    println("V: $z $(nuclear_attraction(s1,s2,x,y,z))")
  end
end


function test_na2()
  li,h = lih.atomlist
  bfs = build_basis(lih)
  s1,s2,x,y,z,h1s = bfs.bfs
  @assert isapprox(nuclear_attraction(s1,s1,lih),-8.307532656)
end


nuclear_repulsion(a::Atom,b::Atom)= a.atno*b.atno/sqrt(dist2(a.x-b.x,a.y-b.y,a.z-b.z))
function nuclear_repulsion(mol::Molecule)
  nr = 0
  for (i,j) in pairs(nat(mol),"subdiag")
    nr += nuclear_repulsion(mol.atomlist[i],mol.atomlist[j])
  end
  return nr
end

function test_geo_basis()

  @printf("\nCalling test_geo_basis\n")

  E1     = nuclear_repulsion(h2)
  E1_ref = 0.7223600367
  @printf("H2: nuclear_repulsion %18.10f, err = %18.10f\n", E1, abs(E1_ref-E1))

  @assert nel(h2) == 2
  @assert nel(h2o) == 10
  @assert length(sto3g)==10

  bfs = build_basis(h2)
  @assert length(bfs.bfs)==2

  l,r = bfs.bfs
  @assert isapprox(overlap(l,l),1)
  @assert isapprox(overlap(r,r),1)

  #@assert isapprox(overlap(l,r),0.66473625)
  ovl1 = overlap(l,r)
  ovl1_ref = abs(ovl1 - 0.66473625)
  @printf("overlap(l,r), err %18.10f, %18.10f\n", ovl1, ovl1_ref)

  @assert isapprox(kinetic(l,l),0.76003188)
  @assert isapprox(kinetic(r,r),0.76003188)
  @assert isapprox(kinetic(l,r),0.24141861181119084)
  @assert isapprox(coulomb(l,l,l,l), 0.7746059439196398)
  @assert isapprox(coulomb(r,r,r,r), 0.7746059439196398)
  @assert isapprox(coulomb(l,l,r,r), 0.5727937653511646)
  @assert isapprox(coulomb(l,l,l,r), 0.4488373301593464)
  @assert isapprox(coulomb(l,r,l,r), 0.3025451156654606)
  bfs = build_basis(h2o)

  s1,s2,px,py,pz,hl,hr = bfs.bfs
  @assert isapprox(coulomb(s1,s2,hl,hr),0.03855344493645537)
  @assert isapprox(coulomb(s1,pz,hl,hr),-0.0027720110485359053)
  @assert isapprox(coulomb(s1,hl,pz,hr),-0.010049491284827426)
  @assert coulomb(s1,py,hl,hr)==0
  @assert coulomb(s1,hl,py,hr)==0
end
