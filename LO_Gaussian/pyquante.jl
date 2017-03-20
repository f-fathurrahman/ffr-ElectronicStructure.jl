"""
PyQuante in Julia
Experimenting with writing quantum chemistry in Julia

"""




# ## Atoms and Molecules

nuclear_repulsion(a::Atom,b::Atom)= a.atno*b.atno/sqrt(dist2(a.x-b.x,a.y-b.y,a.z-b.z))
function nuclear_repulsion(mol::Molecule)
  nr = 0
  for (i,j) in pairs(nat(mol),"subdiag")
    nr += nuclear_repulsion(mol.atomlist[i],mol.atomlist[j])
  end
  return nr
end



function test_geo_basis()
  @assert isapprox(nuclear_repulsion(h2),0.7223600367)
  @assert nel(h2) == 2
  @assert nel(h2o) == 10
  @assert length(sto3g)==10
  bfs = build_basis(h2)
  @assert length(bfs.bfs)==2
  l,r = bfs.bfs
  @assert isapprox(overlap(l,l),1)
  @assert isapprox(overlap(r,r),1)
  @assert isapprox(overlap(l,r),0.66473625)
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

function all_1e_ints(bfs::BasisSet,mol::Molecule)
  n = length(bfs.bfs)
  S = Array(Float64,(n,n))
  T = Array(Float64,(n,n))
  V = Array(Float64,(n,n))
  for (i,j) in pairs(n)
    a,b = bfs.bfs[i],bfs.bfs[j]
    S[i,j] = S[j,i] = overlap(a,b)
    T[i,j] = T[j,i] = kinetic(a,b)
    V[i,j] = V[j,i] = nuclear_attraction(a,b,mol)
  end
  return S,T,V
end

function all_twoe_ints(bflist,ERI=coulomb)
  n = length(bflist.bfs)
  totlen = div(n*(n+1)*(n*n+n+2),8)
  ints2e = Array(Float64,totlen)
  for (i,j,k,l) in iiterator(n)
    ints2e[iindex(i,j,k,l)] = ERI(bflist.bfs[i],bflist.bfs[j],bflist.bfs[k],bflist.bfs[l])
  end
  return ints2e
end

function make2JmK(D::Array{Float64,2},Ints::Array{Float64,1})
  n = size(D,1)
  G = Array(Float64,(n,n))
  D1 = reshape(D,n*n)
  temp = Array(Float64,n*n)
  for (i,j) in pairs(n)
    kl = 1
    for (k,l) in pairs(n,"rect")
      temp[kl] = 2*Ints[iindex(i,j,k,l)]-Ints[iindex(i,k,j,l)]
      kl += 1
    end
    G[i,j] = G[j,i] = dot(D1,temp)
  end
  return G
end

dmat(U::Array{Float64,2},nocc::Int64) = U[:,1:nocc]*U[:,1:nocc]'


function rhf(mol::Molecule,MaxIter::Int64=8,verbose::Bool=false)
  bfs = build_basis(mol)
  S,T,V = all_1e_ints(bfs,mol)
  Ints = all_twoe_ints(bfs)
  h = T+V
  E,U = eig(h,S)
  Enuke = nuclear_repulsion(mol)
  nclosed,nopen = divrem(nel(mol),2)
  Eold = 0
  Energy = 0
  println("Nel=$(nel(mol)) Nclosed=$nclosed")
  if verbose
    println("S=\n$S")
    println("h=\n$h")
    println("T=\n$T")
    println("V=\n$V")
    println("E: $E")
    println("U: $U")
    println("2e ints:\n$Ints")
  end
  for iter in 1:MaxIter
    D = dmat(U,nclosed)
    if verbose
      println("D=\n$D")
    end
    G = make2JmK(D,Ints)
    H = h+G
    E,U = eig(H,S)
    Eone = trace2(D,h)
    Etwo = trace2(D,H)
    Energy = Enuke + Eone + Etwo
    println("HF: $iter  $Energy : $Enuke  $Eone  $Etwo")
    if isapprox(Energy,Eold)
      break
    end
    Eold  = Energy
  end
  return Energy,E,U
end

function test_h2()
  @time Energy, E, U = rhf(h2)
  @assert isapprox(Energy,-1.1170996)
end

function test_lih()
  @time Energy, E, U = rhf(lih)
  @assert isapprox(Energy,-7.86073270525799)
end

function test_h2o()
  @time Energy,E,U = rhf(h2o)
  @assert isapprox(Energy,-74.9597609118851)
end

function test()
  test_utils()
  test_pgbf()
  test_cgbf()
  test_overlap()
  test_kinetic()
  test_a_terms()
  test_gamma()
  test_na()
  test_fgamma()
  test_one()
  test_na2()
  test_two_terms()
  test_coul1()
  test_vrr()
  test_hrr()
  test_geo_basis()
  test_h2()
  test_lih()
  test_h2o()
end

test()
