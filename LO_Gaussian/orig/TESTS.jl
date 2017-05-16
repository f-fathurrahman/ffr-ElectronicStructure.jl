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

function test_overlap()
  s = pgbf(1.0)
  px = pgbf(1.0,0,0,0,1,0,0)
  @assert overlap1d(0,0,0.,0.,1.) == 1
  @assert gaussian_product_center(s,s) == [0,0,0]
  @assert isapprox(overlap(s,s),1)
  @assert isapprox(overlap(px,px),1)
  @assert isapprox(overlap(s,px),0)
  @assert binomial_prefactor(0,0,0,0.,0.) == 1
end

function test_utils()
  @assert factorial2(6)==48
  @assert collect(pairs(3)) == Any[(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)]
  @assert collect(pairs(3,"subdiag")) == Any[(1,2),(1,3),(2,3)]
  @assert collect(pairs(2,"rect")) == Any[(1,1),(1,2),(2,1),(2,2)]
  @assert iindex(1,1,1,1) == 1
  @assert iindex(1,1,1,2) == iindex(1,1,2,1) == iindex(1,2,1,1) == iindex(2,1,1,1) == 2
  @assert iindex(1,1,2,2) == iindex(2,2,1,1) == 4
  @printf "test_utils is passed\n"
end
