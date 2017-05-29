include("PGBF.jl")
include("CGBF.jl")
include("overlap.jl")
include("kinetic.jl")
include("utils.jl")

function test_kinetic()
  s = init_PGBF(1.0)
  c = init_CGBF(0.0,0.0,0.0)
  push!(c,1,1)
  @assert isapprox(amplitude(c,0,0,0),0.71270547)
  @assert isapprox(kinetic(1.,0.,0.,0.,0,0,0,1.,0.,0.,0.,0,0,0),2.9530518648229536)
  @assert isapprox(kinetic(s,s),1.5)
  @assert isapprox(kinetic(c,c),1.5)
  @printf "test_kinetic is passed\n"
end

test_kinetic()
