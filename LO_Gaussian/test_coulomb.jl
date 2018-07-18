using Printf
using SpecialFunctions

include("constants.jl")
include("alias.jl")
include("Atoms.jl")
include("utils.jl")
include("PGBF.jl")
include("CGBF.jl")
include("overlap.jl")
include("kinetic.jl")
include("nuclear.jl")
include("coulomb.jl")
include("sto3g.jl")
include("BasisSet.jl")

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
    s = PGBF(1.0)
    px = PGBF( 1.0, (0.0,0.0,0.0),  (1,0,0) )
    @assert coulomb(s,s,s,px)==0 # 0
    @assert isapprox(coulomb(s,s,px,px), 0.9403159725793305 )
    @printf("Pass test_coul1\n")
end

test_two_terms()

test_coul1()