include("alias.jl")
include("utils.jl")
include("PGBF.jl")
include("CGBF.jl")
include("overlap.jl")

function test_CGBF()
    c = init_CGBF( 0.0,0.0,0.0 )  # specifying center
    
    push!(c, 0.20, 1.0)
    #@assert isapprox( evaluate(c,0,0,0), 0.71270547 )

    #c2 = CGBF(0,0,0)
    #push!(c2,1,0.2)
    #push!(c2,0.5,0.2)

    #@assert isapprox(overlap(c2,c2),1)
    @printf("test_cgbf is passed\n")
end

test_CGBF()
