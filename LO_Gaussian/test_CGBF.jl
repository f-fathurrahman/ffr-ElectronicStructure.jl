if VERSION > v"0.6.3"
    using Printf
end

include("alias.jl")
include("utils.jl")
include("PGBF.jl")
include("CGBF.jl")
include("overlap.jl")

function test_CGBF()
    c1 = init_CGBF( 0.0,0.0,0.0 )  # specifying center
    push!(c1, 1.0, 1.0)
    @printf("Test CGBF evaluate: %f\n", abs(evaluate(c1, 0.0, 0.0, 0.0)-0.71270547))

    c2 = init_CGBF( 0.0, 0.0, 0.0 )
    push!( c2, 1.0, 0.2 )
    push!( c2, 0.5, 0.2 )
    @printf("Test CGBF overlap, <c2|c2>: %f\n", abs(overlap(c2,c2)-1.0))

    @printf("Test CGBF overlap, <c1|c2>: %f\n", overlap(c1,c2))
end

test_CGBF()
