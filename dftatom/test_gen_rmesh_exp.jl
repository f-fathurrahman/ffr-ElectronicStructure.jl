using Printf

include("gen_rmesh_exp.jl")

function main()
    r_min = 0.0
    r_max = 50.0
    a = 1e9
    N = 10
    
    rmesh = gen_rmesh_exp(r_min, r_max, a, N)

    for i in 1:N+1
        @printf("%8d %18.10e\n", i, rmesh[i])
    end
end

main()