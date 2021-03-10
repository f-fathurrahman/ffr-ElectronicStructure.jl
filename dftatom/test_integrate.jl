using Printf

include("gen_rmesh_exp.jl")
include("gen_drmesh_exp.jl")
include("integrate.jl")

function my_func(r)
    return exp(-0.2*r)
    #return exp(-0.2*r*r)
end

function main()
    r_min = 0.0
    r_max = 200.0
    a = 1e9
    N = 500
    
    rmesh = gen_rmesh_exp(r_min, r_max, a, N)
    drmesh = gen_drmesh_exp(r_min, r_max, a, N)

    f = my_func.(rmesh)
    s = integrate_trapz_7(f, drmesh)

    println("s = ", s)

end

main()