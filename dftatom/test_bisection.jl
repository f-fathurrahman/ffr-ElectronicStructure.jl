using Printf

include("gen_rmesh_exp.jl")
include("lagrange_interp.jl")
include("radial_interp.jl")
include("rsch_integ_rk4.jl")
include("sch_rk4_step.jl")
include("get_min_idx.jl")
include("get_n_nodes.jl")

function main()
    # Build radial function
    r_min = 1e-8
    r_max = 50.0
    a = 1e6
    Nr = 3001
    rmesh = gen_rmesh_exp(r_min, r_max, a, Nr)

    Z = 3
    n = 1
    l = 0

    P = zeros(Float64,Nr)
    Q = zeros(Float64,Nr)
    V = -Z ./ rmesh

    Emin = -5000.0
    Emax = 0.0
    #E = 0.5*(Emin + Emax)
    E = -1000.0

    for iterShoot in 1:30

        println("\niterShoot = ", iterShoot)
        println("E = ", E)
        dE = abs(Emax - Emin)
        if dE < 1e-10
            println("Converged by dE")
            println("dE = ", dE)
            break
        end

        ctp = Nr # no perturbation correction

        imax = rsch_integ_rk4!( E, Z, l, rmesh, V, P, Q )
        println("imax = ", imax)

        minidx = get_min_idx(imax, P)
        println("minidx = ", minidx)

        #@views P[minidx:Nr] .= 0.0
        #@views Q[minidx:Nr] .= 0.0

        Nnodes = get_n_nodes( imax, P )
        println("Nnodes = ", Nnodes)
        println("minidx = ", minidx)


        cond1 = Nnodes != n - l - 1
        cond2 = ctp == Nr
        cond3 = imax < ctp

        if cond1 || cond2 || cond3
            is_big = Nnodes > (n - l - 1)
            println("n = ", n, " l = ", l)
            println("is_big = ", is_big)
            @printf("before: %18.10f %18.10f\n", Emin, Emax)
            if is_big
                Emax = E
            else
                Emin = E
            end
            @printf("After: %18.10f %18.10f\n", Emin, Emax)
            #
            E = 0.5*(Emin + Emax)
            println("E = ", E)
        end

    end

end

@time main()
