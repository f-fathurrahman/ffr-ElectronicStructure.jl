if VERSION > v"0.6.2"
    using Printf
end

include("alias.jl")
include("utils.jl")
include("PGBF.jl")
include("CGBF.jl")
include("overlap.jl")

if VERSION <= v"0.6.2"
    using PyPlot
    const plt = PyPlot
end

function do_plot_x( bf::PGBF, filesave::String )
    const N = 100
    x = linspace(-5.0,5.0,N)
    y = zeros(N)
    for i = 1:N
      y[i] = evaluate( bf, x[i], 0.0, 0.0 )
    end
    plt.clf()
    plt.plot( x, y, marker="o" )
    plt.grid()
    plt.savefig(filesave, dpi=300)
end

function test_main()
    bf1 = PGBF(0.3, (0.0,0.0,0.0), (1,0,0))
    println(bf1)  
    #do_plot_x(bf1,"bf1.png")

    bf2 = PGBF(0.4,(0.1,0.0,0.0), (1,0,0))
    println(bf2)
    #do_plot_x(bf2,"bf2.png")

    S1 = overlap(bf1,bf1)
    @printf("overlap(bf1,bf1) = %18.10f\n", S1)

    S1 = overlap(bf2,bf2)
    @printf("overlap(bf2,bf2) = %18.10f\n", S1)

    S1 = overlap(bf1,bf2)
    @printf("overlap(bf1,bf2) = %18.10f\n", S1)
end

@time test_main()
