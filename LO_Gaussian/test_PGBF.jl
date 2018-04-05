include("alias.jl")
include("utils.jl")
include("PGBF.jl")
include("overlap.jl")

using PyPlot

const plt = PyPlot

function do_plot_x(bf::PGBF)
    const N = 100
    x = linspace(-5.0,5.0,N)
    y = zeros(N)
    for i = 1:N
      y[i] = evaluate( bf, x[i], 0.0, 0.0 )
    end
    plt.clf()
    plt.plot( x, y, marker="o" )
    plt.grid()
    plt.savefig("func1_x.png", dpi=300)
end

function test_main()
  bf1 = PGBF(0.3, (0.0,0.0,0.0), (6,0,0))
  println(bf1)
  do_plot_x(bf1)
end

@time test_main()
