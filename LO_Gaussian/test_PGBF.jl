include("utils.jl")
include("PGBF.jl")
include("overlap.jl")

using PyPlot

const plt = PyPlot

function test_main()

  #fun1 = init_PGBF(1.0)
  fun1 = init_PGBF(0.3, I=6)

  basis_info(fun1)

  const N = 100
  x = linspace(-5.0,5.0,N)
  y = zeros(N)
  for i = 1:N
    y[i] = amplitude( fun1, x[i], 0.0, 0.0 )
  end

  plt.clf()
  plt.plot( x, y, marker="o" )
  plt.grid()
  plt.savefig("func1_x.png", dpi=300)

end

function basis_info(f::PGBF)
  @printf("\n")
  @printf("Info for PGBF\n")
  @printf("Exponent: %18.10f\n", f.expn)
  @printf("Center: (%18.10f,%18.10f,%18.10f)\n", f.x, f.y, f.z)
  @printf("Angular momentum: (%2d,%2d,%2d)\n", f.I, f.J, f.K)
  @printf("Norm: %18.10f\n", f.norm)
  @printf("\n")
end

test_main()
