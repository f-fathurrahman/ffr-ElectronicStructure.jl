import PyPlot
const plt = PyPlot

include("lagrange_interp.jl")

function myfunc(x)
    return log(x) + 1
end

function main()
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    y = myfunc.(x)

    NptsPlot = 100
    xo = collect(range(x[begin], stop=x[end], length=NptsPlot))
    yo = zeros(NptsPlot)
    for i in 1:NptsPlot
        @views yo[i] = lagrange_interp(x[1:4], y[1:4], xo[i])
    end

    plt.clf()
    plt.plot(x, y, marker="o", label="data")
    plt.plot(xo, yo, label="interp")
    plt.grid(true)
    plt.legend()
    plt.savefig("IMG_test_lagrange_interp.png", dpi=150)
end

@time main()
@time main()