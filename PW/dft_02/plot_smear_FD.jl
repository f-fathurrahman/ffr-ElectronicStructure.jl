include("smear_FD.jl")

using PyPlot

const plt = PyPlot

function test_plot(kT)

    mu = -2.0

    @printf("kT = %f\n", kT)

    xmin = -5.0
    xmax =  5.0

    NptsPlot = 500
    x = Array{Float64}(NptsPlot)
    f = Array{Float64}(NptsPlot)
    for i = 1:NptsPlot
        x[i] = xmin + (i-1)*(xmax-xmin)/(NptsPlot-1)
        f[i] = smear_FD(x[i], mu, kT, is_spinpol=true)
    end

    plt.plot(x, f, marker="o", label="kT="*string(kT))

end

function test_main()
    plt.clf()
    test_plot(0.001)
    test_plot(0.01)
    test_plot(0.03)
    test_plot(0.1)
    test_plot(0.2)
    test_plot(0.3)
    test_plot(0.4)
    plt.legend()
    plt.grid()
    plt.savefig("plot_smear_FD.png", dpi=300)
end

test_main()
