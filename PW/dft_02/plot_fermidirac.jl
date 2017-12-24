include("fermidirac.jl")

using PyPlot

const plt = PyPlot


function test_plot(Tbeta)

    mu = -2.0
    
    #kT = mu/100
    #Tbeta = abs(1/kT)

    @printf("Tbeta = %f\n", Tbeta)

    xmin = -5.0
    xmax =  5.0

    NptsPlot = 500
    x = Array{Float64}(NptsPlot)
    f = Array{Float64}(NptsPlot)
    for i = 1:NptsPlot
        x[i] = xmin + (i-1)*(xmax-xmin)/(NptsPlot-1)
        f[i] = fermidirac(x[i], mu, Tbeta)
    end

    plt.plot(x, f, marker="o", label="Tbeta="*string(Tbeta))

end

function test_main()
    plt.clf()
    test_plot(1.0)
    test_plot(3.0)
    test_plot(5.0)
    test_plot(8.0)
    test_plot(10.0)
    test_plot(100.0)
    plt.legend()
    plt.grid()
    plt.savefig("fermidirac.png", dpi=300)
end

test_main()
