using Printf
import PyPlot

const plt = PyPlot

include("gen_rmesh_exp.jl")
include("lagrange_interp.jl")
include("radial_interp.jl")

function myfunc(x)
    return exp(-x)
end

function main()
    r_min = 0.0
    r_max = 10.0
    a = 1e9
    N = 50
    
    rmesh = gen_rmesh_exp(r_min, r_max, a, N)
    fr = myfunc.(rmesh)

    #NptsPlot = 5*N
    #rmesh_dense = gen_rmesh_exp(r_min, r_max, a, NptsPlot)
    #fr_interp = zeros(NptsPlot)
    #for i in 1:NptsPlot
    #    fr_interp[i] = radial_interp(rmesh, fr, )
    #end

    Vmid = zeros(N-1)
    rmid = zeros(N-1)
    for i in 1:N-1
        rmid[i] = 0.5*(rmesh[i] + rmesh[i+1])
        Vmid[i] = radial_interp( rmesh, fr, rmid[i], i+1)
        Vtrue = myfunc(rmid[i])
        err = abs(Vmid[i] - Vtrue)
        @printf("%18.10f %18.10f %18.10e\n", Vmid[i], Vtrue, err)
    end

    plt.clf()
    plt.plot(rmesh, fr, marker="o", label="data")
    plt.plot(rmid, Vmid, marker="o", label="midval")
    plt.xlim(0.0, 2.0)
    plt.legend()
    plt.savefig("IMG_test_radial_interp.png")
end

main()