include("smear_FD.jl")
include("calc_Focc.jl")
include("calc_entropy.jl")

function test_main(kT::Float64)
    Nstates = 8
    Nelectrons = 6
    evals = Array{Float64}(Nstates)
    evals = [-2.4, -1.0, -0.5, -0.2, -0.19, -0.10, -0.05, -0.05]
    println("\nkT = ", kT)
    spinpol = false
    println("\nspinpol = ", spinpol)
    Focc, E_fermi = calc_Focc(evals, Nelectrons, kT)
    for ist = 1:Nstates
        @printf("%18.10f %18.10f\n", evals[ist], Focc[ist])
    end
    @printf("E_fermi = %18.10f\n", E_fermi)
    integFocc = 0.0
    for ist = 1:Nstates
        if evals[ist] <= E_fermi 
            integFocc = integFocc + Focc[ist]
        end
    end
    @printf("integFocc = %18.10f\n", integFocc)
    @printf("sum(Focc) = %18.10f\n", sum(Focc))
    @printf("Entropy (-TS) = %18.10f\n", calc_entropy(Focc, kT, is_spinpol=spinpol))
end

test_main(0.01)
