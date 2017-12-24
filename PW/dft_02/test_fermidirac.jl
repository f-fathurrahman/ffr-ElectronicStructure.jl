include("fermidirac.jl")

function test_fermidirac()
    Nstates = 4
    ev = Array{Float64}(Nstates)
    Tbeta = 0.1
    ev = [-2.4, -1.0, -0.5, -0.2]
    efermi = ev[end]
    Focc = fermidirac(ev, efermi, Tbeta)
    for ist = 1:Nstates
        @printf("%f %f\n", ev[ist], Focc[ist])
    end
    @printf("sum(Focc) = %f\n", sum(Focc))
end

test_fermidirac()