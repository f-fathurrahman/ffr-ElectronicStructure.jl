include("fermidirac.jl")

function test_fermidirac()
    Nstates = 5
    ev = Array{Float64}(Nstates)
    Tbeta = 100.0  # large Tbeta, small width
    ev = [-2.4, -1.0, -0.5, -0.2, -0.19]
    efermi = 0.5*(ev[3] + ev[4])
    Focc = fermidirac(ev, efermi, Tbeta)
    for ist = 1:Nstates
        @printf("%f %f\n", ev[ist], Focc[ist])
    end
    @printf("sum(Focc) = %f\n", sum(Focc))
end

function test_fermidirac_v2()
    Nstates = 5
    ev = Array{Float64}(Nstates)
    swidth = 0.01  # Ha (same units as energy)
    ev = [-2.4, -1.0, -0.5, -0.2, -0.19]
    efermi = 0.5*(ev[3] + ev[4])    
    Focc = fermidirac_v2(ev, efermi, swidth)
    for ist = 1:Nstates
        @printf("%f %f\n", ev[ist], Focc[ist])
    end
    @printf("sum(Focc) = %f\n", sum(Focc))
end

test_fermidirac()
test_fermidirac_v2()
