function fermidirac(ev::Array{Float64,1}, efermi::Float64, Tbeta::Float64)
    Nstates = length(ev)
    f = zeros(Nstates)
    for ist = 1:length(ev)
        #f[ist] = 1.0/( 1.0 + exp( Tbeta*(ev[ist] - efermi)) )
        f[ist] = 2.0/( 1.0 + exp( Tbeta*(ev[ist] - efermi)) ) # spin-degenerate case
    end
    return f
end

function fermidirac(ev::Float64, efermi::Float64, Tbeta::Float64)
    return 1.0/( 1.0 + exp( Tbeta*(ev - efermi)) )
end

