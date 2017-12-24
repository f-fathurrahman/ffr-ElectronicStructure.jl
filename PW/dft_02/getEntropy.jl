function getEntropy( Focc, Tbeta )
    const SMALL = 1.e-10
    Nstates = length(Focc)
    e = 0.0
    for ist = 1:Nstates
        if Focc[ist] > SMALL
            e = e + Focc[ist]*log(Focc[ist])
        else
            e = e + (2.0 - Focc[ist])*log(2.0 - Focc[ist]) # spin-degenerate case
        end
    end
    return e/Tbeta
end
