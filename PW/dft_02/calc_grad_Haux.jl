function calc_grad_Haux(pw::PWGrid, Potentials, Focc, evals, psi::Array{ComplexF64,2}, kT)
    Nstates = size(Focc)[1]
    g = zeros(ComplexF64,Nstates,Nstates)
    #
    denum_dmudeta = 0.0
    for ist = 1:Nstates
        denum_dmudeta = denum_dmudeta + Focc[ist]*(1 - Focc[ist])
    end
    #println("denum_dmudeta =", denum_dmudeta)
    #
    # Diagonal elements
    Hpsi = op_H(pw, Potentials, psi)
    dFdmu = 0.0 + im*0.0
    #
    for ist = 1:Nstates
        term1 = (dot( psi[:,ist], Hpsi[:,ist] ) - evals[ist])/kT*Focc[ist]*(1 - Focc[ist])
        #println("term1 = ", term1)
        g[ist,ist] = -term1
        dFdmu = dFdmu + term1
    end
    for ist = 1:Nstates
        g[ist,ist] = g[ist,ist] + Focc[ist]*(1 - Focc[ist])/denum_dmudeta * dFdmu
    end
    #
    # Off diagonals
    for ist = 1:Nstates
        for jst = (ist+1):Nstates
            pref = (Focc[jst] - Focc[ist])/(evals[jst] - evals[ist])
            #@printf("pref = %2d %2d %18.10f\n", ist, jst, pref)
            g[ist,jst] = pref*dot( psi[:,ist], Hpsi[:,jst] )
            g[jst,ist] = conj(g[ist,jst])
        end
    end
    #println("Pass 36 sum g Haux = ", sum(g))
    return g
end
