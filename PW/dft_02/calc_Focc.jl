# spin-unpolarized version
function calc_Focc(evals, Nelectrons, kT::Float64; is_spinpol=false)

    Nstates = length(evals)
    const TOL = 1e-15
    const MAXITER = 100

    Focc = zeros(Nstates)
    Nocc = round(Int64,Nelectrons/2)  # normally
    @printf("Nelectrons = %d\n", Nelectrons)
    @printf("Nocc = %d\n", Nocc)

    # use bisection to find E_fermi such that 
    #  sum_i Focc(i) = Nelectrons
    if Nstates > Nocc
        ilb = Nocc - 1
        iub = Nocc + 1
        lb = evals[ilb]
        ub = evals[iub]
        # make sure flb < Nelectrons and fub > Nelectrons
        @printf("lb = %f, ub = %f\n", lb, ub)
        Focc_lb = smear_FD(evals, lb, kT, is_spinpol=is_spinpol)
        Focc_ub = smear_FD(evals, ub, kT, is_spinpol=is_spinpol)
        #println("Focc_lb = ", Focc_lb)
        #println("Focc_ub = ", Focc_ub)
        flb = sum(Focc_lb)
        fub = sum(Focc_ub)
        while ( (Nelectrons-flb)*(fub-Nelectrons) < 0 )
            @printf("getocc: initial bounds are off:\n");
            @printf("flb = %11.3e, fub = %11.3e, Nelectrons = %d\n", flb,fub,Nelectrons)
            if (flb > Nelectrons)
                if (ilb > 1)
                    ilb = ilb - 1
                    lb = evals[ilb]
                    flb = sum( smear_FD(evals, lb, kT, is_spinpol=spinpol) )
                else
                    @printf("getocc: cannot find a lower bound for E_fermi, something is wrong\n")
                    exit()
                end
            end
            #
            if (fub < Nelectrons)
                if (iub < Nstates)
                    iub = iub + 1
                    ub  = evals[iub]
                    fub = sum( smear_FD(evals, ub, kT, is_spinpol=spinpol))
                else
                    @printf("getocc: cannot find an upper bound for E_fermi\n")
                    @printf("something is wrong, try increasing the number of wavefunctions in X0\n")
                    exit()
                end
            end
        end  # while
        
        @printf("flb = %11.3e, fub = %11.3e\n", flb, fub)
        
        E_fermi = (lb + ub)/2
        Focc = smear_FD(evals, E_fermi, kT, is_spinpol=is_spinpol)
        occsum = sum(Focc)
        
        iter = 1
        
        while ( abs(occsum-Nelectrons) > TOL && iter < MAXITER )
            @printf("iter = %d, E_fermi = %11.3e, sum = %11.3e\n", iter, E_fermi, occsum)
            @printf("lb = %11.3e, ub = %11.3e\n", lb, ub)
            if (occsum < Nelectrons)
                lb = E_fermi
            else
                ub = E_fermi
            end
            E_fermi = (lb + ub)/2
            Focc = smear_FD(evals, E_fermi, kT, is_spinpol=is_spinpol)
            occsum = sum(Focc)
            iter = iter + 1
        end #
    
    # 
    elseif (Nstates == Nocc)
        @printf("Nstates is equal to Nocc\n")
        if spinpol
            Focc = 2.0*ones(Nstates)
            E_fermi = evals[Nstates]
        else
            Focc    = ones(Nstates)
            E_fermi = evals[Nstates]
        end
    
    else
        @printf("ERROR: The number of eigenvalues in evals should be larger than Nelectrons")
        exit()
    end

    return Focc, E_fermi

end
