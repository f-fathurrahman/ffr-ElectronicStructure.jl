# spin-unpolarized version
function getocc(evals, Nelectrons, Tbeta)

    Nstates = length(evals)
    const TOL = 1e-15
    const MAXITER = 100

    Focc = zeros(Nstates)
    Nocc = round(Int64,Nelectrons/2)

    #
    # use bisection to find efermi such that 
    #       sum_i fermidirac(evals(i)) = Nelectrons
    if Nstates > Nocc
        ilb = Nocc - 1
        iub = Nocc + 1
        lb = evals[ilb]
        ub = evals[iub]
        # make sure flb < Nelectrons and fub > Nelectrons
        #@printf("lb = %f, ub = %f\n", lb, ub)
        Focc_lb = fermidirac(evals,lb,Tbeta)
        Focc_ub = fermidirac(evals,ub,Tbeta)
        #println("Focc_lb = ", Focc_lb)
        #println("Focc_ub = ", Focc_ub)
        flb = sum(Focc_lb)
        fub = sum(Focc_ub)
        while ( (Nelectrons-flb)*(fub-Nelectrons) < 0 )
            #@printf("getocc: initial bounds are off:\n");
            #@printf("flb = %11.3e, fub = %11.3e, Nelectrons = %d\n", flb,fub,Nelectrons)
            if (flb > Nelectrons)
                if (ilb > 1)
                    ilb = ilb - 1
                    lb = evals[ilb]
                    flb = sum( fermidirac(evals,lb,Tbeta) )
                else
                    @printf("getocc: cannot find a lower bound for efermi, something is wrong\n")
                    exit()
                end
            end
            #
            if (fub < Nelectrons)
                if (iub < Nstates)
                    iub = iub + 1
                    ub  = evals[iub]
                    fub = sum(fermidirac(evals,ub,Tbeta))
                else
                    @printf("getocc: cannot find an upper bound for efermi\n")
                    @printf("something is wrong, try increasing the number of wavefunctions in X0\n")
                    exit()
                end
            end
        end  # while
        
        #@printf("flb = %11.3e, fub = %11.3e\n", flb, fub)
        
        efermi = (lb + ub)/2
        Focc = fermidirac(evals,efermi,Tbeta)
        occsum = sum(Focc)
        
        iter = 1
        
        while ( abs(occsum-Nelectrons) > TOL && iter < MAXITER )
            #@printf("iter = %d, efermi = %11.3e, sum = %11.3e\n", iter, efermi, occsum)
            #@printf("lb = %11.3e, ub = %11.3e\n", lb, ub)
            if (occsum < Nelectrons)
                lb = efermi
            else
                ub = efermi
            end
            efermi = (lb + ub)/2
            Focc = fermidirac(evals,efermi,Tbeta)
            occsum = sum(Focc)
            iter = iter + 1
        end #
    
    # 
    elseif (Nstates == Nocc)
        @printf("Nstates is equal to Nocc\n")
        Focc    = ones(Nelectrons)
        efermi = evals[Nelectrons]
    
    else
        @printf("ERROR: The number of eigenvalues in evals should be larger than Nelectrons")
    
    end

    return Focc, efermi

end
