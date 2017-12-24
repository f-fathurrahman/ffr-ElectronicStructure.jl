function getocc(evals, Nocc, Tbeta)

    Nstates = length(evals)
    const TOL = 1e-15
    const MAXITER = 100

    occ = zeros(Nstates)

    #
    # use bisection to find efermi such that 
    #       sum_i fermidirac(evals(i)) = Nocc
    if Nstates > Nocc 
        ilb = Nocc - 1
        iub = Nocc + 1
        lb = evals[ilb]
        ub = evals[iub]
        # make sure flb < Nocc and fub > Nocc
        #@printf("lb = %f, ub = %f\n", lb, ub)
        Focc_lb = fermidirac(evals,lb,Tbeta)
        Focc_ub = fermidirac(evals,ub,Tbeta)
        #println("Focc_lb = ", Focc_lb)
        #println("Focc_ub = ", Focc_ub)
        flb = sum(Focc_lb)
        fub = sum(Focc_ub)
        while ( (Nocc-flb)*(fub-Nocc) < 0 )
            #@printf("getocc: initial bounds are off:\n");
            #@printf("flb = %11.3e, fub = %11.3e, Nocc = %d\n", flb,fub,Nocc)
            if (flb > Nocc)
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
            if (fub < Nocc)
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
        
        efermi = (lb+ub)/2
        occ = fermidirac(evals,efermi,Tbeta)
        occsum = sum(occ)
        
        iter = 1
        
        while ( abs(occsum-Nocc) > TOL && iter < MAXITER )
            #@printf("iter = %d, efermi = %11.3e, sum = %11.3e\n", iter, efermi, occsum)
            #@printf("lb = %11.3e, ub = %11.3e\n", lb, ub)
            if (occsum < Nocc)
                lb = efermi
            else
                ub = efermi
            end
            efermi = (lb + ub)/2
            occ = fermidirac(evals,efermi,Tbeta)
            occsum = sum(occ)
            iter = iter + 1
        end #
    
    # 
    elseif (Nstates == Nocc)
        @printf("Nstates is equal to Nocc\n")
        occ    = ones(Nocc)
        efermi = evals[Nocc]
    
    else
        @printf("ERROR: The number of eigenvalues in evals should be larger than Nocc")
    
    end

    return occ, efermi

end
