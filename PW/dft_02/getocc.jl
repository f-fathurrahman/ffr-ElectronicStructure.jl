function getocc(ev, nocc, Tbeta)

    Nstates = length(ev)
    const TOL = 1e-15
    const MAXITER = 100

    occ = zeros(Nstates)

    #
    # use bisection to find efermi such that 
    #       sum_i fermidirac(ev(i)) = nocc
    if Nstates > nocc 
        ilb = nocc - 1
        iub = nocc + 1
        lb = ev[ilb]
        ub = ev[iub]
        # make sure flb < nocc and fub > nocc
        #@printf("lb = %f, ub = %f\n", lb, ub)
        Focc_lb = fermidirac(ev,lb,Tbeta)
        Focc_ub = fermidirac(ev,ub,Tbeta)
        #println("Focc_lb = ", Focc_lb)
        #println("Focc_ub = ", Focc_ub)
        flb = sum(Focc_lb)
        fub = sum(Focc_ub)
        while ( (nocc-flb)*(fub-nocc) < 0 )
            #@printf("getocc: initial bounds are off:\n");
            #@printf("flb = %11.3e, fub = %11.3e, nocc = %d\n", flb,fub,nocc)
            if (flb > nocc)
                if (ilb > 1)
                    ilb = ilb - 1
                    lb = ev[ilb]
                    flb = sum( fermidirac(ev,lb,Tbeta) )
                else
                    @printf("getocc: cannot find a lower bound for efermi, something is wrong\n")
                    exit()
                end
            end
            #
            if (fub < nocc)
                if (iub < Nstates)
                    iub = iub + 1
                    ub  = ev[iub]
                    fub = sum(fermidirac(ev,ub,Tbeta))
                else
                    @printf("getocc: cannot find an upper bound for efermi\n")
                    @printf("something is wrong, try increasing the number of wavefunctions in X0\n")
                    exit()
                end
            end
        end  # while
        @printf("flb = %11.3e, fub = %11.3e\n", flb, fub)
        
        efermi = (lb+ub)/2
        occ = fermidirac(ev,efermi,Tbeta)
        occsum = sum(occ)
        
        iter = 1
        
        while ( abs(occsum-nocc) > TOL && iter < MAXITER )
            #@printf("iter = %d, efermi = %11.3e, sum = %11.3e\n", iter, efermi, occsum)
            #@printf("lb = %11.3e, ub = %11.3e\n", lb, ub)
            if (occsum < nocc)
                lb = efermi
            else
                ub = efermi
            end
            efermi = (lb + ub)/2
            occ = fermidirac(ev,efermi,Tbeta)
            occsum = sum(occ)
            iter = iter + 1
        end #
    
    # 
    elseif (Nstates == nocc)
        @printf("Nstates is equal to nocc\n")
        occ    = ones(nocc)
        efermi = ev[nocc]
    
    else
        @printf("ERROR: The number of eigenvalues in ev should be larger than nocc")
    
    end

    return occ, efermi

end
