function broydenmix1!(vin, vout, beta, df, dv, cdf, iter, mixdim, brank)

    # [vnew,df,dv,cdf] =

    if (brank > mixdim)
        @printf("brank = %d, mixdim = %d\n", brank, mixdim)
        @printf("The rank of the Broyden update must be larger than the mixing memory size")
        exit()
    end
    Npoints = size(vin)[1]
    #
    # Note that the sign convention used in this function is different
    # from that used in other mixing functions!!!
    #
    # function evaluation overwrites vout
    #
    vout = vin - vout
    bragg = zeros(Float64,Npoints)
    vupd  = beta*vout

    if (iter > 1)
        # compute the changes in function evaluations and steps (changes 
        # in potentials)
        if iter < mixdim + 1
            dv[:,iter-1] = -dv[:,iter-1] + vin
        else
            dv[:,mixdim] = -dv[:,mixdim] + vin
        end
        #
        # the difference between the current and the previous function value
        #
        df = -df + vout
        #
        # Construct C*df
        #
        yvec = beta*df
        if (iter > brank+1)
            @printf("Rank-%d update...\n",brank);
            for jbr = 1:iter-brank-1
                jsel = jbr:jbr+brank-1;
                dvmcdf  = dv[:,jsel] - cdf[:,jbr];
                G = dv[:,jsel]'*cdf[:,jbr]
                yvec = yvec + dvmcdf*(G\(dv[:,jsel]'*yvec))
            end
        
            if (brank > 1)
                cdf0 = cdf[:,jbr]
                cdf0 = cdf0[:,2:brank]
                cdf0 = cdf0 + dvmcdf*(G\(dv[:,jsel]'*cdf0))
                #cdf[:,jbr+1] = reshape([cdf0 yvec],n123*brank,1) ???
            else
                cdf[:,jbr+1] = yvec
            end
            #
            # Compute the potential update
            #
            for jbr = 1:iter-brank
                jsel = jbr:jbr+brank-1
                dvmcdf  = dv[:,jsel] - cdf[:,jbr]
                G = dv[:,jsel]'*cdf[:,jbr]
                vupd = vupd + dvmcdf*(G\(dv[:,jsel]'*vupd));
            end
        else
            #
            # perform rank-1 update until there are enough potential vectors to mix
            #
            @printf("Rank-1 update...\n");
            #
            # Save all initial df vectors in the first column of cdf
            #
            for jbr = 1:iter-2
                ibeg = (jbr-1)*Npoints + 1
                iend = jbr*Npoints
                dvmcdf = dv[:,jbr] - cdf[ibeg:iend,1]
                G = dv[:,jbr]'*cdf[ibeg:iend,1]
                yvec = yvec + dvmcdf*(G\(dv[:,jbr]'*yvec))
            end
            ibeg = Npoints*(iter-2)+1;
            iend = Npoints*(iter-1);
            cdf[ibeg:iend,1] = yvec;
            for jbr = 1:iter-1
                ibeg = (jbr-1)*Npoints + 1
                iend = jbr*Npoints
                dvmcdf = dv[:,jbr] - cdf[ibeg:iend,1]
                G = dv[:,jbr]'*cdf[ibeg:iend,1]
                vupd = vupd + dvmcdf*(G\(dv[:,jbr]'*vupd))
            end
        end
    end

    if iter < mixdim+1
        dv[:,iter] = vin
    else
        # delete the first column of dv
        #dv(:,1)=[]; 
        # append a new column at the end;
        dv[:,mixdim] = vin
    end
    df = vout
    vnew = vin - vupd

    return vnew

end