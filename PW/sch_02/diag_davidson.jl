function diag_davidson( pw::PWGrid, Vpot, X0;
                        tol=1e-5, tol_avg=1e-7, NiterMax=200, verbose=false )

    # get size info
    Ncols = size(X0)[2]
    Ngwx  = size(X0)[1]

    if Ncols <= 0
        @printf("diag_davidson requires at least one initial wave function!\n");
        return
    end
    # orthonormalize the initial wave functions.
    X = ortho_gram_schmidt(X0)  # normalize (again)?

    evals    = zeros(Float64, Ncols)
    R        = zeros(ComplexF64, Ngwx, Ncols)
    Hred     = zeros(ComplexF64, 2*Ncols, 2*Ncols)
    Sred     = zeros(ComplexF64, 2*Ncols, 2*Ncols)
    res      = zeros(Float64, Ncols)
    res_norm = zeros(Float64, Ncols)

    HX = op_H( pw, Vpot, X )

    # Initial eigenvalues
    for ic = 1:Ncols
        for ig = 1:Ngwx
            evals[ic] = evals[ic] + real( conj(X[ig,ic]) * HX[ig,ic] )
        end
    end

    # Calculate residuals
    for ic = 1:Ncols
        for ig = 1:Ngwx
            R[ig,ic] = evals[ic]*X[ig,ic] - HX[ig,ic]
        end
    end

    for ic = 1:Ncols
        res[ic] = 0.0
        for ig = 1:Ngwx
            res[ic] = res[ic] + real( R[ig,ic] * conj(R[ig,ic]) )
        end
        res[ic] = sqrt( res[ic] )
    end

    EPS = eps()

    set1 = 1:Ncols
    set2 = Ncols+1:2*Ncols

    sum_evals = sum(evals)
    sum_evals_old = sum_evals

    for iter = 1:NiterMax

        res_norm[:] .= 1.0

        for ic = 1:Ncols
            if EPS < res[ic]
                res_norm[ic] = 1.0/res[ic]
            end
        end

        for ic = 1:Ncols
            for ig = 1:Ngwx
                R[ig,ic] = res_norm[ic] * R[ig,ic]
            end
        end

        R = Kprec(pw, R)

        HR = op_H( pw, Vpot, R )

        # FIXME: Pull this outside the loop?
        if iter == 1
            Hred[set1,set1] = X' * HX
        else
            Hred[1:Ncols,1:Ncols] .= 0.0 + im*0.0
            for ic = 1:Ncols
                Hred[ic,ic] = evals[ic] + im*0.0 # use diagm ?
            end
            #Hred[set1,set1] = diagm(evals[:])
        end

        Hred[set1,set2] = X' * HR
        Hred[set2,set2] = R' * HR
        Hred[set2,set1] = Hred[set1,set2]'

        Sred[set1,set1] = Matrix(Diagonal(ones(ComplexF64,Ncols)))
        Sred[set1,set2] = X' * R
        Sred[set2,set2] = R' * R
        Sred[set2,set1] = Sred[set1,set2]'

        Hred = (Hred + Hred')/2.0
        Sred = (Sred + Sred')/2.0

        λ_red, X_red = eigen( Hermitian(Hred), Hermitian(Sred) )

        evals = λ_red[1:Ncols]

        #for ic = 1:2*Ncols
        #  print(ic); print(" "); println(λ_red[ic])
        #end

        X  = X  * X_red[set1,set1] + R  * X_red[set2,set1]
        HX = HX * X_red[set1,set1] + HR * X_red[set2,set1]

        # Calculate residuals
        for ic = 1:Ncols
            for ig = 1:Ngwx
                R[ig,ic] = evals[ic]*X[ig,ic] - HX[ig,ic]
            end
        end

        for ic = 1:Ncols
            res[ic] = 0.0
            for ig = 1:Ngwx
                res[ic] = res[ic] + real( R[ig,ic] * conj(R[ig,ic]) )
            end
            res[ic] = sqrt( res[ic] )
            @printf("%4d %18.10f %18.10f\n", ic, evals[ic], res[ic] )
        end
        
        sum_evals = sum(evals)
        diffSum = abs(sum_evals - sum_evals_old)
        @printf("iter %d tol_avg = %18.10e\n", iter, abs(diffSum))
        if diffSum < tol_avg
            if verbose
                @printf("Davidson convergence: tol_avg\n")
            end
            break
        end
        sum_evals_old = sum_evals
    end

    return evals, X

end
