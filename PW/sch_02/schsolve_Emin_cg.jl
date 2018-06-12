function schsolve_Emin_cg( pw::PWGrid, Vpot, psi::Array{ComplexF64,2};
                           NiterMax=1000 )
    #
    Npoints = size(psi)[1]
    Nstates = size(psi)[2]

    g = zeros(ComplexF64, Npoints, Nstates)
    Kg = zeros(ComplexF64, Npoints, Nstates)
    Kg_old = zeros(ComplexF64, Npoints, Nstates)
    psic = zeros(ComplexF64, Npoints, Nstates )

    d = zeros(ComplexF64, Npoints, Nstates)
    g_old = zeros(ComplexF64, Npoints, Nstates)
    d_old = zeros(ComplexF64, Npoints, Nstates)
    gt = zeros(ComplexF64, Npoints, Nstates)

    yk1 = zeros(ComplexF64, Npoints, Nstates) # for testing Hager-Zhang beta
    v1 = zeros(ComplexF64, Npoints, Nstates)

    alphat = 1.e-5
    beta = 0.0
    Etot_old = 0.0
    Etot = 0.0
    #
    for iter = 1:NiterMax
        
        g = calc_grad( pw, Vpot,  psi)
        Kg = Kprec(pw,g)

        if iter != 1
          
            #beta = real(tr(g'*Kg))/real(tr(g_old'*Kg_old))
            #beta = real(tr((g-g_old)'*Kprec(pw,g)))/real(tr(g_old'*Kprec(pw,g_old)))
            #beta = real(tr((g-g_old)'*Kg))/real(tr((g-g_old)'*d))
          
            yk1 = g - g_old
            tt = real(tr(d_old'*yk1))
            v1 = yk1 - 2*d_old*norm(yk1)/tt
            beta = real(tr(v1'*g)/tt)
        end
        #println("beta HZ = ", beta)

        if beta < 0.0
            println("beta is negative, resetting beta to zero")
            beta = 0.0
        end

        d = -Kprec(pw, g) + beta * d_old
        psic = ortho_gram_schmidt(psi + alphat*d)
        gt = calc_grad( pw, Vpot, psic )

        denum = real(tr((g-gt)'*d))
        if denum != 0.0
          alpha = abs(alphat*real(tr(g'*d))/denum )
        else
          alpha = 0.0
        end

        psi = psi[:,:] + alpha*d[:,:]

        psi = ortho_gram_schmidt(psi)
        Etot = calc_Etot( pw, Vpot, psi )

        diff = abs(Etot-Etot_old)
        @printf("E step %8d = %18.10f %18.10e\n", iter, Etot, diff)
        if diff < 1e-6
            @printf("CONVERGENCE ACHIEVED\n")
            break
        end
        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Etot_old = Etot
    end
    return psi, Etot
end
