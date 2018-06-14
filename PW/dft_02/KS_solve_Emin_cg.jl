function KS_solve_Emin_cg( pw::PWGrid, V_ionic, Focc, Nstates::Int;
                           psi0=nothing, Potentials0 = nothing, E_NN = 0.0,
                           α_t = 3e-5, NiterMax=1000, verbose=false )

    Ns = pw.Ns
    Npoints = prod(Ns)
    Ngwx = pw.gvecw.Ngwx

    if psi0 == nothing
        srand(2222)
        psi = rand(ComplexF64,Ngwx,Nstates)
        psi = ortho_gram_schmidt(psi)
    else
        psi = copy(psi0)
    end

    if Potentials0 == nothing
        Potentials = PotentialsT( V_ionic, zeros(Npoints), zeros(Npoints) )
        rho = calc_rho( pw, Focc, psi )
        Potentials.Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
        Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )
    else
        Potentials = PotentialsT( Potentials0.Ionic,
                                  Potentials0.Hartree,
                                  Potentials0.XC )
    end

    d = zeros(ComplexF64, Ngwx, Nstates)
    g_old = zeros(ComplexF64, Ngwx, Nstates)
    d_old = zeros(ComplexF64, Ngwx, Nstates)
    Kg = zeros(ComplexF64, Ngwx, Nstates)
    Kg_old = zeros(ComplexF64, Ngwx, Nstates)

    β        = 0.0
    Etot_old = 0.0
    Etot     = 0.0
    if( E_NN != nothing )
        Energies = EnergiesT( 0.0, 0.0, 0.0, 0.0, 0.0, E_NN )
    else
        Energies = EnergiesT( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
    end

    CONVERGED = 0

    for iter = 1:NiterMax

        g = calc_grad( pw, Potentials, Focc, psi)
        Kg = Kprec(pw,g)

        if iter != 1
            # Fletcher-Reeves
            #β = real(sum(conj(g).*Kg))/real(sum(conj(g_old).*Kg_old))

            # Polak-Ribiere
            β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g_old).*Kg_old))
            
            # Hestenes-Stiefel
            #β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g-g_old).*d))

            # Dai-Yuan
            #β = real(sum(conj(g).*Kg))/real(sum((g-g_old).*conj(d_old)))

            # CG descent
            #β = -real(sum(conj(g).*Kg))/real(sum(conj(d_old).*Kg_old))

            # Liu-Storey
            #β = -real(sum(conj(g-g_old).*Kg))/real(sum(conj(d_old).*Kg_old))

            # Muhamed-Rivaie-Mustafa
            #gknorm = real( sum(conj(g).*g) )
            #gk1norm = real( sum(conj(g_old).*g_old) )
            #gknorm = real( sum(conj(g).*Kg) )
            #gk1norm = real( sum(conj(g_old).*Kg_old) )
            #
            #gg = g .- gknorm/gk1norm*g_old
            #β = real(sum( conj(g).*gg )) / ( gk1norm^2 + real( sum(conj(g).*d_old) ) )
            #β = real(sum( conj(gg).*Kg )) / ( gk1norm^2 + real( sum(conj(d_old).*Kg) ) )
        end
        
        if β < 0.0
            @printf("β is smaller than 0, setting it to zero\n")
            β = 0.0
        end

        d = -Kprec(pw, g) + β * d_old

        psic = ortho_gram_schmidt(psi + α_t*d)
        rho = calc_rho( pw, Focc, psic )
        Potentials.Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
        Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )
        gt = calc_grad( pw, Potentials, Focc, psic )

        denum = real(sum(conj(g-gt).*d))
        if denum != 0.0
            α = abs( α_t*real(sum(conj(g).*d))/denum )
        else
            α = 0.0
        end

        # Update wavefunction
        psi = psi[:,:] + α*d[:,:]

        # Update potentials
        psi = ortho_gram_schmidt(psi)
        rho = calc_rho( pw, Focc, psi )

        Potentials.Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
        Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )

        Energies = calc_Energies( pw, Potentials, Focc, psi, Energies.NN )
        Etot = Energies.Total

        diff = abs(Etot-Etot_old)
        @printf("CG step %8d = %18.10f %10.7e\n", iter, Etot, diff)
        
        if diff < 1e-6
            CONVERGED = CONVERGED + 1
        else
            CONVERGED = 0
        end

        if CONVERGED >= 2
            @printf("CONVERGENCE ACHIEVED\n")
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Etot_old = Etot
    end
    return psi, Energies, Potentials
    #
end
