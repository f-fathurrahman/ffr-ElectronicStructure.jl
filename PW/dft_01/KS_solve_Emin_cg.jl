function KS_solve_Emin_cg( pw::PWGrid, V_ionic, Focc, Nstates::Int;
                           psi0=nothing, Potentials0 = nothing,
                           α_t = 3e-5, NiterMax=1000, verbose=false )

    Ns = pw.Ns
    Npoints = prod(Ns)

    if psi0 == nothing
        srand(2222)
        psi = rand(ComplexF64,Npoints,Nstates)
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

    d      = zeros(ComplexF64, Npoints, Nstates)
    g_old  = zeros(ComplexF64, Npoints, Nstates)
    d_old  = zeros(ComplexF64, Npoints, Nstates)
    Kg     = zeros(ComplexF64, Npoints, Nstates)
    Kg_old = zeros(ComplexF64, Npoints, Nstates)

    yk1 = zeros(ComplexF64, Npoints, Nstates) # for testing Hager-Zhang beta
    v1 = zeros(ComplexF64, Npoints, Nstates)

    β     = 0.0
    Etot_old = 0.0
    Etot     = 0.0
    Energies = EnergiesT( 0.0, 0.0, 0.0, 0.0, 0.0 )

    CONVERGED = 0

    for iter = 1:NiterMax

        g = calc_grad( pw, Potentials, Focc, psi)
        Kg = Kprec(pw,g)
        if iter != 1
            #β = real(sum(conj(g).*Kg))/real(sum(conj(g_old).*Kg_old))
            #β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g_old).*Kg_old))
            #β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g-g_old).*d))
            β = real(sum(conj(g).*Kg))/real(sum((g-g_old).*conj(d_old)))

            #yk1 = g - g_old
            #tt = real(tr(d_old'*yk1))
            #v1 = yk1 - 2*d_old*norm(yk1)/tt
            #β = real(tr(v1'*Kg)/tt)

            #@printf("\nβ = %f\n", β)
        end

        if β < 0.0
            println("Encountered β less than 0 setting it to zero")
            β = 0.0
        end
        #
        d = -Kg + β * d_old
        #
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

        Energies = calc_Energies( pw, Potentials, Focc, psi )
        Etot = Energies.Total

        diff = abs(Etot-Etot_old)
        @printf("Emin CG: %8d = %18.10f %18.10f\n", iter, Etot, diff)

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
