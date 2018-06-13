function KS_solve_Emin_CG_smearing( pw::PWGrid, V_ionic, Focc, Nstates::Int;
                           psi0=nothing, Potentials0 = nothing, E_NN = 0.0,
                           α_t = 3e-5, NiterMax=1000, verbose=false, kT=0.01 )

    Ns = pw.Ns
    Npoints = prod(Ns)
    Ngwx = pw.gvecw.Ngwx

    if psi0 == nothing
        srand(2222)
        psi = randn(Ngwx,Nstates) + im*randn(Ngwx,Nstates)
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

    d_Haux = zeros(ComplexF64, Nstates, Nstates)
    g_Haux_old = zeros(ComplexF64, Nstates, Nstates)
    d_Haux_old = zeros(ComplexF64, Nstates, Nstates)
    Kg_Haux = zeros(ComplexF64, Nstates, Nstates)
    Kg_Haux_old = zeros(ComplexF64, Nstates, Nstates)

    β        = 0.0
    β_Haux   = 0.0
    Etot_old = 0.0
    Etot     = 0.0
    if( E_NN != nothing )
        Energies = EnergiesT( 0.0, 0.0, 0.0, 0.0, 0.0, E_NN )
    else
        Energies = EnergiesT( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
    end

    Nelectrons = sum(Focc)
    @printf("Nelectrons = %18.10f\n", Nelectrons)

    # Get eigenvalues
    mu = psi' * op_H(pw, Potentials, psi)
    evals = sort(real(eigvals(mu)))
    #
    kT = 0.01
    Focc, E_fermi = calc_Focc(evals, Nelectrons, kT, is_spinpol=false)
    for ist = 1:Nstates
        @printf("%3d %18.10f %18.10f\n", ist, evals[ist], Focc[ist])
    end
    @printf("E_fermi = %18.10f\n", E_fermi)

    const κ = 0.1

    Haux = zeros(ComplexF64,Nstates,Nstates)
    for ist = 1:Nstates
        Haux[ist,ist] = evals[ist]
    end

    for iter = 1:NiterMax

        g = calc_grad(pw, Potentials, Focc, psi)
        g_Haux = calc_grad_Haux(pw, Potentials, Focc, evals, psi, kT)

        Kg = Kprec(pw,g)
        Kg_Haux = κ*g_Haux

        if iter != 1
            #β = real(sum(conj(g).*Kg))/real(sum(conj(g_old).*Kg_old))
            β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g_old).*Kg_old))
            β_Haux = real(sum(conj(g_Haux-g_Haux_old).*Kg_Haux))/real(sum(conj(g_Haux_old).*Kg_Haux_old))
            #β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g-g_old).*d))
            #β = real(sum(conj(g).*Kg))/real(sum((g-g_old).*conj(d_old)))
        end
        if β < 0.0
            @printf("β is smaller than 0, setting it to zero\n")
            β = 0.0
        end
        if β_Haux < 0.0
            @printf("β_Haux is smaller than 0, setting it to zero\n")
            β_Haux = 0.0
        end

        d = -Kg + β * d_old
        d_Haux = -Kg_Haux + β_Haux * d_Haux_old

        psic = ortho_gram_schmidt(psi + α_t*d)
        rho = calc_rho( pw, Focc, psic )
        Potentials.Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
        Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )
        gt = calc_grad( pw, Potentials, Focc, psic )

        exit()

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

        Energies = calc_Energies(pw, Potentials, Focc, psi, Energies.NN)
        Entropies = calc_entropies(Focc, kT, is_spinpol=false)
        Etot = Energies.Total + Entropies

        diff = abs(Etot-Etot_old)
        @printf("CG step %8d = %18.10f %10.7e\n", iter, Etot, diff)
        if diff < 1e-6
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
