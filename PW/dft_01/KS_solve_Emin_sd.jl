function KS_solve_Emin_sd( pw::PWGrid, V_ionic, Focc, Nstates::Int;
                           psi0=nothing, Potentials0 = nothing,
                           α_t = 3e-5, NiterMax=1000, verbose=false )
    Ns = pw.Ns
    Npoints = prod(Ns)

    if psi0 == nothing
        srand(2222)
        psi = randn(Npoints,Nstates) + im*randn(Npoints,Nstates)
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

    Etot = 0.0
    Etot_old = 0.0
    Energies = EnergiesT( 0.0, 0.0, 0.0, 0.0, 0.0 )

    for iter = 1:NiterMax
        #
        psi = psi - α_t*calc_grad( pw, Potentials, Focc, psi )

        psi = ortho_gram_schmidt(psi)

        # Update potentials
        rho = calc_rho( pw, Focc, psi)
        Potentials.Hartree = real( G_to_R( Ns, Poisson_solve( pw, rho ) ) )
        Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )

        Energies = calc_Energies( pw, Potentials, Focc, psi )
        Etot = Energies.Total

        conv = abs(Etot-Etot_old)
        @printf("Emin SD: (s, step, conv) %8d %18.10f %18.10e\n", iter, Etot, conv )
        if conv < 1e-6
            print("Convergence achieved\n")
            break
        end
        Etot_old = Etot
    end
    return psi, Energies, Potentials
end
