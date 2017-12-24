function KS_solve_scf( pw::PWGrid,
                       V_ionic, Focc, Ncols::Int64;
                       β = 0.5, E_NN = nothing,
                       rho0 = nothing,
                       Potentials0 = nothing,
                       NiterMax=100,
                       verbose=false )

    Ngwx = pw.gvecw.Ngwx
    Ns = pw.Ns
    Npoints = prod(Ns)
    ΔV = pw.Ω/Npoints

    # starting rho
    if rho0 == nothing
        srand(1234)
        v = rand( Complex128, Ngwx, Ncols )
        v = ortho_gram_schmidt( v )
        rho = calc_rho( pw, Focc, v )
    else
        srand(1234)
        v = rand( Complex128, Ngwx, Ncols ) # What if v0 is provided ?
        rho = copy( rho )
    end

    # starting Potentials
    if Potentials0 == nothing
        # initialize potential: Ionic, Hartree, XC
        V_Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
        V_xc = excVWN( rho ) + rho .* excpVWN( rho )
        Potentials = PotentialsT( V_ionic, V_Hartree, V_xc )
    else
        Potentials = PotentialsT( Potentials0.Ionic,
                                  Potentials0.Hartree,
                                  Potentials0.XC )
    end

    Etot_old = 0.0
    rho_old = copy(rho)

    if( E_NN != nothing )
        Energies = EnergiesT( 0.0, 0.0, 0.0, 0.0, 0.0, E_NN )
    else
        Energies = EnergiesT( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
    end

    λ = zeros(Ncols)

    for iter = 1:NiterMax

        λ, v = diag_lobpcg( pw, Potentials, v, verbose_last=false )

        # Calculate energies and update potetials
        Energies = calc_Energies( pw, Potentials, Focc, v, Energies.NN )

        Etot = Energies.Total

        diffE = abs( Etot - Etot_old )
        Etot_old = Etot

        if diffE < 1e-6
            @printf("SCF is is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
            break
        end
        #
        rho_new = calc_rho( pw, Focc, v )
        diffRho = norm(rho_new - rho)
        #
        @printf("SCF: %8d %18.10f %18.10e %18.10e\n", iter, Etot, diffE, diffRho )
        #
        rho = β*rho_new[:] + (1-β)*rho[:]
        #
        #integRho = sum(rho)*ΔV
        #@printf("integRho = %18.10f\n", integRho)
        #
        V_Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
        V_xc = excVWN( rho ) + rho .* excpVWN( rho )
        Potentials = PotentialsT( V_ionic, V_Hartree, V_xc )
    end

    return Energies, Potentials, v, λ

end
