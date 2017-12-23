function KS_solve_scf( LF::LF3dGrid, Gv::GvectorsT,
                      ∇2, precH,
                      V_ionic, Focc, Nstates::Int64;
                      rho0 = nothing,
                      Potentials0 = nothing,
                      α_t=3e-5, Niter=1000,
                      verbose=false )
    #
    Npoints = LF.Nx * LF.Ny * LF.Nz
    ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h

    # starting rho
    if rho0 == nothing
        srand(1234)
        v = rand( Npoints, Nstates )
        v = orthonormalize( LF, v )
        rho = calc_rho( Focc, v )
    else
        srand(1234)
        v = rand( Npoints, Nstates ) # for initial guess in
        rho = copy( rho )
    end

    # starting Potentials
    if Potentials0 == nothing
        # initialize potential: Ionic, Hartree, XC
        V_Hartree = solve_poisson_FFT( Gv, rho )
        V_xc = excVWN( rho ) + rho .* excpVWN( rho )
        Potentials = PotentialsT( V_ionic, V_Hartree, V_xc )
    else
        Potentials = copy(Potentials0)
    end

    Etot_old = 0.0
    rho_old = copy(rho)
    Energies = EnergiesT(0.0, 0.0, 0.0, 0.0, 0.0)
    λ = zeros(Nstates)
    #
    for iter = 1:150
        #
        λ, v = diag_lobpcg( LF, ∇2, precH, Potentials, v, verbose_last=true )
        v = v[:,:]/sqrt(ΔV)  # renormalize eigenvectors
        #
        rho_new = calc_rho( Focc, v )
        diffRho = norm(rho_new-rho)
        #
        # Linear mixing
        #
        rho = 0.7*rho_new[:] + 0.3*rho[:]
        #
        # Check for integrated Rho
        #
        integRho = sum(rho)*ΔV
        @printf("integRho = %18.10f\n", integRho)        
        #
        # Calculate energies and update potetials
        #
        V_Hartree = solve_poisson_FFT( Gv, rho )
        V_xc = excVWN( rho ) + rho .* excpVWN( rho )
        Potentials = PotentialsT( V_ionic, V_Hartree, V_xc )
        #        
        Energies = calc_Energies( LF, ∇2, Potentials, rho, Focc, v )
        Etot = Energies.Total
        #
        diffE = abs(Etot-Etot_old)
        @printf("%8d %18.10f %18.10e %18.10e\n", iter, Etot, diffE, diffRho )
        #        
        if diffE < 1e-7
            break
        end
        #
        Etot_old = Etot        
    end

    return Energies, λ, v, Potentials
end
