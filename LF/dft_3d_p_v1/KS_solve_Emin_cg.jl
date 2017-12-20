function KS_solve_Emin_cg( LF::LF3dGrid, Gv::GvectorsT,
                          V_ionic, Focc, Nstates::Int64;
                          v0 = nothing,
                          Potentials0 = nothing,
                          α_t=3e-5, Niter=1000,
                          verbose=false )
    #
    Npoints = LF.Nx * LF.Ny * LF.Nz
    ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
    #
    # Setup trial solution
    #
    if v0 == nothing
        # random trial solution
        srand(1234)
        v = rand( Npoints, Nstates )
        v = orthonormalize( LF, v )
    else
        v = copy(v0)
        # v0 is assumed to be ortonormalized
    end
    #
    # Setup potentials
    #
    if Potentials0 == nothing
        # initialize potential: Ionic, Hartree, XC
        Potentials = PotentialsT( V_ionic, zeros(Npoints), zeros(Npoints) )
    else
        Potentials = copy(Potentials0)
    end

    # Initialize energy components: Total, Kinetic, Ionic, Hartree, XC
    Energies = EnergiesT( 0.0, 0.0, 0.0, 0.0, 0.0 )
    #
    # Calculate energies and update potetials
    #
    calc_Energies!( LF, Gv, Energies, Potentials, Focc, v )
    #
    Etot_old = Energies.Total
    print_Energies( Energies )
    #
    g     = zeros( Float64, Npoints, Nstates )
    g_old = zeros( Float64, Npoints, Nstates )
    d     = zeros( Float64, Npoints, Nstates )
    d_old = zeros( Float64, Npoints, Nstates )
    #
    β = 0
    α = 0
    for iter = 1:Niter
        #
        g = calc_grad( LF, Potentials, Focc, v )
        #
        if iter != 1
            β = trace( g' * g)/trace( g_old'*g_old )
            #β = trace( (g-g_old)'*g ) / trace( g_old'*g_old )
            #β = trace( (g-g_old)'*g )/ trace( (g-g_old)'*d )
        end
        d = -g + β*d_old
        #
        # compute gradient at trial step
        #
        v2 = orthonormalize(LF, v + α_t*d)
        rho = calc_rho( Focc, v2 )
        V_Hartree = solve_poisson_FFT( Gv, rho )
        V_xc = excVWN( rho ) + rho .* excpVWN( rho )
        wPot = PotentialsT( Potentials.Ionic, V_Hartree, V_xc )
        g_t = calc_grad( LF, wPot, Focc, v2 )
        #
        # compute estimate of best step and update current trial vectors
        #
        denum = trace( (g - g_t)' * d )
        if denum != 0.0
            α = abs( α_t * trace( g'*d )/denum )
        else
            α = 0.0
        end
        #
        v = v + α*d
        #
        # Calculate energies and update potetials
        #
        v = orthonormalize( LF, v )
        calc_Energies!( LF, Gv, Energies, Potentials, Focc, v )
        Etot = Energies.Total
        #
        if verbose
            @printf("%8d %16.9f %16.9e\n", iter, Etot, abs(Etot-Etot_old))
        end
        if abs(Etot-Etot_old) < 1.e-7
            if verbose
                @printf("Emin CG converges in %8d iterations\n", iter)
            end
            break
        end
        Etot_old = Etot
        #
        g_old = g[:,:]
        d_old = d[:,:]
    end
    return Energies, v, Potentials
end
