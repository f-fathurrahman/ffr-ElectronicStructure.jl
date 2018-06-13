function KS_solve_SCF_smearing( pw::PWGrid,
                       V_ionic, Focc_in, Nstates::Int64, Nelectrons::Int64;
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
        v = rand( ComplexF64, Ngwx, Nstates )
        v = ortho_gram_schmidt( v )
        rho = calc_rho( pw, Focc_in, v )
    else
        srand(1234)
        v = rand( ComplexF64, Ngwx, Nstates ) # What if v0 is provided ?
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

    λ = zeros(Nstates)

    const kT = 0.01

    MIXDIM = 4
    df = zeros(Float64,Npoints,MIXDIM)
    dv = zeros(Float64,Npoints,MIXDIM)

    for iter = 1:150

        λ, v = diag_lobpcg( pw, Potentials, v, verbose_last=false )

        Focc, E_fermi = calc_Focc(λ, Nelectrons, kT)
        Entropy = calc_entropy(Focc, kT)

        for ist = 1:Nstates
            @printf("%5d %18.10f %18.10f\n", ist, λ[ist], Focc[ist] )
        end
        
        println("sum(Focc) = ", sum(Focc))
        println("integFocc = ", sum_upto_E_fermi(Focc, λ, E_fermi))
        println("E_fermi = ", E_fermi)

        rho_new = calc_rho( pw, Focc, v )
        diffRho = norm(rho_new - rho)
        
        #@printf("Using simple mix with β = %10.5f\n", β)
        #rho = β*rho_new[:] + (1-β)*rho[:]

        @printf("Using andersonmix with β = %f, MIXDIM = %d\n", β, MIXDIM)
        rho = andersonmix!( rho, rho_new, β, df, dv, iter, MIXDIM )
        
        for ip = 1:Npoints
            if rho[ip] < 1e-12
                rho[ip] = 1e-12
            end
        end
        #
        integRho = sum(rho)*ΔV
        #@printf("integRho = %18.10f\n", integRho)
        diffNelec = abs(integRho - Nelectrons)
        if diffNelec > 1e-12
            println("WARNING: Rescaling electron density, diffNelec = ", diffNelec)
            rho[:] = rho[:]*Nelectrons/integRho
        end

        V_Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
        V_xc = excVWN( rho ) + rho .* excpVWN( rho )
        Potentials = PotentialsT( V_ionic, V_Hartree, V_xc )

        # Calculate energies and update potetials
        Energies = calc_Energies( pw, Potentials, Focc, v, Energies.NN )

        Etot = Energies.Total + Entropy
        diffE = abs( Etot - Etot_old )

        #
        @printf("SCF: %8d %18.10f %18.10e %18.10e\n", iter, Etot, diffE, diffRho )        

        if diffE < 1e-7
            @printf("SCF is is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
            @printf("Etot           = %18.10f\n", Etot)
            @printf("Entropy        = %18.10f\n", Entropy)
            @printf("Energies.Total = %18.10f\n", Etot)        
            break
        end
        #
        Etot_old = Etot
    end

    exit()

    return Energies, Potentials, v, λ

end
