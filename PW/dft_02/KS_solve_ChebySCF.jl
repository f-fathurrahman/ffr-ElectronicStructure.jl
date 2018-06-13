function KS_solve_ChebySCF( pw::PWGrid,
                            V_ionic, Focc, Nstates::Int64;
                            cheby_degree=5,
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
        rho = calc_rho( pw, Focc, v )
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

    # eigenvalues
    λ = zeros(Nstates)
    resnrm = zeros(Nstates)

    for iter = 1:NiterMax

        ub, lb = get_ub_lb_lanczos(pw, Potentials, Nstates*2)
        #@printf("lb = %f, ub = %f\n", lb, ub)

        v = chebyfilt(pw, Potentials, v, cheby_degree, lb, ub)
        v = ortho_gram_schmidt(v)

        #
        rho_new = calc_rho( pw, Focc, v )
        diffRho = norm(rho_new - rho)
        #
        rho = β*rho_new[:] + (1-β)*rho[:]
        #
        #integRho = sum(rho)*ΔV
        #@printf("integRho = %18.10f\n", integRho)
        #
        V_Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
        V_xc = excVWN( rho ) + rho .* excpVWN( rho )
        Potentials = PotentialsT( V_ionic, V_Hartree, V_xc )

        # Calculate energies and update potetials
        Energies = calc_Energies( pw, Potentials, Focc, v, Energies.NN )

        Etot = Energies.Total
        diffE = abs( Etot - Etot_old )

        #
        @printf("ChebySCF: %8d %18.10f %18.10e %18.10e\n", iter, Etot, diffE, diffRho )        

        if diffE < 1e-7
            @printf("ChebySCF is is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
            break
        end
        #
        Etot_old = Etot
    end

    # Calculate eigenvalues
    Hv = op_H(pw,Potentials,v)
    G = v' * Hv
    R = Hv - v*G
    λ = real(eigvals(G))
    #for j = 1:Nstates
    #   resnrm[j] = norm(R[:,j]);
    #   @printf("λ[%2d] = %11.3e, resnrm = %11.3e\n", j, λ[j], resnrm[j]);
    #end

    return Energies, Potentials, v, λ

end
