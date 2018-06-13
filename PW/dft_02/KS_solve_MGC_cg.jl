function KS_solve_MGC_cg( pw::PWGrid, V_ionic, Focc, Nstates::Int;
                           psi0=nothing, Potentials0 = nothing,
                           α_t = 3e-5, NiterMax=1000, verbose=false )

    Ns = pw.Ns
    Npoints = prod(Ns)
    Ngwx = pw.gvecw.Ngwx
    Ω = pw.Ω

    if psi0 == nothing
        srand(2222)
        psi = randn(Ngwx,Nstates) + im*randn(Ngwx,Nstates)
        psi = ortho_gram_schmidt(psi)
    else
        psi = copy(psi0)
    end

    Potentials = PotentialsT( V_ionic, zeros(Npoints), zeros(Npoints) )

    d = zeros(ComplexF64, Ngwx, Nstates)
    g_old = zeros(ComplexF64, Ngwx, Nstates)
    d_old = zeros(ComplexF64, Ngwx, Nstates)
    Kg = zeros(ComplexF64, Ngwx, Nstates)
    Kg_old = zeros(ComplexF64, Ngwx, Nstates)

    β        = 0.0
    Etot     = calc_E_MGC( pw, Potentials, psi )
    Etot_old = Etot

    ctmp = zeros(ComplexF64, Npoints, Nstates)
    idx = pw.gvecw.idx_gw2r

    S = zeros(ComplexF64,Nstates,Nstates)
    Q = zeros(ComplexF64,Nstates,Nstates)

    rho = zeros(Float64,Npoints)

    for iter = 1:NiterMax

        g = calc_grad_MGC( pw, Potentials, Focc, psi)
        Kg = Kprec(pw,g)

        println("sum(g) = ", sum(g))
        println("sum(Kg) = ", sum(Kg))

        if iter != 1
            #β = real(sum(conj(g).*Kg))/real(sum(conj(g_old).*Kg_old))
            β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g_old).*Kg_old))
            #β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g-g_old).*d))
            #β = real(sum(conj(g).*Kg))/real(sum((g-g_old).*conj(d_old)))
        end

        d = -Kprec(pw, g) + β * d_old

        Δ = 2.0*real( sum( conj(d).*g ) )
        psi_t = ortho_gram_schmidt( psi + α_t*d )
        E_trial = calc_E_MGC( pw, Potentials, psi_t )
        curvature = ( E_trial - ( Etot + α_t*Δ ) ) /α_t^2
        α = -Δ/(2*curvature)

        println("β = ", β)
        println("α = ", α)

        # Update wavefunction
        psi = psi[:,:] + α*d[:,:]

        Etot = calc_E_MGC( pw, Potentials, psi )

        diff = abs(Etot-Etot_old)
        @printf("CG step %8d = %18.10f %18.10f\n", iter, Etot, diff)
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


function calc_E_MGC( pw::PWGrid, Potentials, ϕ::Array{ComplexF64,2}; Nel=1.0, η=1.0 )

    Ngwx    = size(ϕ)[1]
    Nstates = size(ϕ)[2]
    Npoints = prod(pw.Ns)
    Ω = pw.Ω

    S = zeros(ComplexF64,Nstates,Nstates)
    Q = zeros(ComplexF64,Nstates,Nstates)

    for j = 1:Nstates
    for i = 1:Nstates
        S[i,j] = sum( conj(ϕ[:,i]).*ϕ[:,j] )
        if i == j
            Q[i,j] = 2 - S[i,j]
        end
    end
    end

    E_kin = 0.0
    for j = 1:Nstates
    for i = 1:Nstates
        E_kin = E_kin + real(Q[i,j]*sum(conj(ϕ[:,j] .* op_K(pw, ϕ[:,i]))))
    end
    end

    ctmp = zeros(ComplexF64, Npoints, Nstates)
    idx = pw.gvecw.idx_gw2r
    for ic = 1:Nstates
        ctmp[idx,ic] = ϕ[:,ic]
    end
    ϕ_r = G_to_R( pw.Ns, ctmp )
    # orthonormalization in real space
    ortho_gram_schmidt!( Nstates, ϕ_r )
    scale!( sqrt(Npoints/Ω), ϕ_r )

    ρ = zeros(Float64,Npoints)

    for j = 1:Nstates
    for i = 1:Nstates
        for ip = 1:Npoints
            # XXX need to use conj??? Use factor 2 ??
            ρ[ip] = real( Q[i,j] * conj(ϕ_r[ip,i]) * ϕ_r[ip,j] )
        end
    end
    end

    V_Hartree = real( G_to_R( pw.Ns, Poisson_solve(pw, ρ) ) )
    E_Hartree = 0.5*dot( V_Hartree, ρ ) * Ω/Npoints

    V_XC = excVWN( ρ ) + ρ .* excpVWN( ρ )

    E_xc = dot( excVWN(ρ), ρ ) * Ω/Npoints

    E_Ionic = dot( Potentials.Ionic, ρ ) * Ω/Npoints

    intRho = sum(ρ)*Ω/Npoints
    #println("intRho = ", intRho)

    E_MGC = E_kin + E_Hartree + E_Ionic + E_xc + η*(Nel - intRho)

    #println("E_MGC Kin     = ", E_kin)
    #println("E_MGC Ionic   = ", E_Ionic)
    #println("E_MGC Hartree = ", E_Hartree)
    #println("E_MGC XC      = ", E_xc)
    #println("E_MGC Total   = ", E_MGC)

    return E_MGC

end


function calc_grad_MGC( pw::PWGrid, Potentials, Focc, psi::Array{ComplexF64,2} )

    # Focc are assumed to be 2s
    Ngwx    = size(psi)[1]
    Nstates = size(psi)[2]

    g = zeros(ComplexF64,Ngwx,Nstates)
    H = zeros(ComplexF64,Nstates,Nstates)
    S = zeros(ComplexF64,Nstates,Nstates)
    Q = zeros(ComplexF64,Nstates,Nstates)

    Hpsi = op_H( pw, Potentials, psi )

    for m = 1:Nstates
    for n = 1:Nstates
        H[m,n] = sum( conj(psi[:,m]).*Hpsi[:,n] )
        S[m,n] = sum( conj(psi[:,m]).*psi[:,n] )
        if m == n
            Q[m,n] = 2 - S[m,m]
        end
    end
    end
    #
    η = 1.0  # ????
    #
    for n = 1:Nstates
        for m = 1:Nstates
            g[:,n] = Hpsi[:,m]*Q[m,n] - psi[:,m]*Q[m,n]*η - psi[:,m]*(H[m,n] - η*S[m,n])
        end
        g[:,n] = 2*g[:,n]
    end
    #
    return g
end
