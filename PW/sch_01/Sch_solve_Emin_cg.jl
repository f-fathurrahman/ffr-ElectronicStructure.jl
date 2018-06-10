include("linmin1.jl")
include("linmin2.jl")

function Sch_solve_Emin_cg( pw::PWGrid, Vpot, psi::Array{ComplexF64,2};
                            NiterMax=1000 )
    #
    Npoints = size(psi)[1]
    Nstates = size(psi)[2]
    d = zeros(ComplexF64, Npoints, Nstates)
    g_old  = zeros(ComplexF64, Npoints, Nstates)
    d_old  = zeros(ComplexF64, Npoints, Nstates)
    Kg     = zeros(ComplexF64, Npoints, Nstates)
    Kg_old = zeros(ComplexF64, Npoints, Nstates)
    #
    α_t = 3e-5
    β = 0.0

    Etot = calc_Etot( pw, Vpot, psi )
    Etot_old = Etot

    #
    for iter = 1:NiterMax
        g = calc_grad( pw, Vpot,  psi)
        nrm = 0.0
        for is = 1:Nstates
          nrm = nrm + real( dot( g[:,is], g[:,is] ) )
        end

        Kg = Kprec(pw,g)
        #Kg = copy(g)  # test for no preconditioner

        if iter != 1
          #β = real(sum(conj(g).*Kg))/real(sum(conj(g_old).*Kg_old))
          #β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g_old).*Kg_old))
          #β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g-g_old).*d_old))
          β = real(sum(conj(g).*Kg))/real(sum((g-g_old).*conj(d_old)))
        end

        if β < 0.0
            println("β is negative: set to zero")
            β = 0.0
        end

        d = -Kg + β * d_old

        #α = linmin2( pw, Vpot, α_t, psi, d, g )

        α = linmin1( pw, Vpot, α_t, psi, d, Etot_old, g )

        # Update wavefunction
        psi = psi[:,:] + α*d[:,:]

        psi = ortho_gram_schmidt(psi)
        Etot = calc_Etot( pw, Vpot, psi )

        diff = abs(Etot-Etot_old)
        @printf("E step %8d = %18.10f %18.10f %18.10f\n", iter, Etot, diff, nrm/Nstates)
        if diff < 1e-6
            @printf("CONVERGENCE ACHIEVED\n")
            break
        end
        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Etot_old = Etot
    end
    return psi, Etot
    #
end
