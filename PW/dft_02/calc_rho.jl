function calc_rho( pw::PWGrid, Focc, psi::Array{ComplexF64,2} )
    Ω  = pw.Ω
    Ns = pw.Ns
    Npoints = prod(Ns)
    Nstates = size(psi)[2]
    #
    rho = zeros(ComplexF64,Npoints)

    # Transform to real space
    cpsi = zeros( ComplexF64, Npoints, Nstates )
    cpsi[pw.gvecw.idx_gw2r,:] = psi[:,:]
    psiR = G_to_R(Ns, cpsi)

    # orthonormalization in real space
    ortho_gram_schmidt!( Nstates, psiR )
    psiR = sqrt(Npoints/Ω)*psiR

    for is = 1:Nstates
        for ip = 1:Npoints
            rho[ip] = rho[ip] + Focc[is]*conj(psiR[ip,is])*psiR[ip,is]
        end
    end
    return real(rho)
end
