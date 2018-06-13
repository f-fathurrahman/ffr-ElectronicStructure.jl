"""
An implementation of simple method to calculate Ewald energy
"""
function calc_ewald( pw::PWGrid, Xpos, Sf; sigma=0.25 )

    @printf("\n***************************************************************\n")
    @printf("Calculating Ewald energy.\n")
    @printf("WARNING: This function should be called only for hydrogen atom.\n")
    @printf("***************************************************************\n")
    #
    Ω  = pw.Ω
    r  = pw.r
    Ns = pw.Ns
    G2 = pw.gvec.G2
    Npoints = prod(Ns)
    #
    # Generate array of distances
    center = sum(pw.LatVecs,dims=2)/2
    dr = gen_dr( r, center )
    #
    # Generate charge density
    rho = gen_rho( Ns, dr, sigma, Sf )
    intrho = sum(rho)*Ω/Npoints
    @printf("sigma, int_rho: %10.5f %18.10f\n", sigma, intrho)
    #
    # Solve Poisson equation and calculate Hartree energy
    ctmp = 4.0*pi*R_to_G( Ns, rho )
    ctmp[1] = 0.0
    for ip = 2:Npoints
        ctmp[ip] = ctmp[ip] / G2[ip]
    end
    phi = real( G_to_R( Ns, ctmp ) )
    Ehartree = 0.5*dot( phi, rho ) * Ω/Npoints
    #
    Eself = 1.0/(2*sqrt(pi))*(1.0/sigma)*size(Xpos,2)
    E_nn = Ehartree - Eself
    #@printf("Ehartree, Eself, E_nn = %20.16f %20.16f %20.16f\n", Ehartree, Eself, E_nn)
    return E_nn
end


function gen_rho( Ns, dr, sigma, Sf )
    Npoints = size(dr)[1]
    g1 = Array{Float64}(undef,Npoints)
    c1 = 2*sigma^2
    cc1 = sqrt(2*pi*sigma^2)^3
    for ip=1:Npoints
        g1[ip] = exp(-dr[ip]^2/c1)/cc1
    end
    ctmp = R_to_G(Ns,g1)
    for ip=1:Npoints
        ctmp[ip] = ctmp[ip]*Sf[ip]
    end
    return real(G_to_R(Ns,ctmp))
end
