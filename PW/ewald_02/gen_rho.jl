function gen_rho( Ns, dr, sigma, Sf )
    Npoints = size(dr)[1]
    g1 = Array{Float64}(Npoints)
    c1 = 2*sigma^2
    cc1 = sqrt(2*pi*sigma^2)^3
    for ip=1:Npoints
        g1[ip] = exp(-dr[ip]^2/c1)/cc1
    end
    ctmp = R_to_G(Ns,g1)
    for ip=1:Npoints
        ctmp[ip] = ctmp[ip]*Sf[ip]
    end
    return real( G_to_R(Ns,ctmp) )
end
