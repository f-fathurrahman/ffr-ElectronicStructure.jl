# Probably this should not be used in the periodic system
function init_pot_Hps_HGH_G( Gvec::GvectorsT, center )

    Ns = Gvec.Ns
    Ngvec = prod(Ns)
    Npoints = Ngvec
    Ω = Gvec.Ω
    G2 = Gvec.G2

    Vpot = zeros( Float64,Npoints )

    const zval = 1
    const rloc = 0.2
    const c1 = -4.0663326
    const c2 = 0.6678322

    Vg = zeros(Complex128,Npoints)
    pre1 = -4*pi*zval/Ω
    pre2 = sqrt(8*pi^3)*rloc^3/Ω
    #
    for ig=2:Ngvec
        Gr = sqrt(G2[ig])*rloc
        expGr2 = exp(-0.5*Gr^2)
        Vg[ig] = pre1/G2[ig]*expGr2 + pre2*expGr2 * (c1 + c2*(3-Gr^2))
    end
    # limiting value, with minus sign ?
    Vg[1] = 2*pi*zval*rloc^2 + (2*pi)^1.5 * rloc^3 * (c1 + 3.0*c2)

    V_ionic = real( G_to_R(Ns,Vg) ) * Npoints

    return V_ionic
end


function init_pot_Hps( LF, center )
    Npoints = LF.Nx * LF.Ny * LF.Nz
    Vpot = zeros( Float64,Npoints )

    const rc1 = 0.25
    const rc2 = 0.284
    const aa  = -1.9287
    const bb  = 0.3374

    # TODO Add journal reference for this pseudopotential
    for ip=1:Npoints
        r = norm( LF.lingrid[:,ip] - center[:] )
        Vpot[ip] = -1.0/r * erf( r/rc1 ) + (aa + bb*r^2)*exp(-(r/rc2)^2)
    end
    return Vpot
end


# Probably this should not be used in the periodic system
function init_pot_Hps_HGH( LF, center )
    Npoints = LF.Nx * LF.Ny * LF.Nz
    Vpot = zeros( Float64,Npoints )

    const Zion = 1
    const rloc = 0.2
    const C1 = -4.0663326
    const C2 = 0.6678322

    # TODO Add journal reference
    for ip = 1:Npoints
        r = norm( LF.lingrid[:,ip] - center[:] )
        if r < 1e-7
            println("WARNING: small r:", r);
            println("Using limiting value:")
            Vpot[ip] = -2*Zion/(sqrt(2*pi)*rloc) + C1
        else
            rrloc = r/rloc
            Vpot[ip] = -Zion/r * erf( r/(sqrt(2.0)*rloc) ) +
                     (C1 + C2*rrloc^2)*exp(-0.5*(rrloc)^2)
        end
    end
    return Vpot
end
