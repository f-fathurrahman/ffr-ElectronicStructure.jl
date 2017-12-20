function init_pot_Hcoul( LF, center )
    Npoints = LF.Nx * LF.Ny * LF.Nz
    Vpot = zeros( Float64,Npoints )
    for ip = 1:Npoints
        r = norm( LF.lingrid[:,ip] - center[:] )
        Vpot[ip] = -1.0/r
    end
    return Vpot
end
