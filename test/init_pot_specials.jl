"""
Initialize harmonic external potential for electronic system.
"""
function init_pot_harm_3d( LF::LF3dGrid, ω, center )
  Npoints = LF.Nx * LF.Ny * LF.Nz
  Vpot = zeros( Float64,Npoints )
  for ip = 1:Npoints
    Vpot[ip] = 0.5 * ω^2 * norm(LF.lingrid[:,ip] - center[:])^2
  end
  return Vpot
end



"""
Initialize ion-electron potential according to pure Coulombic
potential.
This potential is singular at r = 0.
"""
function init_pot_Hcoul( LF::LF3dGrid, center )
  Npoints = LF.Nx * LF.Ny * LF.Nz
  Vpot = zeros( Float64,Npoints )
  for ip = 1:Npoints
    r = norm( LF.lingrid[:,ip] - center[:] )
    Vpot[ip] = -1.0/r
  end
  return Vpot
end


"""
Initialize ion-electron potential for hydrogen atom using
analytic form and parameters used in Gygi's paper.
(need reference)
"""
function init_pot_Hps( LF::LF3dGrid, center )
  Npoints = LF.Nx * LF.Ny * LF.Nz
  Vpot = zeros( Float64,Npoints )

  const rc1 = 0.25
  const rc2 = 0.284
  const aa  = -1.9287
  const bb  = 0.3374

  for ip=1:Npoints
    r = norm( LF.lingrid[:,ip] - center[:] )
    Vpot[ip] = -1.0/r * erf( r/rc1 ) + (aa + bb*r^2)*exp(-(r/rc2)^2)
  end
  return Vpot
end



"""
Initialize local ion-electron potential for hydrogen atom according
to HGH pseudopotential parameters.

This function is used only for testing purpose.

This function should not be used to initialize ion-electron potential for
periodic system.
"""
function init_pot_Hps_HGH( LF::LF3dGrid, center )
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
