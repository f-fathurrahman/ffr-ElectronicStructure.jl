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
    rrloc = r/rloc
    Vpot[ip] = -Zion/r * erf( r/(sqrt(2.0)*rloc) ) +
               (C1 + C2*rrloc^2)*exp(-0.5*(rrloc)^2)
  end
  return Vpot
end
