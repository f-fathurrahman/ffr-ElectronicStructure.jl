function init_pot_harm_3d( LF, ω, center )
  Npoints = LF.Nx * LF.Ny * LF.Nz
  Vpot = zeros( Float64,Npoints )
  for ip = 1:Npoints
    Vpot[ip] = 0.5 * ω^2 * norm(LF.lingrid[:,ip] - center[:])^2
  end
  return Vpot
end
