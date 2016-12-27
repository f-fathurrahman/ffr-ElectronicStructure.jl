function gen_rho( dr, σ1, σ2 )
  Npoints = size(dr)[1]
  rho = Array( Float64, Npoints )
  c1 = 2*σ1^2
  c2 = 2*σ2^2
  cc1 = sqrt(2*pi*σ1^2)^3
  cc2 = sqrt(2*pi*σ2^2)^3
  for ip=1:Npoints
    g1 = exp(-dr[ip]^2/c1)/cc1
    g2 = exp(-dr[ip]^2/c2)/cc2
    rho[ip] = g2 - g1
  end
  return rho
end
