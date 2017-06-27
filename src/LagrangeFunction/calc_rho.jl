function calc_rho( Focc, v::Array{Float64,2} )
  Ncols   = size(v)[2]
  Npoints = size(v)[1]
  rho = zeros(Float64,Npoints)
  for ic = 1:Ncols
    rho[:] = rho[:] + v[:,ic] .* v[:,ic] * Focc[ic]
  end
  return rho
end
