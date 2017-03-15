function get_rho( v::Array{Float64,2} )
  Ncols   = size(v)[2]
  Npoints = size(v)[1]
  Rho = zeros(Float64,Npoints)
  for ic = 1:Ncols
    Rho[:] = Rho[:] + v[:,ic] .* v[:,ic]
  end
  return Rho
end
