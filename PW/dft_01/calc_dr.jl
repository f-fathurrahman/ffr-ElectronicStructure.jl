function calc_dr( r, center )
  Npoints = size(r)[2]
  dr = Array(Float64,Npoints)
  #
  for ip=1:Npoints
    dx2 = ( r[1,ip] - center[1] )^2
    dy2 = ( r[2,ip] - center[2] )^2
    dz2 = ( r[3,ip] - center[3] )^2
    dr[ip] = sqrt( dx2 + dy2 + dz2 )
  end
  return dr
end
