function linsolve_cr( LF::LF3dGrid, b::Array{Float64,1};
                      x0 = nothing,
                      NiterMax = 1000, TOL=5.e-10,
                      convmsg=false, showprogress=false )
  #
  Npoints = size(b)[1]
  #
  if x0 == nothing
    x = zeros(Float64,Npoints)
  else
    x = copy(x0)
  end
  #
  L_x = apply_Laplacian( LF, x )
  #
  r = b[:] - L_x[:]
  p = r[:]

  r_old = r[:]
  p_old = p[:]
  #
  L_p = apply_Laplacian( LF, p )
  x_old = x[:]
  #
  for iter = 1 : NiterMax
    #
    L_r = apply_Laplacian( LF, r )
    #
    rLr = dot( r, L_r )
    α = rLr / dot( L_p, L_p )  # ???
    #
    x_old = x[:]
    x[:] = x[:] + α * p[:]  # FIXME use x[:] to force x to be copied, not referenced
    #
    r_old = r[:]
    rLr_old = rLr
    r[:] = r[:] - α * L_p[:]
    #
    norm_res = sqrt( dot( r, r ) )
    norm_p   = sqrt( dot( p, p ) )
    Δx = abs.(x-x_old)
    norm_Δx = sqrt( dot( Δx, Δx ) )
    if showprogress
      @printf("%8d %18.10e %18.10e %18.10e\n", iter, norm_res, norm_p, norm_Δx)
    end
    #
    if norm_res < TOL || norm_p < TOL || norm_Δx < TOL
      if convmsg
        @printf("#Convergence achieved in linsolve_cr: %8d iterations.\n", iter )
      end
      break
    end
    #
    L_r = apply_Laplacian( LF, r ) # calculate new L_r
    β = dot( r, L_r ) / rLr_old
    #
    p_old = p[:]
    p = r[:] + β * p[:]
    #
    L_p = L_r[:] + β * L_p[:]
  end
  #
  return x
  #
end # of function
