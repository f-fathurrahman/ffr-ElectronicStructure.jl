"""
Solve system of linear equations: L*x = b, where L is Laplacian matrix
using conjugate residual method
"""
function linsolve_cr!( LF::LF3dGrid, x::Array{Float64,1}, b::Array{Float64,1};
                          NiterMax=1000, verbose=false, TOL=5.e-10 )
  #
  Npoints = size(x)[1]
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
    if verbose
      @printf("%8d %18.10f %18.10f %18.10f\n", iter, norm_res, norm_p, norm_Δx)
    end
    #
    if norm_Δx < TOL
      if verbose
        @printf("#CONV in linsolve_cr: iter, norm_res: %8d %10.5e\n", iter, norm_res)
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
  return
  #
end # of function
