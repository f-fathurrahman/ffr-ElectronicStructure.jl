"""
Solve system of linear equations: L*x = b, where L is Laplacian matrix
"""
function linsolve_cg_v2!( LF::LF3dGrid, x::Array{Float64,1}, b::Array{Float64,1};
                          NiterMax=1000, verbose=false, TOL=5.e-10 )
  #
  Npoints = size(x)[1]
  #
  L_x = apply_Laplacian( LF, x )
  #
  r = b[:] - L_x[:]
  z = r[:]
  p = z[:]

  r_old = r[:]
  z_old = z[:]
  p_old = p[:]
  for iter = 1 : NiterMax
    #
    L_x = apply_Laplacian( LF, p )
    #
    α = dot( r, z ) / dot( p, L_x )
    #
    x[:] = x[:] + α * p[:]  # FIXME use x[:] to force x to be copied, not referenced
    r[:] = r[:] - α * L_x[:]
    #
    norm_res = sqrt( dot( r, r ) )
    if verbose
      @printf("%8d %20.10f\n", iter, norm_res)
    end
    #
    if norm_res < TOL
      if verbose
        @printf("#CONV in linsolve_cg_v2: iter, norm_res: %8d %10.5e\n", iter, norm_res)
      end
      break
    end
    #
    z = r[:]  #XXX call appropriate preconditioner here
    # Polak-Ribiere
    β = dot( z, r - r_old ) / dot( z_old, r_old )
    p = z[:] + β * p_old[:]
    #
    z_old = z[:]
    r_old = r[:]
    p_old = p[:]
  end
  #
  return
  #
end # of function
