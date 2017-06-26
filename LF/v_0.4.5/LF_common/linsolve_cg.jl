function linsolve_cg( LF::LF3dGrid, b::Array{Float64,1};
                      x0 = nothing,
                      NiterMax = 1000, TOL=5.e-10,
                      convmsg=true, showprogress=false )
  #
  Npoints = size(b)[1]
  if x0 == nothing
    x = zeros(Float64, Npoints)
  else
    x = copy(x0)
  end
  #
  r = zeros( Float64, Npoints )
  p = zeros( Float64, Npoints )
  #
  L_x = apply_Laplacian( LF, x )
  for ip = 1 : Npoints
    r[ip] = b[ip] - L_x[ip]
    p[ip] = r[ip]
  end

  rsold = dot( r, r )
  for iter = 1 : NiterMax
    #
    L_x = apply_Laplacian( LF, p )
    #
    α = rsold/dot( p, L_x )
    #
    x[:] = x[:] + α * p[:]  # FIXME use x[:] to force x to be copied, not referenced
    r[:] = r[:] - α * L_x[:]
    #
    rsnew = dot( r, r )
    # deltars = rsold - rsnew
    if showprogress
      @printf("%8d %20.10f\n", iter, sqrt(rsnew))
    end
    #
    if sqrt(rsnew) < TOL
      if convmsg
        @printf("#Convergence achieved in linsolve_cg: %8d iterations.\n", iter)
      end
      break
    end
    #
    p = r + (rsnew/rsold) * p
    #
    rsold = rsnew
  end
  #
  return x
  #
end
