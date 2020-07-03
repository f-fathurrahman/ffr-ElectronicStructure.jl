# !INPUT/OUTPUT PARAMETERS:
#   m : order of derivative (in,integer)
#   n : number of points (in,integer)
#   x : abscissa array (in,real(n))
#   f : function array (in,real(n))
#   g : (anti-)derivative of f (out,real(n))
# !DESCRIPTION:
#   Given function $f$ defined on a set of points $x_i$ then if $m\ge 0$ this
#   routine computes the $m$th derivative of $f$ at each point. If $m=-1$ the
#   anti-derivative of $f$ given by
#   $$ g(x_i)=\int_{x_1}^{x_i} f(x)\,dx $$
#   is calculated. Both derivatives and integrals are computed by first fitting
#   the function to a clamped cubic spline.
function fderiv!(m, n, x, f, g)

    @assert n > 0

    # automatic arrays
    cf = zeros(Float64,3,n)

    # high accuracy integration/differentiation from spline interpolation
    spline!(n, x, f, cf)
    if m == -1
        ss =0.0
        g[1] = 0.0
        for i in 1:n-1
          dx = x[i+1] - x[i]
          ss = ss + dx*( f[i] + dx*( 0.50*cf[1,i] + dx*( cf[2,i]/3 + dx*0.25*cf[3,i] ) ) )
          g[i+1] = ss
        end
    elseif m == 1
        for i in 1:n
            g[i] = cf[1,i]
        end
    elseif m == 2
        for i in 1:n
            g[i] = 2.0*cf[2,i]
        end
    elseif m == 3
        for i in 1:n
            g[i] = 6.0*cf[3,i]
        end
    else
        for i in 1:n
            g[i] = 0.0
        end
    end
    return
end
