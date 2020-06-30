function polynm(m::Int64, np::Int64, xa, ya, x::Float64)

    res = 0.0
    
    if np == 1
        if m < 0
            res = ya[1]*(x - xa[1])
        elseif m == 0
            res = ya[1]
        else
            res = 0.0
        end
        return res
    end


    if np == 2
        c1 = ( ya[2] - ya[1] )/( xa[2] - xa[1] )
        t1 = x - xa[1]
        if m < 0
          res = t1*( ya[1] + 0.5*c1*t1 )
        elseif m == 0
          res = c1*t1 + ya[1]
        elseif m == 1
          res = c1
        else
          res = 0.0
        end
        return res
    end

    if np == 3
        x0 = xa[1]
        x1 = xa[2] - x0
        x2 = xa[3] - x0
        y0 = ya[1]
        y1 = ya[2] - y0
        y2 = ya[3] - y0
        t0 = 1.0/(x1*x2*(x2-x1))
        t1 = x1*y2
        t2 = x2*y1
        c1 = x2*t2 - x1*t1
        c2 = t1 - t2
        t1 = x - x0
        if m < 0
          res = t1*(y0 + t0*t1*( 0.5*c1 + c2*t1/3) )
        elseif m == 0
          res = y0 + t0*t1*(c1 + c2*t1)
        elseif m == 1
          res = t0*(2.0*c2*t1 + c1)
        elseif m == 2
          res = t0*2.0*c2
        else
          res = 0.0
        end
        return res
    end
    
    if np == 4
        x0 = xa[1]
        x1 = xa[2] - x0
        x2 = xa[3] - x0
        x3 = xa[4] - x0
        y0 = ya[1]
        y1 = ya[2] - y0
        y2 = ya[3] - y0
        y3 = ya[4] - y0
        t4 = x1 - x2
        t5 = x1 - x3
        t6 = x2 - x3
        t1 = x1*x2*y3
        t2 = x2*x3*y1
        t3 = x1*x3
        t0 = 1.0/(x2*t3*t4*t5*t6)
        t3 = t3*y2
        c3 = t1*t4 + t2*t6 - t3*t5
        t4 = x1^2
        t5 = x2^2
        t6 = x3^2
        y1 = t3*t6 - t1*t5
        y2 = t1*t4 - t2*t6
        y3 = t2*t5 - t3*t4
        c2 = y1 + y2 + y3
        c1 = x1*y1 + x2*y2 + x3*y3
        t1 = x - x0
        if m < 0
          res = t1*(y0 + t0*t1*(0.5*c1 + t1*(0.25*c3*t1 - c2/3.0) ) )
        elseif m == 0
          res = y0 + t0*t1*(c1 + t1*(c3*t1 - c2))
        elseif m == 1
          res =t0*(c1+t1*(3.0*c3*t1 - 2.0*c2))
        elseif m == 2
          res = t0*(6.0*c3*t1 - 2.0*c2)
        elseif m == 3
          res = t0*6.0*c3
        else
          res = 0.0
        end
        return res
    end

    if m >= np
      res = 0.0
      return
    end

    c = zeros(Float64, np)

    # find the polynomial coefficients in divided differences form
    c[:] = ya[:]
    for i in 2:np
        for j in np:-1:i
            c[j] = ( c[j] - c[j-1])/( xa[j] - xa[j+1-i] )
        end
    end

    # special case m = 0
    if m == 0
      ss = c[1]
      t1 = 1.0
      for i in 2:np
        t1 = t1*(x - xa[i-1] )
        ss = ss + c[i]*t1
      end
      res = ss
      return res
    end
    
    x0 = xa[1]
    # convert to standard form
    for j in 1:np-1
        for i in 1:np-j
            k = np - i
            c[k] = c[k] + ( x0 - xa[k-j+1] )*c[k+1]
        end
    end
    
    if m > 0
        # take the m th derivative
        for j in 1:m
            for i in m+1:np
                c[i] = c[i]*(i-j)
            end
        end
        t1 = c[np]
        t2 = x - x0
        #for i in (np-1):-1:(m+1)
        for i in range(np-1, m+1, step=-1)
            t1 = t1*t2 + c[i]
        end
        res = t1
    else
        # find the integral
        t1 = c[np]/np
        t2 = x - x0
        for i in np-1:-1:1
          t1 = t1*t2 + c[i]/i
        end
        res = t1*t2
    end

    return res
end

