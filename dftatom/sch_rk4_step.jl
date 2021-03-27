function sch_rk4_step!(
    R, # Grid
    y0, # Initial condition
    C1, C2, C1mid, C2mid, # Coefficients C1 and C2 at grid points and midpoints:
    max_val,  # Maximum value (if y1 > max_val, the integration stops)
    y1, y2 # Solution y1 and y2
)
    # Integrates the following set of equations outwards:
    # dy1/dx =                y2
    # dy2/dx = C1 * y1 + C2 * y2

    Neqn = size(y0,1)  # should be hardcoded at 2 ?
    dym = zeros(Float64, Neqn)
    dyt = zeros(Float64, Neqn)
    yt = zeros(Float64, Neqn)
    dydx = zeros(Float64, Neqn)
    y = zeros(Float64, Neqn)

    N = size(R,1)
    y = copy(y0)
    y1[1] = y[1]
    y2[1] = y[2]
    
    for i in 2:N
        # rk4 step size
        h = R[i] - R[i-1]

        # k1
        dydx[1] = y[2]
        dydx[2] = C1[i-1] * y[1] + C2[i-1] * y[2]
        #
        yt = y + h/2 * dydx  # yi + h/2*k1
        # k2
        dyt[1] = yt[2]
        dyt[2] = C1mid[i-1] * yt[1] + C2mid[i-1] * yt[2]
        #
        yt = y + h/2 * dyt # yi + h/2*k2
        #
        dym[1] = yt[2]
        dym[2] = C1mid[i-1] * yt[1] + C2mid[i-1] * yt[2]
        #
        yt = y + h * dym
        dym = dyt + dym
        dyt[1] = yt[2]
        dyt[2] = C1[i] * yt[1] + C2[i] * yt[2]
        #
        y = y + h/6 * (dydx + dyt + 2*dym)
        #
        y1[i] = y[1]
        y2[i] = y[2]

        # The integration stops at R(imax)
        if abs(y[1]) >= max_val
            imax = i
            return imax
        end
    end
    imax = n
    return imax
end