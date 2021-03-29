# Integrates the following set of equations outwards:
# dy1/dx =                y2
# dy2/dx = C1 * y1 + C2 * y2
function sch_rk4_step!(
    rmesh::Array{Float64,1}, # Grid
    y0::Array{Float64,1}, # Initial condition
    C1::Array{Float64,1},
    C2::Array{Float64,1},
    C1mid::Array{Float64,1},
    C2mid::Array{Float64,1},
    y1::Array{Float64,1},
    y2::Array{Float64,1}, # Solution y1 and y2
    maxval_stop::Float64,  # Maximum value (if y1 > max_val, the integration stops)
)

    Neqn = size(y0,1)  # should be hardcoded at 2 ?
    @assert Neqn == 2

    dym = zeros(Float64,2)
    dyt = zeros(Float64,2)
    yt = zeros(Float64,2)
    dydx = zeros(Float64,2)
    y = zeros(Float64, Neqn)

    N = size(rmesh,1)
    y = copy(y0)
    y1[1] = y[1]
    y2[1] = y[2]
    
    for i in 2:N
        # rk4 step size
        h = rmesh[i] - rmesh[i-1]

        # k1, evaluate F
        dydx[1] = y[2]
        dydx[2] = C1[i-1] * y[1] + C2[i-1] * y[2]
        #
        yt = y + h/2 * dydx  # yi + h/2*k1
        # k2, evaluate F at x_i + h/2 and yi + k1/2
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
        if abs(y[1]) >= maxval_stop
            imax = i
            return imax
        end
    end
    imax = n
    return imax
end