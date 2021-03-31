function rsch_integ_rk4!(
    E::Float64,
    Z::Int64, l::Int64, rmesh, V, P, Q;
    MAXVAL_STOP=1e6
)
    Nr = size(rmesh,1)

    # Calculate Vmid
    Vmid = zeros(Float64, Nr-1)
    for i in 1:Nr-1
        rmid = 0.5*( rmesh[i] + rmesh[i+1] )
        Vmid[i] = radial_interp(rmesh, V, rmid, i+1)
    end

    #@printf("Some Vmid\n")
    #for i in 1:4
    #    @printf("%8d %20.10e\n", i, Vmid[i])
    #end
    #@printf("...\n")
    #for i in Nr-4:Nr-1
    #    @printf("%8d %20.10e\n", i, Vmid[i])
    #end

    # BC at r -> 0
    y0 = zeros(Float64,2)
    if l == 0
        y0[1] = 1.0 - Z*rmesh[1]
        y0[2] = -Z
    else
        y0[1] = rmesh[1]^l
        y0[2] = l*rmesh[1]^(l-1)
    end

    C1 = zeros(Float64,Nr)
    C2 = zeros(Float64,Nr)
    for i in 1:Nr
        C1[i] = 2*(V[i] - E) + l*(l + 1)/rmesh[i]
        C2[i] = -2/rmesh[i]
    end

    C1mid = zeros(Float64,Nr-1)
    C2mid = zeros(Float64,Nr-1)
    for i in 1:Nr-1
        rmid = 0.5*(rmesh[i] + rmesh[i+1])
        C1mid[i] = 2*(Vmid[i] - E) + l*(l + 1)/rmid^2
        C2mid[i] = -2/rmid
    end

    y1 = zeros(Float64,Nr)
    y2 = zeros(Float64,Nr)
    imax = sch_rk4_step!( rmesh, y0, C1, C2, C1mid, C2mid, y1, y2, MAXVAL_STOP )

    for i in 1:imax
        P[i] = y1[i]*rmesh[i] # P(r) = r*R(r)
        Q[i] = y2[i]*rmesh[i] + y1[i]
    end

    return imax
end