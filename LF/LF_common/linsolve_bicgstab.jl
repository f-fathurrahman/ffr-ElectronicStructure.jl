function linsolve_bicgstab( LF::LF3dGrid, b::Array{Float64,1};
                            x0 = nothing,
                            NiterMax = 1000, TOL=5.e-10,
                            convmsg=false, showprogress=false )
    #
    Npoints = size(b)[1]
    #
    if x0 == nothing
        srand(1234)
        x = randn( Npoints )
    else
        x = copy(x0)
    end
    #
    r = b - apply_Laplacian(LF, x)
    #
    r_hat0 = copy(r)
    #
    rho = 1.0
    α   = 1.0
    ω   = 1.0
    #
    v = zeros(Float64,Npoints)
    p = zeros(Float64,Npoints)

    r_old = copy(r)
    v_old = zeros(Float64,Npoints)
    p_old = zeros(Float64,Npoints)
    x_old = copy(x)
    #
    rho_old = rho
    ω_old = ω
    for iter = 1:NiterMax
        rho = dot( r_hat0, r_old )
        β = (rho/rho_old)*(α/ω_old)
        p = r_old + β*( p_old - ω_old*v_old )
        v = apply_Laplacian(LF,p)
        α = rho/dot(r_hat0,v)
        h = x_old + α*p
        diffV_1 = norm((h-x))
        #@printf("iter, diffV_1 = %8d %e\n", iter, diffV_1)
        s = r_old - α*v
        t = apply_Laplacian(LF,s)
        ω = dot(t,s)/dot(t,t)
        x[:] = h[:] + ω*s[:]
        #
        diffV_2 = norm(x-x_old) #sum( abs(x-x_old) )
        @printf("iter, diffV_2 = %8d %e\n", iter, diffV_2)
        if diffV_1 < TOL || diffV_2 < TOL
            @printf("Convergence achieved in linsolve_bicgstab: %8d\n", iter)
            break
        end
        r = s - ω*t
        #
        rho_old = rho
        r_old = copy(r)
        v_old = copy(v)
        p_old = copy(p)
        x_old = copy(x)
    end
    return x
end
