function gen_drmesh_exp(r_min, r_max, a, Nr)

    N = Nr - 1

    drmesh = zeros(Float64,N+1)
    SMALL = eps()

    if a < 0
        error("gen_drmesh_exp: a > 0 required")
    #
    elseif abs(a - 1) < SMALL
        error("gen_drmesh_exp: a == 1 not implemented")
    #
    else
        if N > 1
            β = log(a)/(N-1)
            α = (r_max - r_min) / (exp(β*N) - 1)
            for i in 1:N+1
                drmesh[i] = α * β * exp( β*(i-1) )
            end
        else
            error("mesh_exp_deriv: N > 1 required")
        end
    end
    return drmesh
end
