# Generates exponential mesh of N elements on [r_mn, r_max]
#
# The domain [r_min, r_max], the mesh will contain both endpoints:
function gen_rmesh_exp(
    r_min::Float64,
    r_max::Float64,
    a::Float64,
    Nr::Int64
)
    N = Nr - 1

    rmesh = zeros(Float64,N+1)
    SMALL = eps()
    
    if a < 0.0
        error("gen_rmesh_exp: a > 0 required")
    #
    elseif abs(a - 1.0) < SMALL
        α = (r_max - r_min) / N
        for i in 1:N+1
            rmesh[i] = α * (i - 1) + r_min
        end
    #
    else
        if N > 1
            β = log(a)/(N-1)
            α = (r_max - r_min) / (exp(β*N) - 1)
            for i in 1:N+1
                rmesh[i] = α * ( exp(β*(i-1)) - 1 ) + r_min
            end
        elseif N == 1
            rmesh[1] = r_min
            rmesh[2] = r_max
        else
            error("gen_rmesh_exp: N >= 1 required")
        end
    end
    return rmesh
end
