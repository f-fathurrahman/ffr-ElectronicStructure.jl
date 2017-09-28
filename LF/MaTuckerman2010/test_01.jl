function init_s_grid(Ns)
    if any( (Ns .% 2) .== 0 )
        error("Ns must be odd numbers")
    end
    NN = round.(Int, (Ns-1)/2)
    Npoints = prod(Ns)
    s = zeros(3,Npoints)
    ip = 0
    for n = -NN[3]:NN[3]
        for m = -NN[2]:NN[2]
            for l = -NN[1]:NN[1]
                ip = ip + 1
                s[1,ip] = l/(2*NN[1] + 1)
                s[2,ip] = m/(2*NN[2] + 1)
                s[3,ip] = n/(2*NN[3] + 1)
            end
        end
    end
    return s
end


# r = hh * s
function s_to_r(s::Array{Float64,2}, hh)
    Npoints = size(s)[2]
    r = zeros(3,Npoints)
    for ip = 1:Npoints
        for α = 1:3
            r[α,ip] = 0.0
            for β = 1:3
                r[α,ip] = r[α,ip] + hh[α,β]*s[β,ip]
            end
        end
    end
    return r
end


# Function to evaluate LF
function eval_μ_lmn(Ns, hh, s_grid, ip, s)
    #
    NN = round.(Int, (Ns-1)/2)
    metric = 1/sqrt(det(hh))
    #
    s_lmn = s_grid[:,ip]
    #
    f = 0.0 + im*0.0
    # loop structure should be the same as the one used in init_s_grid
    for kz = -NN[3]:NN[3]
        for ky = -NN[2]:NN[2]
            for kx = -NN[1]:NN[1]
                ks = dot( [kx,ky,kz], (s - s_lmn) )
                f = f + cos(2π*ks)
            end
        end
    end
    return metric/sqrt(prod(Ns))*f
end

include("../../PW/common/gen_lattice.jl")

function test_main()
    # Sampling points
    Ns = [21,21,21]

    LL = gen_lattice_hexagonal(10)
    #LL = 10.0*eye(3,3)  # cubic
    hh = LL'

    s_grid = init_s_grid(Ns)

    ip = 2

    #r = [5.0, 5.0, 5.0]
    #s = inv(hh)*r
    s = s_grid[:,1]
    println( eval_μ_lmn(Ns, hh, s_grid, ip, s) )

    s = s_grid[:,2]
    println( eval_μ_lmn(Ns, hh, s_grid, ip, s) )
end

test_main()
