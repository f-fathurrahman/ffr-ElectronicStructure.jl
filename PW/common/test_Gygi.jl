
# r[3,Npoints]
# atpos[3,Natoms]
# A[Natoms]
# R0[Natoms]
function to_curvilinear( r::Array{Float64,2}, atpos, A, R0)
    #
    Npoints = size(r)[2]
    Natoms = size(atpos)[2]
    #
    ξ = zeros(Float64,3,Npoints)
    #
    for ip = 1:Npoints
        ξ[:,ip] = r[:,ip]
        for ia = 1:Natoms
            xx = norm( (r[:,ip] - atpos[:,ia])/R0[ia] )
            ξ[:,ip] = ξ[:,ip] + A[ia]*( r[:,ip] - atpos[:,ia] )/R0[ia] * sech(xx)
        end
    end
end
