# distance to center of the simulation cell
function gen_dr_center( r, LatVecs )
    #
    cx = 0.5*sum(LatVecs[:,1])
    cy = 0.5*sum(LatVecs[:,2])
    cz = 0.5*sum(LatVecs[:,3])
    center = [cx,cy,cz]
    #
    Npoints = size(r)[2]
    dr = Array{Float64}(undef,Npoints)
    for ip=1:Npoints
        dr[ip] = norm( r[:,ip] - center )
    end
    return dr
end
