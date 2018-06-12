# distance to center of the simulation cell
function gen_dr( r, LatVecs )
    #
    cx = 0.5*sum(LatVecs[:,1])
    cy = 0.5*sum(LatVecs[:,2])
    cz = 0.5*sum(LatVecs[:,3])
    center = [cx,cy,cz]
    #@printf("Center: (%15.10f,%15.10f,%15.10f)\n", cx,cy,cz)

    Npoints = size(r)[2]
    dr = Array{Float64}(Npoints)

    for ip=1:Npoints
        dr[ip] = norm( r[:,ip] - center )
        #@printf("%8d (%15.10f,%15.10f,%15.10f): %18.10f\n", ip,
        #        r[1,ip], r[2,ip], r[3,ip], dr[ip])
    end
    return dr
end
