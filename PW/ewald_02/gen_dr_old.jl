# center is (0,0,0)
# not working for general case
function gen_dr( r, LatVecs )
    #
    v1 = LatVecs[:,1]
    v2 = LatVecs[:,2]
    v3 = LatVecs[:,3]
    #
    Npoints = size(r)[2]
    dr = Array{Float64}(Npoints)
    #
    Lx = maximum( abs(LatVecs[:,1]) )
    Ly = maximum( abs(LatVecs[:,2]) )
    Lz = maximum( abs(LatVecs[:,3]) )
    @printf("VecLen: (%18.10f,%18.10f,%18.10f)\n", Lx,Ly,Lz)
    # r are guaranted to contained within the unit cell
    for ip=1:Npoints

        x = r[1,ip] < 0.0 ? r[1,ip] + Lx : r[1,ip]
        y = r[2,ip] < 0.0 ? r[2,ip] + Ly : r[2,ip]
        z = r[3,ip] < 0.0 ? r[3,ip] + Lz : r[3,ip]

        dx = x > 0.5*Lx ? Lx-x : x
        dy = y > 0.5*Ly ? Ly-y : y
        dz = z > 0.5*Lz ? Lz-z : z

        dr[ip] = sqrt(dx^2 + dy^2 + dz^2)
        @printf("\nWrapped: (%10.5f,%10.5f,%10.5f)\n", x, y, z);
        @printf("%8d (%10.5f,%10.5f,%10.5f) (%10.5f,%10.5f,%10.5f): %18.10f\n",
                    ip, r[1,ip], r[2,ip], r[3,ip],
                    dx, dy, dz, dr[ip] )
    end
    return dr
end
