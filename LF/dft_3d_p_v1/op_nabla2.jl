# Apply ∇^2 operator
#
# Global vars:
#   Npoints, LF

function op_nabla2( LF, v::Array{Float64,1} )
    Npoints = LF.Nx * LF.Ny * LF.Nz
    #
    nabla2_v = zeros(Float64,Npoints)
    #
    for ip = 1 : Npoints
        i = LF.lin2xyz[1,ip]
        j = LF.lin2xyz[2,ip]
        k = LF.lin2xyz[3,ip]
        #
        for ii = 1 : LF.Nx
            nabla2_v[ip] = nabla2_v[ip] + LF.LFx.D2jl[ii,i] * v[LF.xyz2lin[ii,j,k]]
        end
        #
        for jj = 1 : LF.Ny
            nabla2_v[ip] = nabla2_v[ip] + LF.LFx.D2jl[jj,j] * v[LF.xyz2lin[i,jj,k]]
        end
        #
        for kk = 1 : LF.Nz
            nabla2_v[ip] = nabla2_v[ip] + LF.LFx.D2jl[kk,k] * v[LF.xyz2lin[i,j,kk]]
        end
    end
    #
    return nabla2_v
    #
end  # function
