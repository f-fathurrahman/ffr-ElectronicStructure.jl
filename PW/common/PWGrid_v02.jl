const PWGRID_VERSION = 2

mutable struct GVectorsW
    Ngwx::Int
    idx_gw2r::Array{Int}
end

mutable struct GVectors
    Ng::Int
    G::Array{Float64,2}
    G2::Array{Float64}
end

mutable struct PWGrid
    Ns::Array{Int64}
    LatVecs::Array{Float64,2}
    RecVecs::Array{Float64,2}
    Ω::Float64
    r::Array{Float64,2}
    gvec::GVectors
    gvecw::GVectorsW
end

function PWGrid( Ns, LatVecs::Array{Float64,2} )

    if any( Ns .% 2 .== 1 )
        @printf("Error: Ns must be even numbers\n")
        exit()
    end

    RecVecs = 2*pi*inv(LatVecs')
    Ω = det(LatVecs)
    #
    LatVecsLen = Array{Float64}(undef,3)
    LatVecsLen[1] = norm(LatVecs[:,1])
    LatVecsLen[2] = norm(LatVecs[:,2])
    LatVecsLen[3] = norm(LatVecs[:,3])

    Npoints = prod(Ns)
    r = init_grid_R( Ns, LatVecs )

    gvec = init_grid_G( Ns, RecVecs )
    gvecw = init_gvecw( Ns, gvec.G2 )

    return PWGrid( Ns, LatVecs, RecVecs, Ω, r, gvec, gvecw )
end

function mm_to_nn(mm::Int,S::Int)
    if mm > S/2
        return mm - S
    else
        return mm
    end
end


function init_grid_G( Ns, RecVecs )

    Ng = prod(Ns)

    G  = Array{Float64}(undef,3,Ng)
    G2 = Array{Float64}(undef,Ng)

    ig = 0
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ig = ig + 1
        gi = mm_to_nn( i, Ns[1] )
        gj = mm_to_nn( j, Ns[2] )
        gk = mm_to_nn( k, Ns[3] )
        G[1,ig] = RecVecs[1,1]*gi + RecVecs[1,2]*gj + RecVecs[1,3]*gk
        G[2,ig] = RecVecs[2,1]*gi + RecVecs[2,2]*gj + RecVecs[2,3]*gk
        G[3,ig] = RecVecs[3,1]*gi + RecVecs[3,2]*gj + RecVecs[3,3]*gk
        G2[ig] = G[1,ig]^2 + G[2,ig]^2 + G[3,ig]^2
    end
    end
    end

    return GVectors( Ng, G, G2 )
end

function init_gvecw( Ns, G2 )
    edges = find_edges( Ns )
    G2mx = minimum( G2[edges] )
    #@printf("G2mx = %f\n", G2mx)
    @printf("ecutwfc_Ry = %f\n", G2mx/4)
    idx_gw2r = findall( G2 .< G2mx/4 )
    Ngwx = length(idx_gw2r)
    return GVectorsW( Ngwx, idx_gw2r )
end

function find_edges( Ns )
    eS = Ns/2 .+ 0.5
    edges = []
    ip = 1
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ei = abs(i - eS[1])
        ej = abs(j - eS[2])
        ek = abs(k - eS[3])
        # if any of i, j, or k is equal to Ns/2 or Ns/2+1
        if ei < 1.0 || ej < 1.0 || ek < 1.0
            push!(edges,ip)
        end
        ip = ip + 1
    end
    end
    end
    return edges
end


function init_grid_R( Ns, LatVecs )
    #
    Npoints = prod(Ns)
    #
    R = Array{Float64}(undef,3,Npoints)
    ip = 0
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ip = ip + 1
        R[1,ip] = LatVecs[1,1]*i/Ns[1] + LatVecs[1,2]*j/Ns[2] + LatVecs[1,3]*k/Ns[3]
        R[2,ip] = LatVecs[2,1]*i/Ns[1] + LatVecs[2,2]*j/Ns[2] + LatVecs[2,3]*k/Ns[3]
        R[3,ip] = LatVecs[3,1]*i/Ns[1] + LatVecs[3,2]*j/Ns[2] + LatVecs[3,3]*k/Ns[3]
    end
    end
    end
    #
    return R
end
