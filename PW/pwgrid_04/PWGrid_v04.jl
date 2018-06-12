const PWGRID_VERSION = 4

mutable struct GkVectors
    kpts::Array{Float64,2}
    Ngwx::Int
    Ngw::Array{Int,1}
    idx_gk::Array{Array{Int64,1},1}
end

mutable struct GVectors
    Ng::Int
    G::Array{Float64,2}
    G2::Array{Float64}
end

mutable struct PWGrid
    ecutwfc::Float64
    ecutrho::Float64
    Ns::Array{Int64}
    LatVecs::Array{Float64,2}
    RecVecs::Array{Float64,2}
    Ω::Float64
    r::Array{Float64,2}
    gvec::GVectors
    gkvec::GkVectors
end

function PWGrid( ecutwfc::Float64, LatVecs::Array{Float64,2}, kpts_red )
    ecutrho = 4.0*ecutwfc
    #
    println("LatVecs = ", LatVecs)
    RecVecs = 2*pi*inv(LatVecs)'
    println("RecVecs = ", RecVecs)
    Ω = det(LatVecs)
    #
    LatVecsLen = Array{Float64}(undef,3)
    LatVecsLen[1] = norm(LatVecs[1,:])
    LatVecsLen[2] = norm(LatVecs[2,:])
    LatVecsLen[3] = norm(LatVecs[3,:])

    Ns = Array{Int64}(undef,3)
    Ns[1] = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[1]/pi ) + 1
    Ns[2] = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[2]/pi ) + 1
    Ns[3] = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[3]/pi ) + 1

    # Use even sampling numbers
    Ns[1] = Ns[1] % 2 == 1 ? Ns[1] + 1 : Ns[1]
    Ns[2] = Ns[2] % 2 == 1 ? Ns[2] + 1 : Ns[2]
    Ns[3] = Ns[3] % 2 == 1 ? Ns[3] + 1 : Ns[3]

    Npoints = prod(Ns)
    r = init_grid_R( Ns, LatVecs )

    gvec = init_grid_G( Ns, RecVecs )

    #kpts = (kpts_red' * RecVecs)'

    kpts = RecVecs*kpts_red # new following convention of PWDFT.jl

    #kpts = (kpts_red' * RecVecs')'

    #kpts = RecVecs' * kpts_red
    #kpts = RecVecs * kpts_red

    #Nkpts = size(kpts_red)[2]
    #kpts = zeros(Float64, 3, Nkpts)
    #for ik = 1:Nkpts
    #    kpts[:,ik] = RecVecs' * kpts_red[:,ik]
    #end
    #for ik = 1:Nkpts
    #    for i = 1:3
    #        kpts[i,ik] = kpts_red[i,ik] * RecVecs[1,:]
    #    end
    #end
    gkvec = init_gkvec( ecutwfc, gvec.G, kpts )

    return PWGrid( ecutwfc, ecutrho, Ns, LatVecs, RecVecs, Ω, r, gvec, gkvec )
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

function init_gkvec( ecutwfc, G, kpts )
    Nkpts = size(kpts)[2]
    idx_gk = Array{Array{Int64,1},1}(undef,Nkpts)
    Ng = size(G)[2]
    Gk = zeros( 3, Ng )
    Gk2 = zeros( Ng )
    Ngw = zeros( Int64, Nkpts )
    for ik = 1:Nkpts
        #@printf("\n")
        for ig = 1:Ng
            Gk[:,ig] = G[:,ig] + kpts[:,ik]
            Gk2[ig]  = Gk[1,ig]^2 + Gk[2,ig]^2 + Gk[3,ig]^2
            #@printf("%8d %18.10f %18.10f\n", ig, Gk2[ig], ecutwfc)
        end
        idx_gk[ik] = findall( 0.5*Gk2 .< ecutwfc )
        Ngw[ik] = length(idx_gk[ik])
        #@printf("kpt = %d, Ngw = %d\n", ik, Ngw[ik])
        #@printf("kpts[%d]: %f %f %f\n", ik, kpts[1,ik], kpts[3,ik], kpts[3,ik] )
    end
    Ngwx = maximum(Ngw)
    return GkVectors( kpts, Ngwx, Ngw, idx_gk )
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
