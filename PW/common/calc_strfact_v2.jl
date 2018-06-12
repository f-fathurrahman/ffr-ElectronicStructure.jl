# Calculate structure factor
function calc_strfact( Xpos::Array{Float64,2}, Nspecies::Int,
    atm2species::Array{Int}, G::Array{Float64,2} )
    #
    Ng = size(G)[2]
    Na = size(Xpos)[2]
    Sf = zeros(ComplexF64,Ng,Nspecies)
    for ia = 1:Na
        isp = atm2species[ia]
        for ig = 1:Ng
            GX = Xpos[1,ia]*G[1,ig] + Xpos[2,ia]*G[2,ig] + Xpos[3,ia]*G[3,ig]
            Sf[ig,isp] = Sf[ig,isp] + cos(GX) - im*sin(GX)
        end
    end
    return Sf
end
