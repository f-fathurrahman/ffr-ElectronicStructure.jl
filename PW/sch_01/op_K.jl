function op_K( pw_grid::PWGrid, psi::Array{ComplexF64,2} )
    out = zeros(ComplexF64,size(psi))
    Ncol = size(psi,2)
    Ω = pw_grid.Ω
    G2 = pw_grid.G2
    Npoints = pw_grid.Npoints
    for is = 1:Ncol
        for ip = 1:Npoints
            out[ip,is] = psi[ip,is]*G2[ip]
        end
    end
    return 0.5*out # two minus signs
end

function op_K( pw_grid::PWGrid, psi::Array{ComplexF64,1} )
    G2 = pw_grid.G2
    Npoints = pw_grid.Npoints
    out = zeros( size(psi) )
    for ip = 1:Npoints
        out[ip] = psi[ip]*G2[ip]
    end
    return 0.5*out
end
