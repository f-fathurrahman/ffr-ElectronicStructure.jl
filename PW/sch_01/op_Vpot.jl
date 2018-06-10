function op_Vpot( pw_grid::PWGrid, Vpot, psi::Array{ComplexF64,2} )
    #
    Ns = pw_grid.Ns
    Ω = pw_grid.Ω
    Npoints = prod(Ns)
    # get values of psi in real space grid via forward transform
    ctmp = G_to_R( Ns, psi )
    return R_to_G( Ns, Diagprod(Vpot, ctmp) )
end

# B is usually consists of more than one-column
function Diagprod( a,B )
    Ncol    = size(B)[2]
    Npoints = size(B)[1]
    out = zeros( ComplexF64, size(B) )
    for ic = 1:Ncol
        for ip = 1:Npoints
            out[ip,ic] = a[ip]*B[ip,ic]
        end
    end
    return out
end
