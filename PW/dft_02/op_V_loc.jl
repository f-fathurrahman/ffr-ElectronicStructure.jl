function op_V_loc( pw::PWGrid, V_loc, psi::Array{ComplexF64,2} )
    #
    Ns = pw.Ns
    Ω  = pw.Ω
    Npoints = prod(Ns)
    Ncols   = size(psi)[2]

    ctmp = zeros(ComplexF64, Npoints, Ncols)
    idx = pw.gvecw.idx_gw2r
    for ic = 1:Ncols
        ctmp[idx,ic] = psi[:,ic]
    end

    # get values of psi in real space grid via forward transform
    ctmp = G_to_R( Ns, ctmp )

    cVpsi = R_to_G( Ns, Diagprod(V_loc, ctmp) )
    return cVpsi[idx,:]
end

function op_V_loc( pw::PWGrid, V_loc, psi::Array{ComplexF64,1} )
    #
    Ns = pw.Ns
    Ω  = pw.Ω
    Npoints = prod(Ns)

    ctmp = zeros(ComplexF64, Npoints)
    idx = pw.gvecw.idx_gw2r
    ctmp[idx] = psi[:]

    # get values of psi in real space grid via forward transform
    ctmp = G_to_R( Ns, ctmp )

    cVpsi = R_to_G( Ns, V_loc .* ctmp )
    return cVpsi[idx]
end

# B is usually consists of more than one-column
function Diagprod( a, B::Array{ComplexF64,2} )
    Nstates = size(B)[2]
    Npoints = size(B)[1]
    out = zeros( ComplexF64, size(B) )
    for ist = 1:Nstates
        for ip = 1:Npoints
            out[ip,ist] = a[ip]*B[ip,ist]
        end
    end
    return out
end
