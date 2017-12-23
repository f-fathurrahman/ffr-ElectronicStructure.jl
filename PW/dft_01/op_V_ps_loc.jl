function op_V_ps_loc( pw::PWGrid, V_loc, psi::Array{Complex128,2} )
    #
    Ns = pw.Ns
    Ω  = pw.Ω
    Npoints = prod(Ns)
    # get values of psi in real space grid via forward transform
    ctmp = G_to_R( Ns, psi )
    return R_to_G( Ns, Diagprod(V_loc, ctmp) )
end

# B is usually consists of more than one-column
function Diagprod( a,B )
    Nstates    = size(B)[2]
    Npoints = size(B)[1]
    out = zeros( Complex128, size(B) )
    for ist = 1:Nstates
        for ip = 1:Npoints
            out[ip,ist] = a[ip]*B[ip,ist]
        end
    end
    return out
end
