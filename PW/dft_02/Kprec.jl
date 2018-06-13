function Kprec( pw::PWGrid, psi::Array{ComplexF64,2} )

    Ngwx  = size(psi)[1]
    Nstates = size(psi)[2]
    G2    = pw.gvec.G2[pw.gvecw.idx_gw2r]
    Kpsi  = zeros( ComplexF64, Ngwx, Nstates )

    for ist = 1:Nstates
        for ig = 1:Ngwx
            Kpsi[ig,ist] = psi[ig,ist] / ( 1.0 + G2[ig] )
        end
    end
    return Kpsi
end

# For comparison with non-preconditioned
function Kprec(psi)
    return psi
end
