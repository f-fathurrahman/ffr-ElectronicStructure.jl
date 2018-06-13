# Apply kinetic operator to wave function in reciprocal space

function op_K( pw::PWGrid, psi::Array{ComplexF64,2}, ik::Int64 )
    #
    out = zeros(ComplexF64,size(psi))
    Ncol = size(psi)[2]

    Ngwk = pw.gkvec.Ngw[ik]
    idx_gkk = pw.gkvec.idx_gk[ik]
    G = pw.gvec.G[:,idx_gkk]

    for is = 1:Ncol
        for ig = 1:Ngwk
            Gk = G[:,ig] + pw.gkvec.kpts[:,ik]
            out[ig,is] = psi[ig,is]*( Gk[1]^2 + Gk[2]^2 + Gk[3]^2 )
        end
    end
    return 0.5*out # two minus signs
end
