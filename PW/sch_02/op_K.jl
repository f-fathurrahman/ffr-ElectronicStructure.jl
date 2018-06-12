# Apply kinetic operator to wave function in reciprocal space

function op_K( pw::PWGrid, psi::Array{ComplexF64,2} )
    #
    out = zeros(ComplexF64,size(psi))
    Ncol = size(psi)[2]

    Ngwx = pw.gvecw.Ngwx
    G2   = pw.gvec.G2[pw.gvecw.idx_gw2r]

    for is = 1:Ncol
        for ig = 1:Ngwx
          out[ig,is] = psi[ig,is]*G2[ig]
        end
    end
    return 0.5*out # two minus signs
end
