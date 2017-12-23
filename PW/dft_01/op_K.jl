function op_K( pw::PWGrid, psi::Array{Complex128,2} )
    out = zeros(Complex128,size(psi))
    Ncol = size(psi,2)
    Ω  = pw.Ω
    G2 = pw.G2
    Npoints = pw.Npoints
    for is = 1:Ncol
        for ip = 1:Npoints
            out[ip,is] = psi[ip,is]*G2[ip]
        end
    end
    return 0.5*out # two minus signs
end
