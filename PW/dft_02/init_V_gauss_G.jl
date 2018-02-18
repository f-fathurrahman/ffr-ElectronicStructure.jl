function init_V_gauss_G(strf, Ω, Ns, G2, A, alpha)
    Npoints = size(strf)[1]
    Nspecies = size(strf)[2]
    Vg = zeros(Complex128,Npoints)
    V = zeros(Float64,Npoints)
    for isp = 1:Nspecies
        pf = -A[isp]*(pi/alpha[isp])^1.5/Ω
        @printf("Gaussian parameters:\n")
        @printf("A[%2d]     = %18.10f\n", isp, A[isp])
        @printf("alpha[%2d] = %18.10f\n", isp, alpha[isp])
        for ig = 1:Npoints
            Vg[ig] = pf*exp(-0.25*G2[ig]/alpha[isp])*strf[ig,isp]
        end
        V = V + real( G_to_R(Ns,Vg) )*Npoints
    end
    return V
end
