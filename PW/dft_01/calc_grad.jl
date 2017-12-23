function calc_grad( pw::PWGrid, Potentials, Focc, psi::Array{Complex128,2} )

    Npoints = size(psi)[1]
    Nstates = size(psi)[2]
    Ω = pw.Ω
    Ns = pw.Ns
    #
    grad = zeros( Complex128, Npoints, Nstates )

    H_psi = op_H( pw, Potentials, psi )
    for i = 1:Nstates
        grad[:,i] = H_psi[:,i]
        for j = 1:Nstates
            grad[:,i] = grad[:,i] - dot( psi[:,j], H_psi[:,i] ) * psi[:,j]
        end
        grad[:,i] = Focc[i]*grad[:,i]
    end
    return grad

end
