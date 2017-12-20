function calc_grad( LF, Potentials, Focc, v::Array{Float64,2} )
    #
    Npoints = size(v)[1]
    Ncol    = size(v)[2]
    ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
    grad = zeros( Float64,Npoints,Ncol )
    #
    for ic = 1:Ncol
        Hv = op_H( LF, Potentials, v[:,ic] )
        grad[:,ic] = Hv # copy?
        for icc = 1:Ncol
            grad[:,ic] = grad[:,ic] - dot( v[:,icc], Hv ) * v[:,icc] * ΔV
        end
        grad[:,ic] = Focc[ic] * grad[:,ic]
    end
    return grad
end


function calc_grad( LF, ∇2::SparseMatrixCSC{Float64,Int64},
                Potentials, Focc, v::Array{Float64,2} )
    #
    Npoints = size(v)[1]
    Ncol    = size(v)[2]
    ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
    grad = zeros(Float64,Npoints,Ncol)
    for ic = 1:Ncol
        Hv = op_H( LF, ∇2, Potentials, v[:,ic] )
        grad[:,ic] = Hv # copy?
        for icc = 1:Ncol
            grad[:,ic] = grad[:,ic] - dot( v[:,icc], Hv ) * v[:,icc] * ΔV
        end
        grad[:,ic] = Focc[ic] * grad[:,ic]
    end
    return grad
end
