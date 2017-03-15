function gradE( LF, Vpot, v::Array{Float64,2} )
  Npoints = size(Vpot)[1]
  Ncol = size(v)[2]
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
  grad = zeros( Float64,Npoints,Ncol )
  #
  for ic = 1:Ncol
    Hv = op_H( LF, Vpot, v[:,ic] )
    grad[:,ic] = Hv # copy?
    for icc = 1:Ncol
      grad[:,ic] = grad[:,ic] - dot( v[:,icc], Hv ) * v[:,icc] * ΔV
    end
  end
  return grad
end


function gradE( LF, ∇2::SparseMatrixCSC{Float64,Int64}, Vpot, v::Array{Float64,2} )
  Npoints = size(Vpot)[1]
  Ncol = size(v)[2]
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
  grad = zeros(Float64,Npoints,Ncol)
  for ic = 1:Ncol
    Hv = op_H( LF, ∇2, Vpot, v[:,ic] )
    grad[:,ic] = Hv # copy?
    for icc = 1:Ncol
      grad[:,ic] = grad[:,ic] - dot( v[:,icc], Hv ) * v[:,icc] * ΔV
    end
  end
  return grad
end
