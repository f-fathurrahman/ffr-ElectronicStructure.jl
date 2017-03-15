function calc_evals( LF, Vpot, evecs )
  Ncols   = size(evecs)[2]
  Npoints = size(evecs)[1]
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
  Hv = zeros(Float64,Npoints,Ncols)
  Hred = zeros(Float64,Ncols,Ncols)
  for ic = 1:Ncols
    Hv[:,ic] = op_H( LF, Vpot, evecs[:,ic] )
  end
  # Hred = evecs' * Hv
  BLAS.gemm!('T','N',1.0,evecs,Hv,0.0,Hred)
  #
  D, V = eig(Hred*ΔV)  # XXX deltaV is applied here
  return D
end


function get_evals( LF, ∇2::SparseMatrixCSC{Float64,Int64}, Vpot, evecs )
  Ncols   = size(evecs)[2]
  Npoints = size(Vpot)[1]
  Hv = zeros(Float64,Npoints,Ncols)
  Hred = zeros(Float64,Ncols,Ncols)
  for ic = 1:Ncols
    Hv[:,ic] = op_H( LF, ∇2, Vpot, evecs[:,ic] )
  end
  # Hred = evecs' * Hv
  BLAS.gemm!('T','N',1.0,evecs,Hv,0.0,Hred)
  #
  D, V = eig(Hred*ΔV)
  return D
end
