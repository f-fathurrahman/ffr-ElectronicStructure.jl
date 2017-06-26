function calc_evals( LF, Potentials, evecs )
  Ncols   = size(evecs)[2]
  Npoints = size(evecs)[1]
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
  Hred = zeros(Float64,Ncols,Ncols)
  Hv = op_H( LF, Potentials, evecs )
  # Hred = evecs' * Hv
  BLAS.gemm!('T','N',1.0,evecs,Hv,0.0,Hred)
  #
  D, V = eig(Hred*ΔV)  # XXX deltaV is applied here
  return D
end


function calc_evals( LF, ∇2::SparseMatrixCSC{Float64,Int64}, Potentials, evecs )
  Ncols   = size(evecs)[2]
  Npoints = size(evecs)[1]
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
  Hv = zeros(Float64,Npoints,Ncols)
  Hred = zeros(Float64,Ncols,Ncols)
  Hv = op_H( LF, ∇2, Potentials, evecs )
  # Hred = evecs' * Hv
  BLAS.gemm!('T','N',1.0,evecs,Hv,0.0,Hred)
  #
  D, V = eig(Hred*ΔV)
  return D
end
