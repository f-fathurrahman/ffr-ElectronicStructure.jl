# one-column version
function op_H( LF, Vpot, v::Array{Float64,1} )
  # Is this OK?
  Hv = -0.5*op_nabla2( LF, v ) + Vpot .* v
end

# one-column version, sparse
function op_H( LF, ∇2::SparseMatrixCSC{Float64,Int64}, Vpot, v::Array{Float64,1} )
  # Is this OK?
  Hv = -0.5 * ∇2 * v + Vpot .* v
end

# multi-column version, sparse
function op_H( LF, ∇2::SparseMatrixCSC{Float64,Int64}, Vpot, v::Array{Float64,2} )
  Hv = -0.5 * ∇2 * v + Vpot .* v
end

# multicolumn version
function H_psi( LF, Vpot, psi::Array{Float64,2} )
  Ncol = size(psi)[2]
  Hpsi = zeros( Float64, Npoints, Ncol )
  for ic = 1:Ncol
    Hpsi[:,ic] = -0.5*op_nabla2( LF, psi[:,ic] ) + Vpot .* psi[:,ic]
  end
  return H_psi
end
