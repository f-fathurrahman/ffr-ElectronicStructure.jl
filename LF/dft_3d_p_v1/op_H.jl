# one-column version
function op_H( LF, Potentials, v::Array{Float64,1} )
    #
    V_loc = Potentials.Ionic + Potentials.Hartree + Potentials.XC
    #
    return -0.5*op_nabla2( LF, v ) + V_loc .* v
end

# multi-column version
function op_H( LF, Potentials, v::Array{Float64,2} )
    #
    Npoints = size(v)[1]
    Ncols   = size(v)[2]
    V_loc = Potentials.Ionic + Potentials.Hartree + Potentials.XC
    #
    Hv = zeros( Float64, Npoints, Ncols )
    #
    for ic = 1:Ncols
        Hv[:,ic] = -0.5*op_nabla2( LF, v[:,ic] ) + V_loc .* v[:,ic]
    end
    return Hv
end


# one-column version, sparse
function op_H( LF, ∇2::SparseMatrixCSC{Float64,Int64},
                            Potentials, v::Array{Float64,1} )
    #
    V_loc = Potentials.Ionic + Potentials.Hartree + Potentials.XC
    return -0.5*∇2*v + V_loc.*v
end

# multi-column version, sparse
function op_H( LF, ∇2, Potentials, v::Array{Float64,2} )
    #
    Npoints = size(v)[1]
    Ncols   = size(v)[2]
    V_loc = Potentials.Ionic + Potentials.Hartree + Potentials.XC
    #
    Hv = zeros( Float64, Npoints, Ncols )
    #
    for ic = 1:Ncols
        Hv[:,ic] = -0.5*∇2*v[:,ic] + V_loc.*v[:,ic]
    end
    return Hv
end
