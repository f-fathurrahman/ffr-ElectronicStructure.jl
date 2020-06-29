mutable struct MuffinTinRadialVars
    # scale factor for number of muffin-tin points
    nrmtscf::Float64
    # number of muffin-tin radial points for each species
    nrmt::Vector{Int64}
    # maximum nrmt over all the species
    nrmtmax::Int64
    # optional default muffin-tin radius for all atoms
    rmtall::Float64
    # minimum allowed distance between muffin-tin surfaces
    rmtdelta::Float64
    # muffin-tin radii
    rmt::Vector{Float64}
    # total muffin-tin volume
    omegamt::Float64
    # radial step length for coarse mesh
    lradstp::Int64
    # number of coarse radial mesh points
    nrcmt::Vector{Int64}
    # maximum nrcmt over all the species
    nrcmtmax::Int64
    # coarse muffin-tin radial mesh
    rcmt::Array{Float64,2}
    # r^l on fine radial mesh
    rlmt::Array{Float64,3}
    # r^l on coarse radial mesh
    rlcmt::Array{Float64,3}
    # weights for spline integration on fine radial mesh
    wrmt::Array{Float64,2}
    # weights for spline partial integration on fine radial mesh
    wprmt::Array{Float64,3}
    # weights for spline integration on coarse radial mesh
    wrcmt::Array{Float64,2}
    # weights for spline partial integration on coarse radial mesh
    wprcmt::Array{Float64,3}
    # maximum allowable angular momentum for augmented plane waves
    maxlapw::Int64 # parameter =50
    # maximum angular momentum for augmented plane waves
    lmaxapw::Int64
    # (lmaxapw+1)^2
    lmmaxapw::Int64
    # maximum angular momentum on the outer part of the muffin-tin
    lmaxo::Int64
    # (lmaxo+1)^2
    lmmaxo::Int64
    # maximum angular momentum on the inner part of the muffin-tin
    lmaxi::Int64
    # (lmaxi+1)^2
    lmmaxi::Int64
    # fraction of muffin-tin radius which constitutes the inner part
    fracinr::Float64
    # number of fine/coarse radial points on the inner part of the muffin-tin
    nrmti::Vector{Int64}
    nrcmti::Vector{Int64} 
    # index to (l,m) pairs
    idxlm::Array{Int64,2}
    # inverse index to (l,m) pairs
    idxil::Vector{Int64}
    idxim::Vector{Int64}
    # number of fine/coarse points in packed muffin-tins
    npmti::Vector{Int64}
    npmt::Vector{Int64}
    npcmti::Vector{Int64}
    npcmt::Vector{Int64}
    # maximum number of points over all packed muffin-tins
    npmtmax::Int64
    npcmtmax::Int64
end