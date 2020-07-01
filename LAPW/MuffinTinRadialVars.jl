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

function MuffinTinRadialVars(Nspecies)

    nrmtscf = 0
    nrmt = zeros(Int64,Nspecies)
    nrmtmax = 0
    rmtall = 0
    rmtdelta = 0.05 # default
    rmt = zeros(Nspecies)
    omegamt = 0.0
    lradstp = 4
    nrcmt = zeros(Int64,Nspecies)
    nrcmtmax = 0
    rcmt = zeros(Float64,1,1)
    rlmt = zeros(Float64,1,1,1)
    rlcmt = zeros(Float64,1,1,1) 
    wrmt = zeros(Float64,1,1)
    wprmt = zeros(Float64,1,1,1)
    wrcmt = zeros(Float64,1,1) 
    wprcmt = zeros(Float64,1,1,1)
    maxlapw = 50
    lmaxapw = 8 # default
    lmmaxapw = (lmaxapw+1)^2
    lmaxo  = 6
    lmmaxo = (lmaxo+1)^2
    lmaxi  = 1 # min(lmaxi,lmaxo)
    lmmaxi = (lmaxi+1)^2
    fracinr = 0.01
    
    nrmti = zeros(Int64,Nspecies)
    nrcmti = zeros(Int64,Nspecies)
    
    idxlm = zeros(Int64,1,1)
    idxil = zeros(Int64,1)
    idxim = zeros(Int64,1)

    npmti = zeros(Int64,Nspecies)
    npmt = zeros(Int64,Nspecies)
    npcmti = zeros(Int64,Nspecies)
    npcmt = zeros(Int64,Nspecies)
    
    npmtmax = 0
    npcmtmax = 0

    return MuffinTinRadialVars(
        nrmtscf, nrmt, nrmtmax, rmtall, rmtdelta, rmt, omegamt, lradstp,
        nrcmt, nrcmtmax, rcmt, rlmt, rlcmt, wrmt, wprmt, wrcmt, wprcmt,
        maxlapw, lmaxapw, lmmaxapw, lmaxo, lmmaxo, lmaxi, lmmaxi, fracinr,
        nrmti, nrcmti, idxlm, idxil, idxim, npmti, npmt, npcmti, npcmt,
        npmtmax, npcmtmax,
    )

end


function init_zero!( mtr_vars::MuffinTinRadialVars )
    
    nrmt = mtr_vars.nrmt
    nspecies = size(nrmt)[1]

    nrcmt = mtr_vars.nrcmt
    lradstp = mtr_vars.lradstp
    
    # make the muffin-tin mesh commensurate with lradstp
    for is in 1:nspecies
        nrmt[is] = nrmt[is] - (nrmt[is]-1)%lradstp
        nrcmt[is] =( nrmt[is] - 1)/lradstp + 1
    end
    mtr_vars.nrmtmax = maximum(nrmt)
    mtr_vars.nrcmtmax = maximum(nrcmt)

    @assert mtr_vars.lmaxo <= mtr_vars.lmaxapw

    return
end


# Should be called after genrmesh
function init_packed_mtr!( mtr_vars::MuffinTinRadialVars )
    #
    nrmti = mtr_vars.nrmti
    nspecies = size(nrmti,1)
    nrmt = mtr_vars.nrmt
    #
    nrcmt = mtr_vars.nrcmt
    nrcmti = mtr_vars.nrcmti
    #
    lmmaxo = mtr_vars.lmmaxo
    lmmaxi = mtr_vars.lmmaxi
    #
    mtr_vars.npmtmax  = 1
    mtr_vars.npcmtmax = 1
    #

    for is in 1:nspecies
        #
        mtr_vars.npmti[is] = lmmaxi*nrmti[is]
        mtr_vars.npmt[is] = mtr_vars.npmti[is] + lmmaxo*(nrmt[is] - nrmti[is])
        mtr_vars.npmtmax = max(mtr_vars.npmtmax, mtr_vars.npmt[is])
        #
        mtr_vars.npcmti[is] = lmmaxi*nrcmti[is]
        mtr_vars.npcmt[is] = mtr_vars.npcmti[is] + lmmaxo*(nrcmt[is] - nrcmti[is])
        mtr_vars.npcmtmax = max(mtr_vars.npcmtmax, mtr_vars.npcmt[is])
    end
    return
end