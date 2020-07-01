# Generates the coarse and fine radial meshes for each atomic species in the
# crystal. Also determines which points are in the inner part of the
# muffin-tin using the value of {\tt fracinr}.
function genrmesh!(atm_vars, atsp_vars, mtr_vars)

    nspecies = atm_vars.nspecies

    nrsp = atsp_vars.nrsp
    rminsp = atsp_vars.rminsp
    rmaxsp = atsp_vars.rmaxsp
    
    rmt = mtr_vars.rmt
    nrmt = mtr_vars.nrmt

    # estimate the number of radial mesh points to infinity
    atsp_vars.nrspmax = 1
    for is in 1:nspecies
        #println("rminsp = ", rminsp[is])
        #println("rmaxsp = ", rmaxsp[is])
        #println("nrmt   = ", nrmt[is])
        # logarithmic mesh
        t1 = log(rmaxsp[is]/rminsp[is])/log(rmt[is]/rminsp[is])
        t2 = (nrmt[is] - 1)*t1
        #println("t1 = ", t1)
        #println("t2 = ", t2)
        nrsp[is] = round(Int64,t2) + 1 # XXX compare with nint
        atsp_vars.nrspmax = max(atsp_vars.nrspmax, nrsp[is])
    end

    nrspmax = atsp_vars.nrspmax
    lmaxo = mtr_vars.lmaxo
    nrmtmax = mtr_vars.nrmtmax
    #println("lmaxo = ", lmaxo)
    #println("nrmtmax = ", nrmtmax)

    # The following actually allocate memory

    atsp_vars.rsp = zeros(Float64, nrspmax, nspecies)
    mtr_vars.rlmt = zeros(Float64, nrmtmax, 2*(lmaxo+2), nspecies)
    mtr_vars.wrmt = zeros(Float64, nrmtmax, nspecies)
    mtr_vars.wprmt = zeros(Float64, 4, nrmtmax, nspecies)

    rsp = atsp_vars.rsp
    rlmt = mtr_vars.rlmt
    wrmt = mtr_vars.wrmt
    wprmt = mtr_vars.wprmt

    #println("2*(lmaxo+2) = ", 2*(lmaxo+2))
    #println("size(rlmt) = ", size(rlmt))

    for is in 1:nspecies
        t1 = 1.0/(nrmt[is] - 1)
        # logarithmic mesh
        t2 = log(rmt[is]/rminsp[is])
        for ir in 1:nrsp[is]
            rsp[ir,is] = rminsp[is]*exp( (ir-1) * t1 * t2)
        end
        # calculate r^l on the fine radial mesh
        nr = nrmt[is]
        # XXX We are using lmaxo2idx
        il1 = lmaxo2idx(-1,lmaxo)
        il2 = lmaxo2idx( 0,lmaxo)
        il3 = lmaxo2idx( 1,lmaxo)
        for ir in 1:nr
            rlmt[ir,il1,is] = 1.0/rsp[ir,is]
            rlmt[ir,il2,is] = 1.0
            rlmt[ir,il3,is] = rsp[ir,is]
        end
        #
        for l in range(-2,stop=-lmaxo-1,step=-1)
            il1 = lmaxo2idx(l,lmaxo)
            il2 = lmaxo2idx(l+1,lmaxo)
            for ir in 1:nr
                rlmt[ir,il1,is] = rlmt[ir,il2,is]/rsp[ir,is]
            end
        end
        #
        for l in 2:lmaxo+2            
            il1 = lmaxo2idx(l,lmaxo)
            il2 = lmaxo2idx(l-1,lmaxo)
            for ir in 1:nr
                rlmt[ir,il1,is] = rlmt[ir,il2,is] * rsp[ir,is]
            end
        end
        # determine the weights for spline integration on the fine radial mesh
        @views wsplint!(nr, rsp[:,is], wrmt[:,is])
        # multiply by r^2
        il2 = lmaxo2idx(2,lmaxo)
        for ir in 1:nr
            wrmt[ir,is] = wrmt[ir,is]*rlmt[ir,il2,is]
        end
        # determine the weights for partial integration on fine radial mesh
        @views wsplintp!(nr, rsp[:,is], wprmt[:,:,is])
    end


    # determine the fraction of the muffin-tin radius which defines the inner part
    if mtr_vars.fracinr < 0.0
        # be explicit about conversion to Float64 (not really needed actually for Julia)
        mtr_vars.fracinr = sqrt( Float64(lmmaxi) / Float64(lmmaxo) )
    end

    fracinr = mtr_vars.fracinr
    nrcmtmax = mtr_vars.nrcmtmax
    nrcmt = mtr_vars.nrcmt
    nrmti = mtr_vars.nrmti
    nrcmti = mtr_vars.nrcmti
    lradstp = mtr_vars.lradstp
    println("nrcmtmax = ", nrcmtmax)
    println("nrcmt    = ", nrcmt[1:nspecies])

    # set up the coarse radial meshes and find the inner part of the muffin-tin
    # where rho is calculated with lmaxi
    mtr_vars.rcmt = zeros(Float64, nrcmtmax, nspecies)
    mtr_vars.rlcmt = zeros(Float64, nrcmtmax, 2*(lmaxo+2), nspecies)
    mtr_vars.wrcmt = zeros(Float64, nrcmtmax, nspecies)
    mtr_vars.wprcmt = zeros(Float64, 4, nrcmtmax, nspecies)

    rcmt = mtr_vars.rcmt
    rlcmt = mtr_vars.rlcmt
    wrcmt = mtr_vars.wrcmt
    wprcmt = mtr_vars.wprcmt

    for is in 1:nspecies
        t1 = fracinr*rmt[is]
        nrmti[is] = 1
        nrcmti[is] = 1
        irc = 0
        for ir in range(1,stop=nrmt[is],step=lradstp)
            irc = irc + 1
            rcmt[irc,is] = rsp[ir,is]
            if rsp[ir,is] < t1
                nrmti[is] = ir
                nrcmti[is] = irc
            end
        end
        # store r^l on the coarse radial mesh
        for l in range(-lmaxo-1, stop=lmaxo+2)
            irc = 0
            il1 = lmaxo2idx(l,lmaxo)
            for ir in range(1,stop=nrmt[is],step=lradstp)
                irc = irc + 1
                rlcmt[irc,il1,is] = rlmt[ir,il1,is]
            end
        end
        # determine the weights for spline integration on the coarse radial mesh
        nrc = nrcmt[is]
        @views wsplint!(nrc, rcmt[:,is], wrcmt[:,is])
        # multiply by r^2
        il1 = lmaxo2idx(2,lmaxo)
        for ir in 1:nrc
            wrcmt[ir,is] = wrcmt[ir,is] * rlcmt[ir,il1,is]
        end
        # determine the weights for partial integration on coarse radial mesh
        @views wsplintp!(nrc, rcmt[:,is], wprcmt[:,:,is])
    end

    return
end 


