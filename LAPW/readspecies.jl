function readspecies!( isp::Int64, filename,
    atsp_vars::AtomicSpeciesVars,
    mtr_vars::MuffinTinRadialVars,
    apwlo_vars::APWLOVars;
    verbose=false
)
    f = open(filename, "r")
    
    # species symbol
    line = readline(f)
    atsp_vars.spsymb[isp] = replace(split(line)[1], "'" => "")
    verbose && println("spsymb = ", atsp_vars.spsymb[isp])
    
    # species name
    line = readline(f)
    atsp_vars.spname[isp] = replace(split(line)[1], "'" => "")
    verbose && println("spname = ", atsp_vars.spname[isp])
    
    # atomic number
    line = readline(f)
    atsp_vars.spzn[isp] = parse(Float64, split(line)[1])
    verbose && println("spzn   = ", atsp_vars.spzn[isp])

    # mass
    line = readline(f)
    atsp_vars.spmass[isp] = parse(Float64, split(line)[1])
    verbose && println("spmass = ", atsp_vars.spmass[isp])

    # Radial mesh
    line = readline(f)
    ll = split(line)
    
    atsp_vars.rminsp[isp] = parse(Float64, ll[1])
    
    mtr_vars.rmt[isp] = parse(Float64, ll[2])  # muffin tin
    
    atsp_vars.rmaxsp[isp] = parse(Float64, ll[3])
    
    mtr_vars.nrmt[isp] = parse(Int64, ll[4])
    
    verbose && println("rminsp = ", atsp_vars.rminsp[isp])
    verbose && println("rmt    = ", mtr_vars.rmt[isp])
    verbose && println("rmaxsp = ", atsp_vars.rmaxsp[isp])
    verbose && println("nrmt   = ", mtr_vars.nrmt[isp])

    # Atomic states
    line = readline(f)
    atsp_vars.nstsp[isp] = parse(Int64, split(line)[1])
    #
    for ist in 1:atsp_vars.nstsp[isp]
        #
        line = readline(f)
        ll = split(line)
        #
        atsp_vars.nsp[ist,isp] = parse(Int64, ll[1])
        atsp_vars.lsp[ist,isp] = parse(Int64, ll[2])
        atsp_vars.ksp[ist,isp] = parse(Int64, ll[3]) # whats this?
        atsp_vars.occsp[ist,isp] = parse(Float64, ll[4])
        if ll[5] == "T"
            atsp_vars.spcore[ist,isp] = true
        elseif ll[5] == "F"
            atsp_vars.spcore[ist,isp] = false
        else
            error("Unable to parse spcore")
        end
    end

    # apword
    line = readline(f)
    apwlo_vars.apword[1,isp] = parse(Int64, split(line)[1])
    # set the APW orders for > 0
    lmaxapw = mtr_vars.lmaxapw
    apwlo_vars.apword[2:lmaxapw+1,isp] .= apwlo_vars.apword[1,isp]

    # apwe0, apwdm, apwve
    for iorb in apwlo_vars.apword[1,isp]
        line = readline(f)
        ll = split(line)
        apwlo_vars.apwe0[iorb,1,isp] = parse(Float64,ll[1])
        apwlo_vars.apwdm[iorb,1,isp] = parse(Int64, ll[2])
        if ll[3] == "T"
            apwlo_vars.apwve[iorb,1,isp] = true
        elseif ll[3] == "F"
            apwlo_vars.apwve[iorb,1,isp] = false
        else
            error("Unable to parse apwve")
        end
        # 
        # set the APW linearisation energies, derivative orders and variability for l > 0
        apwlo_vars.apwe0[iorb,2:lmaxapw+1,isp] .= apwlo_vars.apwe0[iorb,1,isp]
        apwlo_vars.apwdm[iorb,2:lmaxapw+1,isp] .= apwlo_vars.apwdm[iorb,1,isp]
        apwlo_vars.apwve[iorb,2:lmaxapw+1,isp] .= apwlo_vars.apwve[iorb,1,isp]
        #
        apwlo_vars.e0min = min( apwlo_vars.e0min, apwlo_vars.apwe0[iorb,1,isp] )
    end

    # nlx (skipped?) XXX
    line = readline(f)
    nlx = parse(Int64, split(line)[1])
    if nlx > 0
        error("Unsupported value of nlx")
    end

    # nlorb
    line = readline(f)
    nlorb = parse(Int64, split(line)[1])
    for iorb in 1:nlorb
        # lorbl, lorbord
        line = readline(f); ll = split(line)
        apwlo_vars.lorbl[iorb,isp] = parse(Int64, ll[1])
        apwlo_vars.lorbord[iorb,isp] = parse(Int64, ll[2])
        #
        for i in 1:apwlo_vars.lorbord[iorb,isp]
            # lorbe0, lorbdm, lorbve
            line = readline(f); ll = split(line)
            apwlo_vars.lorbe0[i,iorb,isp] = parse(Float64, ll[1])
            apwlo_vars.lorbdm[i,iorb,isp] = parse(Float64, ll[2])
            if ll[3] == "T"
                apwlo_vars.lorbve[i,iorb,isp] = true
            elseif ll[3] == "F"
                apwlo_vars.lorbve[i,iorb,isp] = false
            else
                error("Unable to parse lorbve")
            end
        end
    end

    close(f)

end