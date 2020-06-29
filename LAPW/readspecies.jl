function readspecies!( isp::Int64, filename,
    atsp_vars::AtomicSpeciesVars,
    mtr_vars::MuffinTinRadialVars
)
    f = open(filename, "r")
    
    # species symbol
    line = readline(f)
    atsp_vars.spsymb[isp] = replace(split(line)[1], "'" => "")
    println("spsymb = ", atsp_vars.spsymb[isp])
    
    # species name
    line = readline(f)
    atsp_vars.spname[isp] = replace(split(line)[1], "'" => "")
    println("spname = ", atsp_vars.spname[isp])
    
    # atomic number
    line = readline(f)
    atsp_vars.spzn[isp] = parse(Float64, split(line)[1])
    println("spzn   = ", atsp_vars.spzn[isp])

    # mass
    line = readline(f)
    atsp_vars.spmass[isp] = parse(Float64, split(line)[1])
    println("spmass = ", atsp_vars.spmass[isp])

    # Radial mesh
    line = readline(f)
    ll = split(line)
    
    atsp_vars.rminsp[isp] = parse(Float64, ll[1])
    
    mtr_vars.rmt[isp] = parse(Float64, ll[2])  # muffin tin
    
    atsp_vars.rmaxsp[isp] = parse(Float64, ll[3])
    
    mtr_vars.nrmt[isp] = parse(Int64, ll[4])
    
    println("rminsp = ", atsp_vars.rminsp[isp])
    println("rmt    = ", mtr_vars.rmt[isp])
    println("rmaxsp = ", atsp_vars.rmaxsp[isp])
    println("nrmt   = ", mtr_vars.nrmt[isp])

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

    close(f)

end