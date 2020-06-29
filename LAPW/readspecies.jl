function readspecies!(isp::Int64, vars::AtomicSpeciesVars, filename)
    f = open(filename, "r")
    
    # species symbol
    line = readline(f)
    vars.spsymb[isp] = replace(split(line)[1], "'" => "")
    println("spsymb = ", vars.spsymb[isp])
    
    # species name
    line = readline(f)
    vars.spname[isp] = replace(split(line)[1], "'" => "")
    println("spname = ", vars.spname[isp])
    
    # atomic number
    line = readline(f)
    vars.spzn[isp] = parse(Float64, split(line)[1])
    println("spzn   = ", vars.spzn[isp])

    # mass
    line = readline(f)
    vars.spmass[isp] = parse(Float64, split(line)[1])
    println("spmass = ", vars.spmass[isp])

    # Radial mesh
    line = readline(f)
    ll = split(line)
    vars.rminsp[isp] = parse(Float64, ll[1])
    rmt = parse(Float64, ll[2])  # muffin tin
    vars.rmaxsp[isp] = parse(Float64, ll[3])
    nrmt = parse(Int64, ll[4])
    println("rminsp = ", vars.rminsp[isp])
    println("rmt    = ", rmt)
    println("rmaxsp = ", vars.rmaxsp[isp])
    println("nrmt   = ", nrmt)

    # Atomic states
    line = readline(f)
    vars.nstsp[isp] = parse(Int64, split(line)[1])
    for ist in 1:vars.nstsp[isp]
        line = readline(f)
        ll = split(line)
        #
        vars.nsp[ist,isp] = parse(Int64, ll[1])
        vars.lsp[ist,isp] = parse(Int64, ll[2])
        vars.ksp[ist,isp] = parse(Int64, ll[3]) # whats this?
        vars.occsp[ist,isp] = parse(Float64, ll[4])
        if ll[5] == "T"
            vars.spcore[ist,isp] = true
        elseif ll[5] == "F"
            vars.spcore[ist,isp] = false
        else
            error("Unable to parse spcore")
        end
    end

    close(f)

end