mutable struct AtomicSpeciesVars
    # species files path
    sppath::String
    # species filenames
    spfname::Vector{String}
    # species name
    spname::Vector{String}
    # species symbol
    spsymb::Vector{String}
    # species nuclear charge
    spzn::Vector{Float64}
    # ptnucl is .true. if the nuclei are to be treated as point charges, if .false.
    # the nuclei have a finite spherical distribution
    ptnucl::Bool
    # nuclear radius
    rnucl::Vector{Float64}
    # nuclear volume
    volnucl::Vector{Float64}
    # number of radial mesh points to nuclear radius
    nrnucl::Vector{Int64}
    # number of coarse radial mesh points to nuclear radius
    nrcnucl::Vector{Int64}
    # nuclear Coulomb potential
    vcln::Array{Float64,2}
    # species electronic charge
    spze::Vector{Float64}
    # species mass
    spmass::Vector{Float64}
    # smallest radial point for each species
    rminsp::Vector{Float64}
    # effective infinity for species
    rmaxsp::Vector{Float64}
    # number of radial points to effective infinity for each species
    nrsp::Vector{Int64}
    # maximum nrsp over all the species
    nrspmax::Int64
    # maximum allowed states for each species
    maxstsp::Int64 # parameter: 40
    # number of states for each species
    nstsp::Vector{Int64}
    # maximum nstsp over all the species
    nstspmax::Int64
    # core-valence cut-off energy for species file generation
    ecvcut::Float64
    # semi-core-valence cut-off energy for species file generation
    esccut::Float64
    # state principle quantum number for each species
    nsp::Array{Int64,2} # (maxstsp,maxspecies)
    # state l value for each species
    lsp::Array{Int64,2} # (maxstsp,maxspecies)
    # state k value for each species
    ksp::Array{Int64,2} #(maxstsp,maxspecies)
    # spcore is .true. if species state is core
    spcore::Array{Bool,2} #(maxstsp,maxspecies)
    # total number of core states
    nstcr::Int64
    # state eigenvalue for each species
    evalsp::Array{Float64,2} # (maxstsp,maxspecies)
    # state occupancy for each species
    occsp::Array{Float64,2} # (maxstsp,maxspecies)
    # species radial mesh to effective infinity
    rsp::Array{Float64,2}
    # r^l on radial mesh to muffin-tin radius
    rlsp::Array{Float64,3}
    # species charge density
    rhosp::Array{Float64,2}
    # species self-consistent potential
    vrsp::Array{Float64,2}
    # exchange-correlation type for atomic species (the converged ground-state of
    # the crystal does not depend on this choice)
    xctsp::Tuple{Int64,Int64,Int64}
end


function AtomicSpeciesVars( Nspecies::Int64 )
    
    sppath = ""
    spfname = Vector{String}(undef, Nspecies)
    spname = Vector{String}(undef, Nspecies)
    spsymb = Vector{String}(undef, Nspecies)
    
    spzn = zeros(Float64, Nspecies)    

    ptnucl = true

    rnucl = zeros(Float64, Nspecies)

    volnucl = zeros(Float64, Nspecies)

    nrnucl = zeros(Int64, Nspecies)
    nrcnucl = zeros(Int64, Nspecies)
    
    vcln = zeros(Float64,1,1) # XXX will be properly setup later
    
    spze = zeros(Float64, Nspecies)
    spmass = zeros(Float64, Nspecies)
    
    rminsp = zeros(Float64, Nspecies)
    rmaxsp = zeros(Float64, Nspecies)
    nrsp = zeros(Int64, Nspecies)
    nrspmax = 0
    maxstsp = 40
    nstsp = zeros(Int64, Nspecies)
    nstspmax = 0
    
    ecvcut = 0.0
    esccut = 0.0
    
    nsp = zeros(Int64,maxstsp,Nspecies)
    lsp = zeros(Int64,maxstsp,Nspecies)
    ksp = zeros(Int64,maxstsp,Nspecies)
    spcore = zeros(Bool,maxstsp,Nspecies)
    nstcr = 0
    evalsp = zeros(Float64,maxstsp,Nspecies)
    occsp  = zeros(Float64,maxstsp,Nspecies)

    # XXX Will be setup properly later
    rsp = zeros(Float64,1,1)
    rlsp  = zeros(Float64,1,1,1) 
    rhosp = zeros(Float64,1,1)
    vrsp  = zeros(Float64,1,1)
    xctsp  = (3,0,0)

    return AtomicSpeciesVars(
        sppath, spfname, spname, spsymb, spzn,
        ptnucl, rnucl, volnucl, nrnucl, nrcnucl,
        vcln, spze, spmass, rminsp, rmaxsp, nrsp,
        nrspmax, maxstsp, nstsp, nstspmax, ecvcut, esccut,
        nsp, lsp, ksp, spcore, nstcr, evalsp, occsp, rsp, rlsp, rhosp, vrsp, xctsp
    )
end

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

function main()
    atomic_species_vars = AtomicSpeciesVars(2)
    
    println(atomic_species_vars.spsymb) # before reading information from species file

    readspecies!(1, atomic_species_vars, "DATA_species/Si.in")
    readspecies!(2, atomic_species_vars, "DATA_species/Pt.in")

    println(atomic_species_vars.spsymb)
    println(atomic_species_vars.nstsp)
end

main()
