include("r3frac.jl")
include("r3mv.jl")
include("LatticeVars.jl")
include("AtomicVars.jl")
include("AtomicSpeciesVars.jl")
include("MuffinTinRadialVars.jl")
include("readspecies.jl")

function main()
    Nspecies = 2
    atsp_vars = AtomicSpeciesVars(Nspecies)
    mtr_vars = MuffinTinRadialVars(Nspecies)

    readspecies!(1, "DATA_species/Si.in", atsp_vars, mtr_vars)
    readspecies!(2, "DATA_species/Pt.in", atsp_vars, mtr_vars)

    println(atsp_vars.spsymb)
    println(atsp_vars.nstsp)
end

main()
