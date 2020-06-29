
function main()
    atomic_species_vars = AtomicSpeciesVars(2)
    
    println(atomic_species_vars.spsymb) # before reading information from species file

    readspecies!(1, atomic_species_vars, "DATA_species/Si.in")
    readspecies!(2, atomic_species_vars, "DATA_species/Pt.in")

    println(atomic_species_vars.spsymb)
    println(atomic_species_vars.nstsp)
end

main()
