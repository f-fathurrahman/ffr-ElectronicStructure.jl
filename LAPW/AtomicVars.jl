struct AtomicVars
    # number of species
    nspecies::Int64
    # number of atoms for each species
    natoms::Vector{Int64} # of length (maxspecies)
    # maximum number of atoms over all the species
    natmmax::Int64
    # total number of atoms
    natmtot::Int64
    # index to atoms and species
    idxas::Array{Int64,2} # (maxatoms,maxspecies)
    # inverse atoms and species indices
    idxis::Vector{Int64} # (maxatoms*maxspecies)
    idxia::Vector{Int64} # (maxatoms*maxspecies)
    # molecule is .true. is the system is an isolated molecule
    molecule::Bool
    # primcell is .true. if primitive unit cell is to be found automatically
    primcell::Bool
    # atomic positions in lattice coordinates
    atposl::Array{Float64,3} #(3,maxatoms,maxspecies)
    # atomic positions in Cartesian coordinates
    atposc::Array{Float64,3} # (3,maxatoms,maxspecies)
    # magnitude of random displacements added to the atomic positions
    rndatposc::Float64
end

function init_AtomicVars()
    # maximum allowed species
    maxspecies = 8
    # maximum allowed atoms per species
    maxatoms = 200
end