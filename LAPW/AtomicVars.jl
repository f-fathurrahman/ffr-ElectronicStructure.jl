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

function AtomicVars(
    nspecies::Int64,
    natoms_::Vector{Int64},
    atposl_::Array{Float64,3},
    lattice_vars::LatticeVars
)

    epslat = lattice_vars.epslat
    avec = lattice_vars.avec

    natoms = copy(natoms_)

    # FIXME: Should not have this limitation
    maxspecies = 8 # maximum allowed species
    maxatoms = 200 # maximum allowed atoms per species

    atposl = copy(atposl_)
    atposc = zeros(3,maxatoms,maxspecies)

    for is in 1:nspecies
        for ia in 1:natoms[is]
            # map atomic lattice coordinates to [0,1)
            @views r3frac!( epslat, atposl[:,ia,is] )
            # determine atomic Cartesian coordinates
            @views r3mv!( avec,atposl[:,ia,is], atposc[:,ia,is])
        end
    end

    molecule = false
    rndatposc = 0.0
    primcell = false

    idxas = zeros(Int64, maxatoms, maxspecies)
    # inverse atoms and species indices
    idxis = zeros(Int64, maxatoms*maxspecies)
    idxia = zeros(Int64, maxatoms*maxspecies)
    natmmax = 0
    ias = 0
    for is in 1:nspecies
        for ia in 1:natoms[is]
            ias = ias + 1
            idxas[ia,is] = ias
            idxis[ias] = is
            idxia[ias] = ia
        end
        # maximum number of atoms over all species
        natmmax = max(natmmax, natoms[is])
    end
    # total number of atoms
    natmtot = ias

    return AtomicVars(
        nspecies, natoms, natmmax, natmtot,
        idxas, idxis, idxia,
        molecule, primcell, atposl, atposc, rndatposc
    )

end
