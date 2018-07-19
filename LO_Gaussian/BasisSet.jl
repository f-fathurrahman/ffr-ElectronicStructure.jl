const BasisSet = Array{CGBF,1}

include("read_basis_g94.jl")
import .g94basis

function build_basis_new( atoms::Atoms )
    
    data_dict = g94basis.read_basis_set("g94_sto-3g")

    basis = BasisSet()
    
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs

    Nbasis = 0
    for ia = 1:Natoms
        #
        x = atoms.positions[1,ia]
        y = atoms.positions[2,ia]
        z = atoms.positions[3,ia]
        #
        data = data_dict[atsymbs[ia]]
        Nsyms = data.Nsyms
        #        
        for isym = 1:Nsyms
            sym = data.syms[isym]
            for (I,J,K) in sym2power[sym]
                cbf = init_CGBF( (x,y,z), (I,J,K) )
                ncontr = data.Ncontr[isym]
                for i = 1:ncontr
                    expn = data.expns[isym][i]
                    coef = data.coeffs[isym][i]
                    push!(cbf,expn,coef)
                end
                push!(basis,cbf)
                Nbasis = Nbasis + 1
            end
        end
    end

    println("Nbasis = ", Nbasis)

    return basis

end

function build_basis( atoms::Atoms, name="sto3g", verbose::Bool=false )
    
    data = basis_set_data[name]
    basis = BasisSet()
    
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs
    
    Nbasis = 0
    for ia = 1:Natoms
        #
        atno = ZATOMS[atsymbs[ia]]
        x = atoms.positions[1,ia]
        y = atoms.positions[2,ia]
        z = atoms.positions[3,ia]
        #
        for btuple in data[atno]
            sym, primlist = btuple
            for (I,J,K) in sym2power[sym]
                cbf = init_CGBF( (x,y,z), (I,J,K) )
                for (expn,coef) in primlist
                    push!(cbf,expn,coef)
                end
                push!(basis,cbf)
                Nbasis = Nbasis + 1
            end
        end
    end

    println("Nbasis = ", Nbasis)

    return basis
end

const sym2power = Dict{Any,Any}(
  "S" => [(0,0,0)],
  "P" => [(1,0,0),(0,1,0),(0,0,1)],
  "D" => [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)],
  "F" => [(3,0,0),(0,3,0),(0,0,3),(2,1,0),(2,0,1),(1,2,0),
          (1,1,1),(1,0,2),(0,2,1),(0,1,2)]
  )

