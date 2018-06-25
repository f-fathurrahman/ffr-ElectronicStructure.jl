
const BasisSet = Array{CGBF,1}

function push!(basis::BasisSet,cbf::CGBF)
    Base.push!(basis,cbf)
end

function build_basis( atoms::Atoms, name="sto3g" )
    data = basis_set_data[name]
    basis = BasisSet()
    
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs
    
    for ia = 1:Natoms

        atno = ZATOMS[atsymbs[ia]]
        x = atoms.positions[1,ia]
        y = atoms.positions[2,ia]
        z = atoms.positions[3,ia]                
        
        for btuple in data[atno]
            sym,primlist = btuple
            for (I,J,K) in sym2power[sym]
                cbf = init_CGBF( (x,y,z), (I,J,K) )
                push!(basis,cbf)
                for (expn,coef) in primlist
                    push!(cbf,expn,coef)
                end
            end
        end

    end

    return basis
end

const sym2power = Dict{Any,Any}(
  "S" => [(0,0,0)],
  "P" => [(1,0,0),(0,1,0),(0,0,1)],
  "D" => [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]
  )

