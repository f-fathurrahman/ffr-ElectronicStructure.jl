mutable struct BasisSet # list of CGBFs
    bfs::Array{CGBF,1}
end

basisset() = BasisSet(CGBF[])

function push!(basis::BasisSet,cbf::CGBF)
    Base.push!(basis.bfs,cbf)
end

function build_basis( atoms::Atoms,name="sto3g" )
    data = basis_set_data[name]
    basis_set = basisset()
    
    for atom in atoms.atomlist
        
        #atno = get_Zatoms(atoms)
        
        for btuple in data[atoms.atno]
            sym,primlist = btuple
            for (I,J,K) in sym2power[sym]
                cbf = init_CGBF( (atom.x,atom.y,atom.z), (I,J,K) )
                push!(basis_set,cbf)
                for (expn,coef) in primlist
                    push!(cbf,expn,coef)
                end
            end
        end

    end

    return basis_set
end

const sym2power = Dict{Any,Any}(
  'S' => [(0,0,0)],
  'P' => [(1,0,0),(0,1,0),(0,0,1)],
  'D' => [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]
  )
