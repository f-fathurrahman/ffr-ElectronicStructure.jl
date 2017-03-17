type Atom
  atno::Int64
  x::Float64
  y::Float64
  z::Float64
end

type Molecule
  atomlist::Array{Atom,1}
end

function push!(mol::Molecule,at::Atom)
  Base.push!(atomlist,at)
end

tobohr(x::Float64) = x/0.52918
function tobohr!(at::Atom)
  at.x /= 0.52918
  at.y /= 0.52918
  at.z /= 0.52918
end

function tobohr!(mol::Molecule)
  for at in mol.atomlist
    tobohr!(at)
  end
end
