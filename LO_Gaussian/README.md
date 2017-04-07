These came from Rick Mueller (PyQuante's author)

I modified it slightly to remove errors and warnings.

# Defining molecular structures

Molecule is described by simple arrays of Atom:

```julia
type Atom
  atno::Int64
  x::Float64
  y::Float64
  z::Float64
end

type Molecule
  atomlist::Array{Atom,1}
end
```

Probably `Molecule` will be renamed to `Atoms`.

No convenience constructors/functions are provided to
construct an instance of `Molecule`.
