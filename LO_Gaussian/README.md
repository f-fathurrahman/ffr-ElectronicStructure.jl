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

# Basis set

Primitive Gaussian basis function:

```Julia
type PGBF
  expn::Float64
  x::Float64
  y::Float64
  z::Float64
  I::Int64
  J::Int64
  K::Int64
  norm::Float64
end
```

Contracted Gaussian basis function:

```julia
type CGBF
  x::Float64
  y::Float64
  z::Float64
  I::Int64
  J::Int64
  K::Int64
  norm::Float64
  pgbfs::Array{PGBF,1}
  coefs::Array{Float64,1}
end
```
