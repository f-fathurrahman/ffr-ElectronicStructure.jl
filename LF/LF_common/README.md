## File: `m_LF1d.jl`

In this file, a module named `m_LF1d` is defined.
This module describes one-dimensional LF basis.
A data type is defined here: `LF1dGrid`.

```julia
type LF1dGrid
    N::Int64
    A::Float64
    B::Float64
    h::Float64
    grid::Array{Float64,1}
    D1jl::Array{Float64,2}
    D2jl::Array{Float64,2}
end
```

`LF1dGrid` describes an LF basis defined at domain [A,B] with
spacing `h`. Number of basis functions is `N`.
The grid points are given in array `grid`.
Matrices `D1jl` and `D2jl` are required to calculate first
and second derivative operator, respectively.
Depending on the type of LF, we have several constructors
for `LF1dGrid`.

- Periodic LF:

```julia
function init_LF1d_p( N::Int64, A::Float64, B::Float64, verbose=false )
```

- Cluster (box) LF:

```julia
function init_LF1d_c( N::Int64, A::Float64, B::Float64, verbose=false )
```

- Sinc LF:

```julia
function init_LF1d_sinc( N::Int64, h::Float64, verbose=false )
```


## File: `m_LF3d.jl`

Describes LF basis functions in 3D.

Data type

```julia
type LF3dGrid
    LFx::LF1dGrid
    LFy::LF1dGrid
    LFz::LF1dGrid
    #
    Nx::Int64
    Ny::Int64
    Nz::Int64
    #
    Lx::Float64
    Ly::Float64
    Lz::Float64
    #
    lingrid::Array{Float64,2}
    xyz2lin::Array{Int64,3}
    lin2xyz::Array{Int64,2}
end
```

Constructors:

- Periodic LF:

```julia
function init_LF3d_p( NN::Array{Int64,1}, AA::Array{Float64,1}, BB::Array{Float64,1};
                      verbose=false )
```

- Cluster (box) LF:

```julia
function init_LF3d_c( NN::Array{Int64,1}, AA::Array{Float64,1}, BB::Array{Float64,1};
                      verbose=false )
```

- Sinc LF:

```julia
function init_LF3d_sinc( NN::Array{Int64,1}, hh::Array{Float64,1}; verbose=false )
```
