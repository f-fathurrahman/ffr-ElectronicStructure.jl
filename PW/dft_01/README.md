# Density functional calculation

In this directory, simple DFT calculation is performed for a electronic system
with harmonic potential and hydrogen atom. The main script can be found in file
`main_harm.jl` and `main_H_atom.jl` for harmonic potential and hydrogen atom,
respectively. The programs in this directory are using `PWGrid_v01.jl`.

## Main program flow

We will illutrate the steps taken in the program for the case of
harmonic potential.

1. Define the simulation box and generate plane wave basis for this box.

```julia
Ns = [ns1,ns2,ns3]
const LatVecs = diagm( [6.0, 6.0, 6.0] )
pw = PWGrid( Ns, LatVecs )
```

2. In the second step we setup the potential. The potential is centered at
   the simulation box and the potentials are calculated relative to this
   point (center of the box).
   We calculate the center of the box and generate an array of distance
   between each (real space) grid points and the center of the box.

```julia
center = sum(LatVecs,2)/2
dr = gen_dr( r, center )
```

   After this, we can call `init_pot_harm_3d` to initialize the potential

```julia
V_ionic = init_pot_harm_3d( pw, dr )
```

3. One the potential is set, we can solve the minimize the Kohn-Sham energy functional
   using several minimization algorithms. Before that, we need to set number of
   Kohn-Sham states and occupation number for these states.

```julia
   const Nstates = 4
   Focc = 2.0*ones(Nstates)
```

   Total energy minimization using steepest descent method

```julia
psi, Energies, Potentials = kssolve_Emin_sd( pw, V_ionic, Focc, Nstates, NiterMax=1000 )
```

   Total energy minimization using conjugate gradient method

```julia
psi, Energies, Potentials = kssolve_Emin_cg( pw, V_ionic, Focc, Nstates, NiterMax=1000 )
```

## Steepest descent minization


## Conjugate gradient minimization
