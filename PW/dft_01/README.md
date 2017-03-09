# Density functional calculation

Total energy minimization using steepest descent method

```julia
psi, Energies, Potentials = kssolve_Emin_sd( pw, V_ionic, Focc, Nstates, NiterMax=1000 )
```

Total energy minimization using conjugate gradient method

```julia
psi, Energies, Potentials = kssolve_Emin_cg( pw, V_ionic, Focc, Nstates, NiterMax=1000 )
```
