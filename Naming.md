An attempt to make naming scheme more consistent

## Functions

```julia
#
# Hamiltonian operator, apply Hamiltonian to input vectors
#
op_H()

#
# Kinetic operator
#
op_K()

#
# Local pseudopotential operator
# can also represent other non-atomic local potential
#
op_V_ps_loc()

#
# Calculate structure factor
#
calc_strfact()

#
# Calculate total energy components
#
calc_Energies()

#
# Calculate electronic densities
#
calc_rho
# or
calc_rhoe # ?

#
# Calculating Ewald energy
#
calc_ewald()

#
# Solving Poisson equation
#
solve_poisson

#
# Solving KS equation or minimize KS energy functional
#
kssolve_scf
kssolve_Emin_cg
# or ?
solve_ks_scf
solve_ks_Emin_cg

#
# Fourier transform
#
G_to_R()
#
R_to_G()
```

## Variables and data structures

```
PWGrid

PotentialsT

EnergiesT
```
