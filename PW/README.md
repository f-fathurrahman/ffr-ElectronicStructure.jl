Electronic structure calculation using plane wave basis set.


- Directory `common` contains several definitions of type `PWGrid` and
  various commonly used functions.

- Directory `pwgrid_01` contains a test program for `PWGrid_v01.jl`.

- Directory `poisson_01` contains a test program that solve Poisson equation
  for a given charge density. It uses `PWGrid_v01.jl`.

- Directory `pwgrid_02` contains several programs that test the implementation
  of structure factor calculation. The implementation is tested by calculating
  atom-centered quantity and visualizing it in 3D.

- Directory `pwgrid_03` contains several a program to test `PWGrid_v02.jl`.
  In this test, there are three main types: `PWGrid`, `GVectors` and `GVectorsW`.
  Note that `PWGrid_v02.jl` and `PWGrid_v03.jl` is quite similar.

- Directory `pwgrid_04` contains a test program for `PWGrid_v04.jl` which includes
  k-point sampling. The k-point list is supplied by the user.
  FIXME: The file `PWGrid_v04.jl`
  is not yet pulled up to directory `common`.

- Directory `ewald_01` contains several test programs to test the implementation
  for calculating Ewald energy. It uses `PWGrid_v01.jl`.

- Directory `ewald_02` contains several test programs to test the implementation
  of Ewald energy calculator for atoms with Z /= 1. The results are compared to
  PWSCF results.

- Directory `sch_01` contains several test programs to solve Schrodinger equation
  using several methods: energy minimization methods (SD and PCG) and diagonalization
  (LOBPCG). It uses `PWGrid_v01.jl`.

- Directory `sch_02` contains similar test programs as `sch_01` does. The main
  difference is that `sch_02` uses `PWGrid_v02.jl` and `PWGrid_v03.jl`.

- Directory `sch_03` several test programs that solve Schrodinger equation with
  k-points using iterative diagonalization method. *Energy minimization method is
  not yet tested*.

- Directory `dft_01` contains several test programs for simple DFT calculations.
  It uses `PWGrid_v01.jl`.

- Directory `dft_02` contains similar tests as `dft_01`, but uses `PWGrid_v02.jl`
  and `PWGrid_v03.jl` instead.
