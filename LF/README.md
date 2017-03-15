Electronic structure calculations using Lagrange-basis set.

- Directory `poisson_3d_p` contains test programs for solving
  Poisson equation in periodic system. Currently, I am using
  FFT because this is much faster than CG and related methods.

- Directory `sch_3d` contains test programs for solving Schrodinger equation
  in 3D. Three types of LF basis functions are used: periodic, cluster,
  and sinc.

- Directory `extlibs`
  For preconditoning, I am using ILU preconditioner from MKL.
