Electronic structure calculations using Lagrange-basis set.

- Directory `poisson_3d_p` contains test programs for solving
  Poisson equation in periodic system. Currently, I am using
  FFT because this is much faster than CG and related methods.

- Directory `sch_3d` contains test programs for solving Schrodinger equation
  in 3D. Three types of LF basis functions are used: periodic, cluster,
  and sinc.

- Directory `dft_3d_p_v1` contains test programs for solving Kohn-Sham
  equations in periodic systems.

- Directory `extlibs` contains various source files in C and Fortran.
  For preconditoning, ILU0 preconditioner from MKL is currently used, so
  MKL is needed to build the file `libmkl.so`. This is simply a wrapper
  for calling MKL routines to build and apply ILU0 preconditioner.
