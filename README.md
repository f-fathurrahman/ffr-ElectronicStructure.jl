# ffr-ElectronicStructure.jl

Simple electronic structure calculations implemented in
Julia programming language.

Many of the codes here are inspired from [Prof. Tomas Arias' Practical DFT course](http://jdftx.org/PracticalDFT.html)
especially the example problem for Poisson equation and minimization algorithms.

Although written in Julia, I have tried to avoid using advanced language
features of Julia. I mainly used Julia for rapid prototyping.
I have chosen to make the Julia code easy to port to Fortran, which is the
main language I used for the implementation.

I found Julia is very pleasing to work with, as compared to MATLAB, Octave or
Python (Numpy).
The main reason is that the loop structure can be made very similar to
Fortran without having to worry about worse performance. In MATLAB or Octave
I have to avoid this loop and find a way to vectorize the loop operation.


## TODO

- Numerical methods:
  - minimization (SD and CG)
  - iterative diagonalization: Lanczos, band-by-band CG, Davidson, and LOBPCG

- Longer term: Finite difference, Lagrange basis, and Gaussian basis
