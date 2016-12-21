# Solving Poisson equation

This is an example of using FFT to solve Poisson equation, i.e.
```
\nabla^2 V_Ha(r) = -4*pi*rho(r)
```

In G-space, Poisson equation becomes
```
-G^2 V_Ha(G) = -4*pi*rho(G)
```

So, Hartree potential can be calculated as
```
V_Ha(G) = 
```
excluding the G /= 0 term.

Hartree energy is calculated as
```
E_Ha = \int V_Ha(r) * rho(r) * dr
```


## Program flow

First, number of sampling points and box dimension is defined
and an instance of `PWGrid` is initialized.

The charge density `rho` is constructed in real space first, and
then transformed to reciprocal space using FFT. This charge density
is difference between two normalized Gaussian functions with different
`sigma` parameters centered at the calculation box.

The numerical result is then compared to analytical result.


