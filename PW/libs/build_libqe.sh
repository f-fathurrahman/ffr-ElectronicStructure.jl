gfortran -fdefault-integer-8 -O3 -c -fPIC erf.f90
gfortran -fdefault-integer-8 -O3 -c -fPIC sort.f90
gfortran -fdefault-integer-8 -O3 -c -fPIC rgen.f90
gfortran -fdefault-integer-8 -O3 -c -fPIC ewald.f90

gfortran -shared -o libqe.so erf.o sort.o rgen.o ewald.o
