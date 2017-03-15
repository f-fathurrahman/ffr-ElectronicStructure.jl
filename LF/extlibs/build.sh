gfortran -O3 -fPIC -shared m_Laplace3d_sparse.f90 -o Laplace3d_sparse.so
ifort -O3 -fPIC -shared m_Laplace3d_sparse.f90 -o Laplace3d_sparse_intel.so
