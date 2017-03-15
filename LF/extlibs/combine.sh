icc -I/home/efefer/intel/mkl/include -c -fPIC -O3 mkl_ilu0.c
icc -I/home/efefer/intel/mkl/include -c -fPIC -O3 mkl_dfgmres.c
icc -shared -o libmkl.so -Wl,--whole-archive \
   mkl_dfgmres.o mkl_ilu0.o \
   libmkl_core.so libmkl_intel_ilp64.so libmkl_sequential.so \
   libmkl_mc3.so libmkl_def.so \
   -Wl,--no-whole-archive

# No warnings during compilation, but does not work
#icc -DINTEL_ILP64 -I/home/efefer/intel/mkl/include -c -fPIC -O3 mkl_ilu0.c
#icc -DINTEL_ILP64 -shared -o libmkl_ilu0.so -Wl,--whole-archive mkl_ilu0.o \
#    libmkl_core.so libmkl_intel_ilp64.so libmkl_sequential.so -Wl,--no-whole-archive
