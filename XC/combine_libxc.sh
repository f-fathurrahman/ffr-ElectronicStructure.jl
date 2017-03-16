rm -f ffrlibxc.o
gcc -I ~/mysoftwares/libxc-3.0.0/include -c -fPIC -O3 ffrlibxc.c
gcc -shared -o ffrlibxc.so -Wl,--whole-archive \
   ffrlibxc.o libxc.so.4 \
   -Wl,--no-whole-archive
