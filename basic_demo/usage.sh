gfortran libmod.f90 -O3 -shared -fPIC -o demo_lib.so -g
LD_LIBRARY_PATH=. python3 main.py
