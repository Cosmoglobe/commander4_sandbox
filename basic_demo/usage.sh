# Fortran driver
gfortran data_mod.f90 compsep.f90 tod.f90 main.f90 -fbacktrace -O0 -C -fPIC -o main -g
./main

# Python driver
gfortran data_mod.f90 compsep.f90 tod.f90 -fbacktrace -O0 -C -fPIC -shared -o libcommander.so -g
LD_LIBRARY_PATH=. python3 main.py
