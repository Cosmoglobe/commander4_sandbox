# Python driver
#gfortran libmod.f90 -fbacktrace -C -shared -fPIC -o demo_lib.so -g
#LD_LIBRARY_PATH=. python3 main.py

# Fortran driver
gfortran libmod.f90 main.f90 -fbacktrace -Og -C -fPIC -o main -g -fno-underscoring 
./main
