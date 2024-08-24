# Python driver
#gfortran compsep.f90 -fbacktrace -C -shared -fPIC -o libcommander.so -g
#LD_LIBRARY_PATH=. python3 main.py

# Fortran driver
gfortran compsep.f90 main.f90 -fbacktrace -O0 -C -fPIC -o main -g -fno-underscoring 
./main
