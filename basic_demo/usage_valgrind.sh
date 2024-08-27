# Demo for running executables with valgrind
# To get good diagnostics, use "-g" for compilation.
# Valgrind controlled runs tend to be *slow* ... you'll need patience.

# Fortran driver
gfortran data_mod.f90 compsep.f90 tod.f90 main.f90 -fbacktrace -O3 -o main -g
valgrind ./main

# Python driver
# NB Python will cause a few valgrind warnings on startup. Just make sure
# there aren't any warnings after the code proper has started to run.
gfortran data_mod.f90 compsep.f90 tod.f90 -fbacktrace -O3 -C -fPIC -shared -o libcommander.so -g
LD_LIBRARY_PATH=. valgrind python3 main.py
