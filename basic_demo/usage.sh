gfortran lib.f90 -O3 -shared -o demo_lib.so -g
LD_LIBRARY_PATH=. python3 main.py
