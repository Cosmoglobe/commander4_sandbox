g++ -shared -O3 -fPIC -fopenmp mapmaker.cpp -o mapmaker.so
mpirun -n 5 python -u -m mpi4py GibbsMPI.py