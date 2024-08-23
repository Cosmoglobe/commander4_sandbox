import numpy as np
import ctypes as ct

# NB Most of the ugly details of the function calls and will be hidden from
# the "end user" later on

# open the shared library in question
fortlib = ct.CDLL('demo_lib.so')
fortlib.compsep_init.argtypes = [ct.c_int64]

fortlib.compsep_init_band.argtypes = [ct.c_int64, ct.POINTER(ct.c_char), ct.c_int64, ct.c_int64, ct.c_double, ct.c_int64]

fortlib.compsep_compute_rhs.argtypes = [ct.c_int64, ct.POINTER(ct.c_double), ct.c_int64]

fortlib.compsep_compute_Ax.argtypes = [ct.c_int64, ct.POINTER(ct.c_double), ct.c_int64]

nband = 5
nside = 256
lmax = 512
fwhm = 0.42

fortlib.compsep_init(nband)

for i in range(nband):
    name = b"test_str"
    fortlib.compsep_init_band(i, name, nside, lmax, fwhm, len(name))
