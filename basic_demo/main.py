import numpy as np
import ctypes as ct

# NB Most of the ugly details of the function calls and will be hidden from
# the "end user" later on

# open the shared library in question
fortlib = ct.CDLL('demo_lib.so')
fortlib.compsep_init.argtypes = [ct.c_int64]

fortlib.compsep_init_band.argtypes = [ct.c_int64, ct.c_char_p, ct.c_int64, ct.c_int64, ct.c_int64, ct.c_double]

fortlib.compsep_compute_rhs.argtypes = [ct.c_int64, ct.POINTER(ct.c_double), ct.c_int64]

fortlib.compsep_compute_Ax.argtypes = [ct.c_int64, ct.POINTER(ct.c_double), ct.c_int64]

nband = 5
nside = 256
lmax = 512
fwhm = 0.42

fortlib.compsep_init(nband)

for i in range(nband):
    name = b"test_str"
    buf = ct.create_string_buffer(name)
    fortlib.compsep_init_band(i+1, buf, len(name), nside, lmax, fwhm)
    rhs = np.zeros(12*nside**2)
    print(rhs.dtype)
    fortlib.compsep_compute_rhs(i+1, rhs.ctypes.data_as(ct.POINTER(ct.c_double)), rhs.shape[0])
    x = np.empty(12*nside**2)
    fortlib.compsep_compute_Ax(i+1, x.ctypes.data_as(ct.POINTER(ct.c_double)), x.shape[0])
