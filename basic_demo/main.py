import numpy as np
import ctypes as ct

# NB Most of the ugly details of the function calls and will be hidden from
# the "end user" later on

# open the shared library in question
fortlib = ct.CDLL('libcommander.so')
fortlib.compsep_init.argtypes = [ct.c_int64]
fortlib.compsep_init_band.argtypes = [ct.c_int64, ct.c_int64, ct.c_int64, ct.c_double]
fortlib.compsep_compute_rhs.argtypes = [ct.POINTER(ct.c_double), ct.c_int64]
fortlib.compsep_compute_ax.argtypes = [ct.POINTER(ct.c_double), ct.c_int64]

nband = 5
nside = 256
lmax = 512
fwhm = 0.42

fortlib.compsep_init(nband)

# Initialize data
for i in range(nband):
    fortlib.compsep_init_band(i+1, nside, lmax, fwhm)

# Compute RHS of mapmaking equation
rhs = np.zeros(12*nside**2)
print(rhs.dtype)
fortlib.compsep_compute_rhs(rhs.ctypes.data_as(ct.POINTER(ct.c_double)), rhs.shape[0])

# Solve for best-fit map by CG
x = np.empty(12*nside**2)
fortlib.compsep_compute_ax(x.ctypes.data_as(ct.POINTER(ct.c_double)), x.shape[0])
