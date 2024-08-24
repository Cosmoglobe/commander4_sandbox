import numpy as np
import ctypes as ct

# NB Most of the ugly details of the function calls and will be hidden from
# the "end user" later on

# open the shared library in question
fortlib = ct.CDLL('demo_lib.so')

p_i64 = ct.POINTER(ct.c_int64)
p_dbl = ct.POINTER(ct.c_double)

def ref_i64(v):
    return ct.byref(ct.c_int64(v))

def ref_dbl(v):
    return ct.byref(ct.c_double(v))

# returns a pointer to a double array, for passing to Fortran
def ptr_dbl(arr):
    if arr.dtype != np.float64:
        raise RuntimeError("array dtype must be float64")
    if not arr.flags.f_contiguous:
        raise RuntimeError("bad array ordering")
    return arr.ctypes.data_as(p_dbl)

fortlib.compsep_init.argtypes = [p_i64]

fortlib.compsep_init_band.argtypes = [p_i64, p_i64, p_i64, p_dbl]

fortlib.compsep_compute_rhs.argtypes = [p_dbl, p_i64]

fortlib.compsep_compute_ax.argtypes = [p_dbl, p_i64]

nband = 5
nside = 256
lmax = 512
fwhm = 0.42

fortlib.compsep_init(ref_i64(nband))

# Initialize data
for i in range(nband):
    fortlib.compsep_init_band(ref_i64(i+1), ref_i64(nside), ref_i64(lmax), ref_dbl(fwhm))

# Compute RHS of mapmaking equation
rhs = np.zeros(12*nside**2)
print(rhs.dtype)
fortlib.compsep_compute_rhs(ptr_dbl(rhs), ref_i64(rhs.shape[0]))

# Solve for best-fit map by CG
x = np.empty(12*nside**2)
fortlib.compsep_compute_ax(ptr_dbl(x), ref_i64(x.shape[0]))
