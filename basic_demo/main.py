import numpy as np
import ctypes as ct

# NB Most of the ugly details of the function calls and will be hidden from
# the "end user" later on

# open the shared library in question
fortlib = ct.CDLL('libcommander.so')

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

fortlib.data_init.argtypes           = [p_i64]
fortlib.data_init_band.argtypes      = [p_i64, p_i64, p_i64, p_dbl]
fortlib.compsep_compute_rhs.argtypes = [p_dbl, p_i64]
fortlib.compsep_compute_ax.argtypes  = [p_dbl, p_i64]
fortlib.tod_init_band.argtypes       = [p_i64, p_i64]
fortlib.tod_init_scan.argtypes       = [p_i64, p_i64, p_i64] # Need to add float+int array at the end
fortlib.tod_estimate_sigma0.argtypes = [p_i64, p_i64, p_i64] # Need to add double array between the two last
fortlib.tod_mapmaker.argtypes        = [p_i64]

ngibbs = 10
nband  = 5
nscan  = 10
ntod   = 1024
nside  = 256
npix   = 12*nside*nside
lmax   = 512
fwhm   = 0.42

# Initialize basic data
fortlib.data_init(ref_i64(nband))
for i in range(nband):
    fortlib.compsep_init_band(ref_i64(i+1), ref_i64(nside), ref_i64(lmax), ref_dbl(fwhm))

# Initialize TOD data
for i in range(nband):
    fortlib.tod_init_band(ref_i64(i+1), ref_i64(nscan))
    for j in range(nscan):
        d   = np.zeros(ntod)
        pix = np.zeros(ntod)
        for k in range(ntod):
            d[k]   = mod(k,4) + 0.6
            pix[k] = mod(k,npix)
        fortlib.tod_init_scan(ref_i64(i+1), ref_i64(j+1), ntod, d, pix)

# Run Gibbs sampler
for iter in range(ngibbs):

     print('Gibbs iter = '+iter+' of '+ngibbs)
     
     # **********************
     # COMPSEP stage
     # **********************
     
     # Compute RHS of mapmaking equation
     rhs = zeros(npix)
     fortlib.compsep_compute_rhs(ptr_dbl(rhs), ref_i64(rhs.shape[0]))
     
     # Solve for best-fit map by CG
     signal = zeros(npix)
     fortlib.compsep_compute_ax(ptr_dbl(signal), ref_i64(signal.shape[0]))
     
     # **********************
     # TOD stage
     # **********************

     # Estimate white noise rms per scan
     for i in range(nband):
        for j in range(nscan):
           fortlib.tod_estimate_sigma0(ref_i64(i+1), ref_i64(j+1), signal, ref_i64(signal.shape[0]))

     # Make frequency maps
     for i in range(nband):
        fortlib.tod_mapmaker(ref_i64(i+1))
