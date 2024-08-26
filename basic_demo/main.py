import numpy as np
import ctypes as ct

# NB Most of the ugly details of the function calls and will be hidden from
# the "end user" later on

# open the shared library in question
fortlib = ct.CDLL('libcommander.so')

print("START")

i64 = ct.c_int64
i32 = ct.c_int32
dbl = ct.c_double
flt = ct.c_float
p_i32 = ct.POINTER(ct.c_int32)
p_i64 = ct.POINTER(ct.c_int64)
p_flt = ct.POINTER(ct.c_float)
p_dbl = ct.POINTER(ct.c_double)
a_f32_1 = np.ctypeslib.ndpointer(dtype=ct.c_float, ndim=1, flags="F")
a_f64_1 = np.ctypeslib.ndpointer(dtype=ct.c_double, ndim=1, flags="F")
a_i32_1 = np.ctypeslib.ndpointer(dtype=ct.c_int32, ndim=1, flags="F")

def ref_i64(v):
    return ct.byref(ct.c_int64(v))

def ref_dbl(v):
    return ct.byref(ct.c_double(v))

fortlib.data_init_ifc.argtypes           = [i64]
fortlib.data_init_band_ifc.argtypes      = [i64, i64, i64, dbl]
fortlib.compsep_compute_rhs_ifc.argtypes = [a_f64_1, i64]
fortlib.compsep_compute_Ax_ifc.argtypes  = [a_f64_1, i64]
fortlib.tod_init_band_ifc.argtypes       = [i64, i64]
fortlib.tod_init_scan_ifc.argtypes       = [i64, i64, i64, a_f32_1, a_i32_1]
fortlib.tod_estimate_sigma0_ifc.argtypes = [i64, i64, a_f64_1, i64]
fortlib.tod_mapmaker_ifc.argtypes        = [i64]

ngibbs = 10
nband  = 5
nscan  = 10
ntod   = 1024
nside  = 256
npix   = 12*nside*nside
lmax   = 512
fwhm   = 0.42

# Initialize basic data
fortlib.data_init_ifc(nband)
for i in range(nband):
    fortlib.data_init_band_ifc(i+1, nside, lmax, fwhm)

# Initialize TOD data
for i in range(nband):
    fortlib.tod_init_band_ifc(i+1, nscan)
    for j in range(nscan):
        d   = np.zeros(ntod, dtype=np.float32)
        pix = np.zeros(ntod, dtype=np.int32)
        for k in range(ntod):
            d[k]   = k%4 + 0.6
            pix[k] = k%npix
        fortlib.tod_init_scan_ifc(i+1, j+1, ntod, d, pix)

# Run Gibbs sampler
for iter in range(ngibbs):

     print('Gibbs iter = '+str(iter)+' of '+str(ngibbs))
     
     # **********************
     # COMPSEP stage
     # **********************
     
     # Compute RHS of mapmaking equation
     rhs = np.zeros(npix)
     fortlib.compsep_compute_rhs_ifc(rhs, rhs.shape[0])
     
     # Solve for best-fit map by CG
     signal = np.zeros(npix)
     fortlib.compsep_compute_Ax_ifc(signal, signal.shape[0])
     
     # **********************
     # TOD stage
     # **********************

     # Estimate white noise rms per scan
     for i in range(nband):
        for j in range(nscan):
           fortlib.tod_estimate_sigma0_ifc(i+1, j+1, signal, signal.shape[0])

     # Make frequency maps
     for i in range(nband):
        fortlib.tod_mapmaker_ifc(i+1)
