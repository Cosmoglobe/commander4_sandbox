import numpy as np
import ctypes as ct
import healpy as hp
import h5py
from tqdm import tqdm

# NB Most of the ugly details of the function calls and will be hidden from
# the "end user" later on

# open the shared library in question
fortlib = ct.CDLL('libcommander.so')

print("START")

i64 = ct.c_int64
i32 = ct.c_int32
dbl = ct.c_double
flt = ct.c_float
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
fortlib.tod_mapmaker_ifc.argtypes        = [i64, a_f64_1, a_f64_1, i64]

ngibbs = 5
nband  = 5
nscan  = 108
ntod   = 2**16
nside  = 256
npix   = 12*nside*nside
lmax   = 512
fwhm   = 0.42


# Initialize basic data
fortlib.data_init_ifc(nband)
for i in range(nband):
    fortlib.data_init_band_ifc(i+1, nside, lmax, fwhm)


h5file = '../src/python/preproc_scripts/tod_example.h5'
bands = ['030', '070', '100', '217', '353']
# Initialize TOD data
with h5py.File(h5file) as data:
    nscan = len(data.keys())-1
    for i in range(nband):
        fortlib.tod_init_band_ifc(i+1, nscan)
        for j in range(nscan):
            d = data[f'{j+1:06}/{bands[i]}/tod'][()]
            pix = data[f'{j+1:06}/{bands[i]}/pix'][()]
            ntod = data[f'{j+1:06}/common/ntod'][0]
            fortlib.tod_init_scan_ifc(i+1, j+1, ntod, d, pix)

# Run Gibbs sampler
for iter in range(1,ngibbs+1):

     print(f'Gibbs iter = {iter} of {ngibbs}')
     
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
        m_i   = np.zeros(12*nside**2, dtype=dbl)
        rms_i = np.zeros(12*nside**2, dtype=dbl)
        fortlib.tod_mapmaker_ifc(i+1, m_i, rms_i, npix)
        hp.write_map(f'output/map_band_{i:02}_c{iter:06}.fits', m_i,
                overwrite=True, dtype=dbl)
        hp.write_map(f'output/rms_band_{i:02}_c{iter:06}.fits', rms_i,
                overwrite=True, dtype=dbl)
