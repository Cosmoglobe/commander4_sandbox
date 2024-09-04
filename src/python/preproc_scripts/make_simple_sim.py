# TOD simulation script

# Data = CMB + noise
#   Nside  = 2048
#   Lmax   = 6000
#   FWHM   = 10 arcmin
#   sigma0 = 30µK
#
# Scanning strategy = visit each pixel in order; repeat 9 times, such that final noise is 10µK/pix
#
# Split in chunks with 2^22 samples each (except for last one) = ~109 files, total of ~4GB
# 

import numpy as np
import healpy as hp
import astropy.units as u
from commander_tod import commander_tod
import matplotlib.pyplot as plt

import camb
from camb import model, initialpower


nside = 512
lmax = 3*nside
fwhm = 10*u.arcmin
sigma0s = np.array([100, 80, 30, 150, 220])/1e3
freqs = [30, 70, 100, 217, 353]


pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06,
                       As=2e-9, ns=0.965, halofit_version='mead', lmax=lmax)

results = camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars, CMB_unit='muK', raw_cl=True)
totCL=powers['total']

ell = np.arange(lmax+1)
Cl = totCL[ell,0]
Cl_EE = totCL[ell,1]
Cl_BB = totCL[ell,2]
Cl_TE = totCL[ell,3]

Cls = np.array([Cl, Cl_EE, Cl_BB, Cl_TE])




chunk_size = 2**16

np.random.seed(0)
alms = hp.synalm(Cls, lmax=3*nside-1, new=True)
ms   = hp.alm2map(alms, nside, pixwin=False)
hp.write_map("true_sky.fits", ms, overwrite=True)

ms = hp.smoothing(ms, fwhm=fwhm.to('rad').value)

npix = 12*nside**2

ntod = 9*npix

pix = np.arange(ntod) % npix
psi = np.repeat(np.arange(9)*np.pi/9, npix)

T,Q,U = ms
d = T[pix] + Q[pix]*np.cos(2*psi) + U[pix]*np.sin(2*psi)
ds = []
for i in range(len(freqs)):
    ds.append((d + np.random.randn(ntod)*sigma0s[i]).astype('float32'))

pix = pix.astype('int32')


n_chunks = ntod // chunk_size
print(f'Number of scans is {n_chunks}')



output_path = '.'
version = 'hello'
comm_tod = commander_tod(output_path, "", version, overwrite=True)

hdf_filename = 'tod_example'

COMMON_GROUP = "/common"
HUFFMAN_COMPRESSION = ["huffman", {"dictNum": 1}]

comm_tod.init_file(freq=hdf_filename, od="", mode="w")
#
## nside
comm_tod.add_field(COMMON_GROUP + "/nside", [nside])

for pid in range(n_chunks):
    pid_label = f'{pid+1:06}'
    pid_common_group = pid_label + "/common"
    for i, freq in enumerate(freqs):
        pid_data_group = f'{pid_label}/{freq:03}'

        comm_tod.add_field(pid_common_group + "/ntod", [chunk_size])


        tod_chunk_i = ds[i][pid*chunk_size : (pid+1)*chunk_size]
        pix_chunk_i =   pix[pid*chunk_size : (pid+1)*chunk_size]
        psi_chunk_i =   psi[pid*chunk_size : (pid+1)*chunk_size]
        
        comm_tod.add_field(pid_data_group + "/tod", tod_chunk_i)
        comm_tod.add_field(pid_data_group + "/pix", pix_chunk_i)
        comm_tod.add_field(pid_data_group + "/psi", psi_chunk_i)
    comm_tod.finalize_chunk(pid+1)

if (ntod//chunk_size != ntod/chunk_size):
    pid = n_chunks
    pid_label = f'{pid+1:06}'
    pid_common_group = pid_label + "/common"
    for i, freq in enumerate(freqs):
        pid_data_group = f'{pid_label}/{freq:03}'
    
        tod_chunk_i = ds[i][pid*chunk_size : ]
        pix_chunk_i = pix[pid*chunk_size : ]
        psi_chunk_i = psi[pid*chunk_size : ]
        comm_tod.add_field(pid_common_group + "/ntod", [len(tod_chunk_i)])
        comm_tod.add_field(pid_data_group + "/tod", tod_chunk_i)
        comm_tod.add_field(pid_data_group + "/pix", pix_chunk_i)
        comm_tod.add_field(pid_data_group + "/psi", psi_chunk_i)
    comm_tod.finalize_chunk(pid+1)


comm_tod.finalize_file()
