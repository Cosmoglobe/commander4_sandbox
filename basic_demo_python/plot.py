import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import os
import tqdm
import camb
#import cosmoglobe as cg

try:
    os.mkdir('plots')
except FileExistsError:
    pass

NSIDE   = 64
NPIX = 12*NSIDE**2
LMAX = 3*NSIDE-1


iband = 0

map_sky = hp.read_map(f"output/map_band_{iband:02}_c000001.fits")
hp.mollview(map_sky, title=f"Observed sky", min=-250, max=250, cmap="RdBu_r")
plt.savefig(f"plots/map_obs.png")
map_rms = hp.read_map(f"output/rms_band_{iband:02}_c000001.fits")
hp.mollview(map_rms, title=f"rms")
plt.savefig(f"plots/rms.png")

pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06,
                       As=2e-9, ns=0.965, halofit_version='mead', lmax=LMAX)

results = camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
totCL=powers['total']

ell = np.arange(LMAX+1)
Dl = totCL[ell,0]

for iter in tqdm.tqdm(range(1, 1000)):
    if os.path.exists(f"plots/Cl_iter{iter:06}.png"):
        continue


    # Component maps
    m = hp.read_map(f'output/comp_dust_c{iter:06}.fits')
    hp.mollview(m, norm='log', min=75e3, max=1e9, cmap='afmhot')
    plt.savefig(f'plots/comp_dust_c{iter:06}.png')
    m = hp.read_map(f'output/comp_cmb_c{iter:06}.fits')
    hp.mollview(m, min=-250, max=250, cmap='RdBu_r')
    plt.savefig(f'plots/comp_cmb_c{iter:06}.png')

    
    m = hp.read_map(f'output/comp_dust_beta_c{iter:06}.fits')
    hp.mollview(m, min=1.4, max=1.6)
    plt.savefig(f'plots/comp_dust_beta_c{iter:06}.png')
    m = hp.read_map(f'output/comp_dust_T_c{iter:06}.fits')
    hp.mollview(m, min=17, max=23)
    plt.savefig(f'plots/comp_dust_T_c{iter:06}.png')


    # Residual maps
    for iband in range(5):
        res = hp.read_map(f'output/res_band_{iband:02}_c{iter:06}.fits')
        rms = hp.read_map(f'output/rms_band_{iband:02}_c{iter:06}.fits')

        hp.mollview(res/rms, min=-3, max=3, cmap='RdBu_r')
        plt.savefig(f'plots/res_sigma_band_{iband:02}_c{iter:06}.png')
        hp.mollview(rms)
        plt.savefig(f'plots/rms_band_{iband:02}_c{iter:06}.png')

        m = hp.read_map(f'output/map_band_{iband:02}_c{iter:06}.fits')
        hp.mollview(m, min=-250, max=250, cmap='RdBu_r')
        plt.savefig(f'plots/map_band_{iband:02}_c{iter:06}.png')
        plt.close('all')
   

    # CMB Constrained Realizations
    mapx = hp.read_map(f"output/mapx_c{iter:06}.fits") 
    hp.mollview(mapx, title=f"mean field map, iter {iter}", min=-250, max=250, cmap="RdBu_r")
    plt.savefig(f"plots/mapx_iter{iter:06}.png")
    mapy = hp.read_map(f"output/mapy_c{iter:06}.fits")
    hp.mollview(mapy, title=f"fluctuation map, iter {iter}", min=-250, max=250, cmap="RdBu_r")
    plt.savefig(f"plots/mapy_iter{iter:06}.png")
    hp.mollview(mapx+mapy, title=f"Full sky realization, iter {iter}", min=-250, max=250, cmap="RdBu_r")
    plt.savefig(f"plots/mapx+y_iter{iter:06}.png")

    Cl_sample = np.load(f"output/Cl_sample_c{iter:06}.npy")

    Cl_sky = hp.alm2cl(hp.map2alm(map_sky))
    Clx = hp.alm2cl(hp.map2alm(mapx))
    Clxy = hp.alm2cl(hp.map2alm(mapx+mapy))
    plt.figure()
    Z = ell*(ell+1)/(2*np.pi)
    plt.semilogy(ell[2:], Z[2:]*Clx[2:], label="Mean field map Cl")
    plt.semilogy(ell[2:], Z[2:]*Clxy[2:], label="Mean field + fluct Cl")
    plt.semilogy(ell[2:], Z[2:]*Cl_sample[2:], label="Gibbs sampled Cl")
    plt.semilogy(ell[2:], Dl[2:], label="True Cl", color='k', zorder=10)
    plt.semilogy(ell[2:], Z[2:]*Cl_sky[2:], label="Observed sky Cl")
    plt.legend()
    plt.ylim(100, 6500)
    plt.title(f"Dl, iter {iter}")
    plt.savefig(f"plots/Cl_iter{iter:06}.png")
    plt.close('all')

