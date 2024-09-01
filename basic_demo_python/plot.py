import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
NSIDE   = 256
NPIX = 12*NSIDE**2
LMAX = 3*NSIDE-1

maptrue = hp.read_map(f"../src/python/preproc_scripts/true_sky.fits")
hp.mollview(maptrue, title=f"True sky map", min=-10, max=10, cmap="RdBu_r")
plt.savefig(f"plots/maptrue.png")

map_sky = hp.read_map(f"output/map_band_00_c000001.fits")
hp.mollview(map_sky, title=f"Observed sky", min=-10, max=10, cmap="RdBu_r")
plt.savefig(f"plots/map_obs.png")
map_rms = hp.read_map("output/rms_band_00_c000001.fits")
hp.mollview(map_rms, title=f"rms")
plt.savefig(f"plots/rms.png")
map_hits = hp.read_map("output/hits_band_00_c000001.fits")
hp.mollview(map_hits, title=f"hitmap")
plt.savefig(f"plots/hits.png")

for iter in range(1, 5):
    mapx = hp.read_map(f"output/mapx_band_00_c{iter:06}.fits") 
    hp.mollview(mapx, title=f"mean field map, iter {iter}", min=-10, max=10, cmap="RdBu_r")
    plt.savefig(f"plots/mapx_iter{iter:06}.png")
    mapy = hp.read_map(f"output/mapy_band_00_c{iter:06}.fits")
    hp.mollview(mapy, title=f"fluctuation map, iter {iter}", min=-10, max=10, cmap="RdBu_r")
    plt.savefig(f"plots/mapy_iter{iter:06}.png")
    hp.mollview(mapx+mapy, title=f"Full sky realization, iter {iter}", min=-10, max=10, cmap="RdBu_r")
    plt.savefig(f"plots/mapx+y_iter{iter:06}.png")

    Cl_sample = np.load(f"output/Cl_sample_band_00_c{iter:06}.npy")
    ell = np.arange(LMAX+1)
    Cl = np.zeros(len(ell))
    Cl[2:] = 1./ell[2:]**2

    Cl_sky = hp.alm2cl(hp.map2alm(map_sky))
    Clx = hp.alm2cl(hp.map2alm(mapx))
    Clxy = hp.alm2cl(hp.map2alm(mapx+mapy))
    plt.figure()
    plt.loglog(ell[2:], Clx[2:], label="Mean field map Cl")
    plt.loglog(ell[2:], Clxy[2:], label="Mean field + fluct Cl")
    plt.loglog(ell[2:], Cl_sample[2:], label="Gibbs sampled Cl")
    plt.loglog(ell[2:], Cl[2:], label="True Cl")
    plt.loglog(ell[2:], Cl_sky[2:], label="Observed sky Cl")
    plt.legend()
    plt.ylim(1e-8, 1)
    plt.title(f"Cl, iter {iter}")
    plt.savefig(f"plots/Cl_iter{iter}.png")

