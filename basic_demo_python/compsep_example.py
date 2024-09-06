import numpy as np
import healpy as hp
import astropy.units as u
from tqdm import tqdm
import matplotlib.pyplot as plt

def g(nu):
    # From uK_CMB to MJy/sr
    x = nu/56.78
    return np.expm1(x)**2/(x**4*np.exp(x))

from astropy.modeling.physical_models import BlackBody

T = 20
beta = 1.5
bb = BlackBody(temperature=T*u.K)

ms = []
rmss = []
for i in range(5):
    ms.append(hp.read_map(f'output/map_band_{i:02}_c000001.fits'))
    rmss.append(hp.read_map(f'output/rms_band_{i:02}_c000001.fits'))

ms = np.array(ms)
rmss = np.array(rmss)

freqs = np.array([30,70,100,217,353])

T = np.zeros((5, 2))
T[:,0] = 1
T[:,1] = g(freqs)/g(857)*bb(freqs*u.GHz)/bb(857*u.GHz)*(freqs/857)**beta

m_cmb = hp.read_map('../src/python/preproc_scripts/true_sky_cmb_256.fits')
m_dust = hp.read_map('../src/python/preproc_scripts/true_sky_dust857_256.fits')

npix = 12*256**2
b = np.zeros((2,npix))
den = np.zeros((2,2,npix))
x = np.zeros((2,npix))
for i in tqdm(range(npix)):
    mean = T.T.dot((1/np.ones_like(rmss[:,i]**2)*ms[:,i]))
    fluc = T.T.dot(np.random.randn(5)/np.ones_like(rmss[:,i]))
    den = (T.T.dot(np.diag(1/np.ones_like(rmss[:,i]**2))).dot(T))
    try:
        x[:,i] = np.linalg.solve(den, mean + fluc)
    except np.linalg.LinAlgError:
        x[:,i] = 0

hp.mollview(x[0], min=-250, max=250)
hp.mollview(m_cmb, min=-250, max=250)
hp.mollview(x[1], norm='log', min=70e3, max=1e9)
hp.mollview(m_dust, norm='log', min=70e3, max=1e9)
plt.show()
