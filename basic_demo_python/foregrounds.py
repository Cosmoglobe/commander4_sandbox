import astropy.units as u
import astropy.constants as c
import numpy as np

A = (2*c.h*u.GHz**3/c.c**2).to('MJy').value
h_over_k = (c.h/c.k_B/(1*u.K)).to('GHz-1').value
h_over_kTCMB = (c.h/c.k_B/(2.7255*u.K)).to('GHz-1').value
def blackbody(nu, T):
    return nu**3/np.expm1(nu*h_over_k/T)

def g(nu):
    # From uK_CMB to MJy/sr
    x = nu*h_over_kTCMB
    return np.expm1(x)**2/(x**4*np.exp(x))

def cmb_sed(freq):
    # Assuming we are working in uK_CMB units
    return np.ones_like(freq)

def dust_sed(nu, beta, T, nu0 = 857):
    # Modified blackbody, in uK_CMB
    x = nu*h_over_k/T
    return g(nu)/g(nu0) * (nu/nu0)**beta * blackbody(nu, T)/blackbody(nu0, T)


def lnlike_beta_d(beta_d, args):
    '''
    data - element of Gibbs class
    pix_ind - pixel index to evaluate sky signal over
    '''
    data, pix_ind = args

    s_cmb  = cmb_sed(data.freqs)*data.comp_maps[0, pix_ind]
    s_dust = dust_sed(data.freqs, beta_d, data.T_d[pix_ind])*data.comp_maps[1, pix_ind]
    signal = s_cmb + s_dust
    lnlike = -0.5*(((signal - data.map_sky[:,pix_ind])/data.map_rms[:,pix_ind])**2).sum()

    return lnlike

def lnlike_T_d(T_d, args):
    '''
    data - element of Gibbs class
    pix_ind - pixel index to evaluate sky signal over
    '''
    data, pix_ind = args

    s_cmb  = cmb_sed(data.freqs)*data.comp_maps[0, pix_ind]
    s_dust = dust_sed(data.freqs, data.beta_d[pix_ind], T_d)*data.comp_maps[1, pix_ind]
    signal = s_cmb + s_dust
    lnlike = -0.5*(((signal - data.map_sky[:,pix_ind])/data.map_rms[:,pix_ind])**2).sum()

    return lnlike
