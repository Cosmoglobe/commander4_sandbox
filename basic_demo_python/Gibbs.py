import numpy as np
import healpy as hp
import h5py
import ducc0
from pixell import utils
import time
import ctypes as ct
from scipy.sparse.linalg import cg, LinearOperator
import os

# number of threads that ducc0 should use. NOte that this can be varied on a
# call-by-call basis if necessary.
nthreads = 20

# Some global variables - should be read from a parameter file.
FWHM    = 0.16666*np.pi/180.0*np.ones(5)
NSIDE   = 256
LMAX = 3*NSIDE-1
NTOD = 2**16
NSCAN = 108
VERBOSE = False

ell = np.arange(LMAX+1)
Cl_true = np.zeros(len(ell))
Cl_true[2:] = 1./ell[2:]**2  # True Cl used in the "make_simple_sim.py" script.
Cl_true[:2] = 1e-12


def dot_alm(alm1, alm2):
    """ Function calculating the dot product of two alms, given that they follow the Healpy standard,
        where alms are represented as complex numbers, but with the conjugate 'negative' ms missing.
    """
    return np.sum((alm1[:LMAX]*alm2[:LMAX]).real) + np.sum((alm1[LMAX:]*np.conj(alm2[LMAX:])).real*2)


def alm2map(alm, nside, lmax):
    base = ducc0.healpix.Healpix_Base(nside, "RING")
    geom = base.sht_info()
    return ducc0.sht.synthesis(alm=alm.reshape((1,-1)),
                               lmax=lmax,
                               spin=0,
                               nthreads=nthreads, **geom).reshape((-1,))


def alm2map_adjoint(map, nside, lmax):
    base = ducc0.healpix.Healpix_Base(nside, "RING")
    geom = base.sht_info()
    return ducc0.sht.adjoint_synthesis(map=map.reshape((1,-1)),
                                       lmax=lmax,
                                       spin=0,
                                       nthreads=nthreads, **geom).reshape((-1,))


class Gibbs:
    def __init__(self):
        self.ntod   = NTOD
        self.nside  = NSIDE
        self.npix = 12*self.nside**2
        self.lmax = LMAX
        self.alm_len = ((LMAX+1)*(LMAX+2))//2
        self.nscan  = NSCAN
        self.nband = 5
        self.fwhm = FWHM[:self.nband]

        self.map_rms = np.zeros((self.nband, self.npix))
        self.map_inv_var = np.zeros((self.nband, self.npix))
        self.map_sky = np.zeros((self.nband, self.npix))
        self.tod_signal_sample = np.zeros((self.nband, self.nscan, self.ntod))
        self.Cl_sample = np.zeros(self.lmax) + np.inf

    def read_tod_from_file(self, h5_filename, bands):
        """ Reads TOD from a h5 file generated by 'make_simple_sims.py'. Currently reads only one band for debuging purposes.
        """
        self.tod = np.zeros((self.nband, self.nscan, self.ntod), dtype=np.float64)
        self.pix = np.zeros((self.nband, self.nscan, self.ntod), dtype=int)
        with h5py.File(h5_filename) as f:
            for iband in range(self.nband):
                for iscan in range(self.nscan):
                    self.tod[iband, iscan] = f[f'{iscan+1:06}/{bands[iband]}/tod'][:self.ntod]
                    self.pix[iband, iscan] = f[f'{iscan+1:06}/{bands[iband]}/pix'][:self.ntod]


    def LHS_func(self, x):
        """ The LHS of equations 5 and 6 from Eriksen 2004, implemented as a function on the alm-vector x.
            The equation can be written as (C^-1 x + A^T Y^T N Y A x), where Y is a map->alm conversion, and Y^T is map-> alm,
            A is the beam-smoothing, and C is the current C(ell) sample.
        """
        LHS_sum = np.zeros_like(x)
        LHS_sum += hp.almxfl(x, 1.0/self.Cl_sample)
        for iband in range(self.nband):
            Ax = hp.smoothalm(x, self.fwhm[iband], inplace=False)
            YAx = alm2map(Ax, self.nside, self.lmax)
            NYAx = YAx.copy()/self.map_rms[iband]**2
            YTNYAx = alm2map_adjoint(NYAx, self.nside, self.lmax)
            ATYTNYAx = hp.smoothalm(YTNYAx, self.fwhm[iband], inplace=False)
            LHS_sum += ATYTNYAx
        return LHS_sum


    def get_RHS_eqn_mean(self):
        """ Calculates and returns the RHS of the mean-field (Wiener filtered) map equation (eqn 5 from Eriksen 2004).
            This RHS can be written as (A^T Y^T N d), where d is the observed sky, and Y^T is a map->alm conversion,
            N is the noise covariance, and A is the beam.
        """
        RHS_sum = np.zeros(self.alm_len, dtype=np.complex128)
        for iband in range(self.nband):
            Nd = self.map_sky[iband]/self.map_rms[iband]**2
            YTNd = alm2map_adjoint(Nd, self.nside, self.lmax)
            ATYTNd = hp.smoothalm(YTNd, self.fwhm[iband], inplace=False)
            RHS_sum += ATYTNd
        return RHS_sum


    def get_RHS_eqn_fluct(self):
        """ Calculates and returns the RHS of the map fluctuation equation (eqn 6 from Eriksen 2004).
            This RHS can be written as (C^-1/2 Y^T omega0 + A^T Y^T N^-1/2 omega1), where omega0 and omega1 are N(0,1) maps,
            and C is the currentl C(ell) sample.
        """
        # YTomega0 = hp.map2alm(np.random.normal(0, 1, self.npix), iter=0)#*(4*np.pi/self.npix)
        # CYTomega0 = hp.almxfl(YTomega0, np.sqrt(1.0/self.Cl_sample))
        RHS_sum = np.zeros(self.alm_len, dtype=np.complex128)
        CYTomega0 = hp.synalm(1.0/self.Cl_sample, self.lmax)
        RHS_sum += CYTomega0

        for iband in range(self.nband):
            omega1 = np.random.normal(0, 1, self.npix)
            Nomega1 = omega1/self.map_rms[iband]
            YTNomega1 = alm2map_adjoint(Nomega1, self.nside, self.lmax)
            ATYTNomega1 = hp.smoothalm(YTNomega1, self.fwhm[iband], inplace=False)
            RHS_sum += ATYTNomega1
        return RHS_sum


    def compsep_compute_Ax(self, LHS, RHS):
        """ Solves the equation Ax=b for x given A (LHS) and b (RHS) using CG from the pixell package.
            Assumes that both x and b are in alm space.

            Args:
                LHS: A callable taking x as argument and returning Ax.
                RHS: A Numpy array representing b, in alm space.
            Returns:
                m_bestfit: The resulting best-fit solution, in alm space.
        """
        CG_solver = utils.CG(LHS, RHS, dot=dot_alm)
        err_tol = 1e-6
        maxiter = 251
        iter = 0
        while CG_solver.err > err_tol:
            CG_solver.step()
            iter += 1
            if VERBOSE and iter%10 == 1:
                print(f"CG iter {iter:3d} - Residual {CG_solver.err:.3e}")
            if iter >= maxiter:
                print(f"Warning: Maximum number of iterations ({maxiter}) reached in CG.")
                break
        print(f"CG finished after {iter} iterations with a residual of {CG_solver.err:.3e} (err tol = {err_tol})")
        s_bestfit = CG_solver.x

        return s_bestfit


    def tod_estimate_sigma0(self):
        """ Noise estimation (currently only white noise level per scan).
            Calculates sigma0 = std(d_i+1 - d_i)/sqrt(2).
        """
        self.sigma0_est = np.std(self.tod_signalsubtracted[...,1:] - self.tod_signalsubtracted[...,:-1], axis=-1)/np.sqrt(2)


    def tod_mapmaker_purepython(self):
        """ A simple binned mapmaker for making observed sky and rms maps from the TOD scans and rms estimates.
            The only thing that changes in this function between Gibbs iterations is the sigma0-estimate.
        """
        self.map_sky[:] = 0.0
        self.map_inv_var[:] = 0.0
        self.map_rms[:] = 0.0
        for iband in range(self.nband):
            for iscan in range(self.nscan):
                self.map_sky[iband] += np.bincount(self.pix[iband,iscan], weights=self.tod[iband,iscan]/self.sigma0_est[iband,iscan]**2, minlength=self.npix)
                self.map_inv_var[iband] += np.bincount(self.pix[iband,iscan], weights=1.0/self.sigma0_est[iband,iscan]**2*np.ones(self.ntod), minlength=self.npix)
        self.map_rms = 1.0/np.sqrt(self.map_inv_var)
        self.map_sky /= self.map_inv_var
        ipix_mask = hp.query_disc(self.nside, (10,0,0), np.radians(60))  # Quick way of simulating a "mask", aka a region of infinite rms.
        self.map_rms[:,ipix_mask] = np.inf


    def tod_mapmaker(self):
        self.map_sky[:] = 0.0
        self.map_inv_var[:] = 0.0
        self.map_rms[:] = 0.0
        maplib = ct.cdll.LoadLibrary("/home/jonas/github/commander4_sandbox/basic_demo_python/mapmaker.so")
        ct_i64_dim2 = np.ctypeslib.ndpointer(dtype=ct.c_int64, ndim=2, flags="contiguous")
        ct_f64_dim1 = np.ctypeslib.ndpointer(dtype=ct.c_double, ndim=1, flags="contiguous")
        ct_f64_dim2 = np.ctypeslib.ndpointer(dtype=ct.c_double, ndim=2, flags="contiguous")
        # Replace maplib.mapmaker with maplib.mapmakerOMP for parallelized version (worth it if pixels are hit many times, >>10).
        maplib.mapmaker.argtypes = [ct_f64_dim1, ct_f64_dim1, ct_f64_dim2, ct_i64_dim2, ct_f64_dim1, ct.c_int64, ct.c_int64, ct.c_int64]
        for iband in range(self.nband):
            maplib.mapmaker(self.map_sky[iband], self.map_inv_var[iband], self.tod[iband], self.pix[iband], self.sigma0_est[iband], self.ntod, self.nscan, self.npix)
        self.map_rms = 1.0/np.sqrt(self.map_inv_var)
        ipix_mask = hp.query_disc(self.nside, (10,0,0), np.radians(60))  # Quick way of simulating a "mask", aka a region of infinite rms.
        self.map_rms[:,ipix_mask] = np.inf


    def map2tod(self):
        """ From the current realization of CMB mean field + fluctuation maps, create a TOD by reprojecting onto the pointing.
        """
        for iband in range(self.nband):
            for iscan in range(self.nscan):
                self.tod_signal_sample[iband,iscan] = self.map_signal_mean[self.pix[iband,iscan]] + self.map_signal_fluct[self.pix[iband,iscan]]


    def sample_Cl(self):
        """ Sample C(ell) from the current CMB estimate.
        """
        Cl_signal = hp.alm2cl(self.alm_signal_mean + self.alm_signal_fluct)
        rho = np.zeros(self.lmax+1)
        for i in range(self.lmax):
            rho[i+1] = np.sum(np.random.normal(0, 1, 2*(i+1))**2)/(2*(i+1))
        rho[0] = np.inf
        self.Cl_sample = Cl_signal/rho
        self.Cl_sample[:2] = 1e-6
        # Hard-coding the Cl-sampling to the true Cl can be a useful debug tool:
        # self.Cl_sample = Cl_true


    def solve(self, niter):
        for iter in range(1, niter+1):
            print(f'#### Gibbs iter = {iter} of {ngibbs} ####')

            # **********************
            # TOD stage
            # **********************
            # Estimate white noise rms per scan
            t0 = time.time()
            self.tod_signalsubtracted = self.tod - self.tod_signal_sample
            self.tod_estimate_sigma0()
            print(f">TOD sampling finished in {time.time()-t0:.2f}s.")
            for iband in range(self.nband):
                np.save(f"output/sigma0_est_band_{iband:02}_c{iter:06}.npy", self.sigma0_est)

            # **********************
            # Mapmaking stage
            # **********************
            # Make frequency maps
            t0 = time.time()
            self.tod_mapmaker()
            print(f">Mapmaker finished in {time.time()-t0:.2f}s.")
            # Write maps to file.
            for iband in range(self.nband):
                hp.write_map(f'output/map_band_{iband:02}_c{iter:06}.fits', self.map_sky,
                        overwrite=True, dtype=np.float64)
                hp.write_map(f'output/rms_band_{iband:02}_c{iter:06}.fits', self.map_rms,
                        overwrite=True, dtype=np.float64)

            # **********************
            # COMPSEP stage
            # **********************
            t0 = time.time()
            # Compute RHS of comp-sep equation
            compsep_RHS_eqn_mean = self.get_RHS_eqn_mean()
            compsep_RHS_eqn_fluct = self.get_RHS_eqn_fluct()
            # Solve for best-fit map by CG
            compsep_LHS = LinearOperator(shape=((self.alm_len, self.alm_len)), matvec=self.LHS_func, dtype=np.complex128)
            self.alm_signal_mean = self.compsep_compute_Ax(compsep_LHS, compsep_RHS_eqn_mean)
            self.map_signal_mean = alm2map(self.alm_signal_mean, self.nside, self.lmax)
            self.alm_signal_fluct = self.compsep_compute_Ax(compsep_LHS, compsep_RHS_eqn_fluct)
            self.map_signal_fluct = alm2map(self.alm_signal_fluct, self.nside, self.lmax)
            print(f">CompSep finished in {time.time()-t0:.2f}s.")
            # Re-project the sky realization onto the TOD
            t0 = time.time()
            self.map2tod()
            print(f">Map->TOD projection finished in {time.time()-t0:.2f}s.")

            # **********************
            # C(ell) sampling
            # **********************
            t0 = time.time()
            self.sample_Cl()
            print(f">Cl sampling finished in {time.time()-t0:.2f}s.")

            # Write Wiener filtered and fluctuation maps to file.
            hp.write_map(f'output/mapx_c{iter:06}.fits', self.map_signal_mean,
                    overwrite=True, dtype=np.float64)
            hp.write_map(f'output/mapy_c{iter:06}.fits', self.map_signal_fluct,
                    overwrite=True, dtype=np.float64)
            np.save(f"output/Cl_sample_c{iter:06}.npy", self.Cl_sample)


if __name__ == "__main__":
    try:
        os.mkdir('output')
    except FileExistsError:
        pass
    np.random.seed(128)
    ngibbs = 250
    gibbs = Gibbs()
    gibbs.read_tod_from_file('../src/python/preproc_scripts/tod_example.h5', ['030', '070', '100', '217', '353'])
    gibbs.solve(ngibbs)
