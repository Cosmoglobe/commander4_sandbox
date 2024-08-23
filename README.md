# commander4_sandbox

Setup:
- Compsep + TOD
- 2 bands at Nside=2048, lmax=6000, fwhm=10 arcmin
- data = CMB+white noise

- Compsep solve for CMB only, full-sky, brute-force in alm space,
  (S^-1 + sum_nu B^t*N^-1*B) s = sum_nu B^t*N^-1*map
- TOD: d = P*B*s + n, P = diag(Npix), ie., visit each pixel just once, results in uniform rms

- Both compsep and TOD run on two nodes each

- Define minimal but realistic parameter files. Absolute requirements:
  - easily editable with standard text editors, such as emacs and vim
  - human readable 
