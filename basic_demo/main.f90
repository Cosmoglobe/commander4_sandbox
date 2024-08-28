program main
  use iso_c_binding
  use data_mod
  use tod_mod
  use compsep_mod
  implicit none

  integer(c_int64_t) :: nband, nside, npix, lmax, i, j, k, nscan, ntod, ngibbs, iter
  real(c_double)     :: fwhm

  real(c_double),     allocatable  :: rhs(:), signal(:), m_map(:), m_rms(:)
  real(c_float),      allocatable  :: d(:)
  integer(c_int32_t), allocatable  :: pix(:)

  ngibbs = 5
  nband  = 5
  nscan  = 108
  ntod   = 2**16
  nside  = 256
  npix   = 12*nside**2
  lmax   = 512
  fwhm   = 0.42

  write(*,*) nband
  call data_init(nband)
  
  ! Initialize basic data  
  do i = 1, nband
     call data_init_band(i, nside, lmax, fwhm)
  end do

  ! Initialize TOD data
  do i = 1, nband
     call tod_init_band(i, nscan)
     do j = 1, nscan
        allocate(d(ntod), pix(ntod))
        do k = 1, ntod
           d(k)   = mod(k,4) + 0.6d0
           pix(k) = mod(k,npix)
        end do
        call tod_init_scan(i, j, d, pix)
        deallocate(d,pix)
     end do
  end do  


  ! Run Gibbs sampler
  do iter = 1, ngibbs

     write(*,*) 'Gibbs iter = ', iter, ' of ', ngibbs
     
     ! **********************
     ! COMPSEP stage
     ! **********************
     
     ! Compute RHS of mapmaking equation
     allocate(rhs(0:12*nside**2-1))
     call compsep_compute_rhs(rhs)
     
     ! Solve for best-fit map by CG
     allocate(signal(0:12*nside**2-1))
     call compsep_compute_Ax(signal)
     
     ! **********************
     ! TOD stage
     ! **********************

     ! Estimate white noise rms per scan
     do i = 1, nband
        do j = 1, size(data(i)%tod%scans,1)
           call tod_estimate_sigma0(i, j, signal)
        end do
     end do

     ! Make frequency maps
     do i = 1, nband
        call tod_mapmaker(i)
     end do

     ! Clean up
     deallocate(rhs, signal)
  end do

  
end program main
