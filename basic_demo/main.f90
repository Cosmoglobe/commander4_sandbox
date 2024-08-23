program main
  use iso_c_binding
  use compsep
  implicit none

  integer(c_int64_t) :: nband, nside, lmax, i
  real(c_double)     :: fwhm

  real(c_double), allocatable  :: rhs(:), x(:)
  
  nband = 5
  nside = 256
  lmax  = 512
  fwhm  = 0.42

  write(*,*) nband
  call compsep_init(nband)
  
  ! Initialize data  
  do i = 1, nband
     call compsep_init_band(i, nside, lmax, fwhm)
  end do
  
  ! Compute RHS of mapmaking equation
  allocate(rhs(0:12*nside**2-1))
  call compsep_compute_rhs(rhs, size(rhs))

  ! Solve for best-fit map by CG
  allocate(x(0:12*nside**2-1))
  call compsep_compute_Ax(x, size(x))

end program main
