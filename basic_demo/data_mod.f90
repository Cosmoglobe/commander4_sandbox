module data_mod
  use iso_c_binding
  implicit none

  real(c_double) :: pi=3.141592653589793238462643383279502884197_c_double;
  
  type :: scan_class
     integer(c_int64_t) :: ntod
     real(c_float)      :: sigma0
     real(c_float),      allocatable :: d(:)              ! Detector values in time domain
     integer(c_int32_t), allocatable :: pix(:)            ! pointer array of pixels 
  end type scan_class

  type :: tod_class
     integer(c_int64_t)       :: nscan
     type(scan_class), allocatable, dimension(:) :: scans    ! Array of all scans
  end type tod_class
  
  type dataset_class
     integer(c_int64_t)        :: nside, npix, lmax
     real(c_double), allocatable :: map(:)
     real(c_double), allocatable :: rms(:)
     real(c_double), allocatable :: b_l(:)
     type(tod_class) :: tod
  end type dataset_class
  
  integer(c_int64_t) :: numband
  type(dataset_class), allocatable, dimension(:) :: data
  
end module data_mod

subroutine data_init(numband_arg) bind(c, name="data_init")
  use data_mod
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in) :: numband_arg
  print *, "data_init called with ",numband_arg
  numband = numband_arg
  allocate(data(numband))
  
end subroutine data_init

subroutine data_init_band(i, nside, lmax, fwhm) bind(c, name="data_init_band")
  use data_mod
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in) :: i, nside, lmax
  real(c_double),     intent(in) :: fwhm
  
  integer(4) :: l

  print *, "data_init_band called with ", i, nside, lmax, fwhm
  
  data(i)%nside = nside
  data(i)%npix  = 12*nside**2
  data(i)%lmax  = lmax
  allocate(data(i)%map(0:data(i)%npix-1))
  allocate(data(i)%rms(0:data(i)%npix-1))
  allocate(data(i)%b_l(0:lmax))
  
  data(i)%map =   i
  data(i)%rms = 2*i
  do l = 0, lmax
     data(i)%b_l(l) = exp(-0.5*l*(l+1)*(fwhm*180*60*pi/sqrt(8.d0*log(2.d0))))
  end do
  
end subroutine data_init_band
