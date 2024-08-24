module compsep
  use iso_c_binding
  implicit none

  real(c_double) :: pi=3.141592653589793238462643383279502884197_c_double;
  
  type dataset
     integer(c_int64_t)        :: nside, npix, lmax
     real(c_double), allocatable, dimension(:) :: map
     real(c_double), allocatable, dimension(:) :: rms
     real(c_double), allocatable, dimension(:) :: b_l
  end type dataset
  
  integer(c_int64_t) :: numband
  type(dataset), allocatable, dimension(:) :: data
  
end module compsep

subroutine compsep_init(numband_arg) bind(c, name="compsep_init")
  use compsep
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in) :: numband_arg
  print *, "compsep_init called with ",numband_arg
  numband = numband_arg
  allocate(data(numband))
  
end subroutine compsep_init

subroutine compsep_init_band(i, nside, lmax, fwhm) bind(c, name="compsep_init_band")
  use compsep
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in) :: i, nside, lmax
  real(c_double),     intent(in) :: fwhm
  
  integer(4) :: l

  print *, "compsep_init_band called with ", i, nside, lmax, fwhm
  
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
  
end subroutine compsep_init_band

subroutine compsep_compute_rhs(rhs, l_rhs) bind(c, name="compsep_compute_rhs")
  use compsep
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in)  :: l_rhs
  real(c_double),     intent(out)        :: rhs(l_rhs)

  integer(4) :: i
  
  rhs = 0.d0
  do i = 1, numband
     rhs = data(i)%map / data(i)%rms**2
  end do
  
end subroutine compsep_compute_rhs

subroutine compsep_compute_ax(x, l_x) bind(c, name="compsep_compute_ax")
  use compsep
  use iso_c_binding
  implicit none
  integer(c_int64_t),  intent(in) :: l_x
  real(c_double),      intent(inout)       :: x(l_x)

  integer(4) :: i
  real(c_double), allocatable, dimension(:) :: y

  allocate(y(l_x))

  y = 0.d0
  do i = 1, numband
     y = y + x / data(i)%rms**2
  end do
  
end subroutine compsep_compute_Ax
