module compsep
  use iso_c_binding
  implicit none

  real(c_double) :: pi=3.141592653589793238462643383279502884197_c_double;
  
  type dataset
     character(len=512)  :: label
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
  integer(c_int64_t), intent(in), value :: numband_arg
  print *, "compsep_init called with ",numband_arg
  numband = numband_arg
  allocate(data(numband))
  
end subroutine compsep_init

subroutine compsep_init_band(i, label_raw, label_len, nside, lmax, fwhm) bind(c, name="compsep_init_band")
  use compsep
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in), value :: i, nside, lmax, label_len
  real(c_double), intent(in), value :: fwhm
  integer(c_int8_t), dimension(label_len), intent(in) :: label_raw
  character(len=label_len) label
  
  integer(4) :: l
  do l=1,label_len
    label(l:l) = achar(label_raw(l))
  end do

  print *, "compsep_init_band called with ",i, label, nside, lmax, fwhm
  
  data(i)%label = label
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

subroutine compsep_compute_rhs(i, rhs, l_rhs) bind(c, name="compsep_compute_rhs")
  use compsep
  use iso_c_binding
  implicit none
  integer(c_int64_t),               intent(in), value  :: i, l_rhs
  real(c_double), intent(out) :: rhs(l_rhs)
  
  rhs = data(i)%map / data(i)%rms**2
  
end subroutine compsep_compute_rhs

subroutine compsep_compute_Ax(i, x, l_x) bind(c, name="compsep_compute_Ax")
  use compsep
  use iso_c_binding
  implicit none
  integer(c_int64_t),               intent(in), value    :: i, l_x
  real(c_double), intent(inout) :: x(l_x)
  
  x = x / data(i)%rms**2
  
end subroutine compsep_compute_Ax
