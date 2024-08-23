module compsep
  use iso_c_binding
  implicit none
  
  type dataset
     character(len=512)  :: label
     integer(c_int64_t)        :: nside, npix, lmax
     real(c_double), allocatable, dimension(:) :: map
     real(c_double), allocatable, dimension(:) :: rms
     real(c_double), allocatable, dimension(:) :: b_l
  end type dataset
  
  integer(c_int64_t) :: numband
  type(dataset), allocatable, dimension(:) :: data
  
contains

  subroutine init(numband) bind(c, name="init")
    implicit none
    integer(c_int64_t), intent(in) :: numband

    allocate(data(numband))

  end subroutine init

  subroutine init_band(i, label, nside, lmax, fwhm) bind(c, name="init_band")
    implicit none
    character(len=*),   intent(in) :: label
    integer(c_int64_t), intent(in) :: numband

    integer(4) :: l
    
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
    
  end subroutine init_band

  subroutine compute_rhs(i, rhs) bind(c, name="init_band")
    implicit none
    integer(c_int64_t),               intent(in)  :: i
    real(c_double),     dimension(:), intent(out) :: rhs

    rhs = data(i)%map / data(i)%rms**2
    
  end subroutine compute_rhs

  subroutine compute_Ax(i, x) bind(c, name="init_band")
    implicit none
    integer(c_int64_t),               intent(in)    :: i
    real(c_double),     dimension(:), intent(inout) :: x

    x = x / data(i)%rms**2
    
  end subroutine compute_Ax
  
end module compsep
