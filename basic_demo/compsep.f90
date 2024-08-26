module compsep_mod
use data_mod
use iso_c_binding
implicit none

contains

subroutine compsep_compute_rhs(rhs, l_rhs)
  integer(c_int64_t), intent(in)  :: l_rhs
  real(c_double),     intent(out) :: rhs(l_rhs)

  integer(4) :: i
  
  rhs = 0.d0
  do i = 1, numband
     rhs = data(i)%map / data(i)%rms**2
  end do
  
end subroutine compsep_compute_rhs

subroutine compsep_compute_Ax(x, l_x)
  integer(c_int64_t),  intent(in)    :: l_x
  real(c_double),      intent(inout) :: x(l_x)

  integer(4) :: i
  real(c_double), allocatable, dimension(:) :: y

  allocate(y(l_x))

  y = 0.d0
  do i = 1, numband
     y = y + x / data(i)%rms**2
  end do
  
end subroutine compsep_compute_Ax

end module compsep_mod


subroutine compsep_compute_rhs_ifc(rhs, l_rhs) bind(c, name="compsep_compute_rhs_ifc")
  use compsep_mod
  implicit none
  integer(c_int64_t), intent(in)  :: l_rhs
  real(c_double),     intent(out) :: rhs(l_rhs)

  call compsep_compute_rhs(rhs, l_rhs)  
end subroutine compsep_compute_rhs_ifc

subroutine compsep_compute_Ax_ifc(x, l_x) bind(c, name="compsep_compute_Ax_ifc")
  use compsep_mod
  implicit none
  integer(c_int64_t),  intent(in)    :: l_x
  real(c_double),      intent(inout) :: x(l_x)

  call compsep_compute_Ax(x, l_x)
end subroutine compsep_compute_Ax_ifc
