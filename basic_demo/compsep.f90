module compsep_mod
use data_mod
use iso_c_binding
implicit none

contains

subroutine compsep_compute_rhs(rhs)
  real(c_double), dimension(0:), intent(out) :: rhs

  integer :: i

  rhs = 0.d0
  do i = 1, size(data,1)
     rhs = rhs + data(i)%map / data(i)%rms**2
  end do
end subroutine compsep_compute_rhs

subroutine compsep_compute_Ax(x)
  real(c_double), dimension(:), intent(inout) :: x

  integer :: i
  real(c_double), allocatable, dimension(:) :: y

  allocate(y(size(x,1)))

  y = 0.d0
  do i = 1, size(data,1)
     y = y + x / data(i)%rms**2
  end do
end subroutine compsep_compute_Ax

end module compsep_mod

! foreign language interface routines ... simple forwarders to module members

subroutine compsep_compute_rhs_ifc(rhs, l_rhs) bind(c, name="compsep_compute_rhs_ifc")
  use compsep_mod
  implicit none
  integer(c_int64_t), intent(in), value  :: l_rhs
  real(c_double),     intent(out) :: rhs(l_rhs)

  call compsep_compute_rhs(rhs)  
end subroutine compsep_compute_rhs_ifc

subroutine compsep_compute_Ax_ifc(x, l_x) bind(c, name="compsep_compute_Ax_ifc")
  use compsep_mod
  implicit none
  integer(c_int64_t),  intent(in), value    :: l_x
  real(c_double),      intent(inout) :: x(l_x)

  call compsep_compute_Ax(x)
end subroutine compsep_compute_Ax_ifc
