subroutine compsep_compute_rhs(rhs, l_rhs) bind(c, name="compsep_compute_rhs")
  use data_mod
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in)  :: l_rhs
  real(c_double),     intent(out) :: rhs(l_rhs)

  integer(4) :: i
  
  rhs = 0.d0
  do i = 1, numband
     rhs = data(i)%map / data(i)%rms**2
  end do
  
end subroutine compsep_compute_rhs

subroutine compsep_compute_Ax(x, l_x) bind(c, name="compsep_compute_ax")
  use data_mod
  use iso_c_binding
  implicit none
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
