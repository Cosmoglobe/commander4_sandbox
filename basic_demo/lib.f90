subroutine demo_sub(arr, nx, ny) bind(c, name="demo_sub")
  use iso_c_binding
  implicit none

  integer(c_int64_t), intent(in), value :: nx, ny
  real(c_double), intent(inout) :: arr(nx,ny)

  print *, "demo_sub called!"
  print *, "nx=", nx, "ny=", ny
  arr(:,:) = 5.
end subroutine
