module f90mod
implicit none
private

real(kind=8), allocatable, dimension(:,:), public :: mod_arr

end module f90mod


subroutine demo_sub(arr, nx, ny) bind(c, name="demo_sub")
  use iso_c_binding
  implicit none

  integer(c_int64_t), intent(in), value :: nx, ny
  real(c_double), intent(inout) :: arr(nx,ny)

  print *, "demo_sub called!"
  print *, "nx=", nx, "ny=", ny
  arr(:,:) = 5.
end subroutine

subroutine set_mod_arr(arr, nx, ny) bind(c, name="set_mod_arr")
  use iso_c_binding
  use f90mod
  implicit none

  integer(c_int64_t), intent(in), value :: nx, ny
  real(c_double), intent(in) :: arr(nx,ny)

  print *, "set_mod_arr called!"
  print *, "nx=", nx, "ny=", ny
  allocate(mod_arr(nx,ny))
  mod_arr = arr
end subroutine

subroutine get_mod_arr(arr, nx, ny) bind(c, name="get_mod_arr")
  use iso_c_binding
  use f90mod
  implicit none

  integer(c_int64_t), intent(in), value :: nx, ny
  real(c_double), intent(out) :: arr(nx,ny)

  print *, "get_mod_arr called!"
  print *, "nx=", nx, "ny=", ny
  arr = mod_arr
end subroutine
