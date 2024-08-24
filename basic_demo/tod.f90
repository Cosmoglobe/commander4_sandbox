subroutine tod_init_band(i, nscan) bind(c, name="tod_init_band")
  use data_mod
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in), value :: i, nscan
  
  integer(4) :: l

  print *, "tod_init_band called with ", i, nscan

  data(i)%tod%nscan = nscan
  allocate(data(i)%tod%scans(nscan))
  
end subroutine tod_init_band

subroutine tod_init_scan(band, scan, ntod, d, pix) bind(c, name="tod_init_scan")
  use data_mod
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in), value :: band, scan, ntod
  real(c_float),      intent(in), value :: d
  real(c_int32_t),    intent(in), value :: pix
  
  print *, "tod_init_scan called with ", band, scan, ntod, d, pix

  data(band)%tod%scans(scan)%ntod = ntod
  allocate(data(band)%tod%scans(scan)%d(ntod))
  allocate(data(band)%tod%scans(scan)%pix(ntod))
  
end subroutine tod_init_scan

subroutine tod_estimate_sigma0(band, scan, s) bind(c, name="tod_estimate_sigma0")
  use data_mod
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in), value  :: band, scan
  real(c_float),      intent(in), value  :: s

  integer(4) :: i
  real(c_double) :: sigmasq
  
  sigmasq = 0.d0
  do i = 1, data(band)%tod%scans(scan)%ntod
     sigmasq = sigmasq + (data(band)%tod%scans(scan)%d(i) -  s)**2
  end do
  data(band)%tod%scans(scan)%sigma0 = sqrt(sigmasq / (data(band)%tod%scans(scan)%ntod-1))
  
end subroutine tod_estimate_sigma0

subroutine tod_mapmaker(band) bind(c, name="tod_mapmaker")
  use data_mod
  use iso_c_binding
  implicit none
  integer(c_int64_t),  intent(in),   value :: band

  integer(4) :: i, pix
  real(c_double), allocatable, dimension(:) :: A, b

  allocate(A(0:data(band)%npix-1), b(0:data(band)%npix-1))
  A = 0.d0
  b = 0.d0
  do i = 1, data(band)%tod%nscan
     pix    = data(band)%tod%scans(i)%pix
     A(pix) = A(pix) + 1.d0
     b(pix) = b(pix) + data(band)%tod%scans(i)%d(pix)
  end do
  
end subroutine tod_mapmaker
