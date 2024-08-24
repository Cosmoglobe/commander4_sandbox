subroutine tod_init_band(band, nscan) bind(c, name="tod_init_band")
  use data_mod
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in) :: band, nscan
  
  print *, "tod_init_band called with ", band, nscan

  data(band)%tod%nscan = nscan
  allocate(data(band)%tod%scans(nscan))
  
end subroutine tod_init_band

subroutine tod_init_scan(band, scan, ntod, d, pix) bind(c, name="tod_init_scan")
  use data_mod
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in) :: band, scan, ntod
  real(c_float),      intent(in) :: d(ntod)
  integer(c_int32_t), intent(in) :: pix(ntod)
  
  print *, "tod_init_scan called with ", band, scan, ntod

  data(band)%tod%scans(scan)%ntod = ntod
  allocate(data(band)%tod%scans(scan)%d(ntod), data(band)%tod%scans(scan)%pix(ntod))
  data(band)%tod%scans(scan)%d   = d
  data(band)%tod%scans(scan)%pix = pix

end subroutine tod_init_scan

subroutine tod_estimate_sigma0(band, scan, signal, npix) bind(c, name="tod_estimate_sigma0")
  use data_mod
  use iso_c_binding
  implicit none
  integer(c_int64_t), intent(in)  :: band, scan, npix
  real(c_double),     intent(in)  :: signal(npix)

  integer(c_int64_t) :: i, pix
  real(c_double) :: sigmasq
  
  sigmasq = 0.d0
  do i = 1, data(band)%tod%scans(scan)%ntod
     pix     = data(band)%tod%scans(scan)%pix(i)
     sigmasq = sigmasq + (data(band)%tod%scans(scan)%d(i) -  signal(pix))**2
  end do
  data(band)%tod%scans(scan)%sigma0 = sqrt(sigmasq / (data(band)%tod%scans(scan)%ntod-1))
  
end subroutine tod_estimate_sigma0

subroutine tod_mapmaker(band) bind(c, name="tod_mapmaker")
  use data_mod
  use iso_c_binding
  implicit none
  integer(c_int64_t),  intent(in) :: band

  integer(4) :: i, pix
  real(c_double), allocatable, dimension(:) :: A, b

  allocate(A(0:data(band)%npix-1), b(0:data(band)%npix-1))
  A = 0.d0
  b = 0.d0
  ! Collect contributions from all samples
  do i = 1, data(band)%tod%nscan
     pix    = data(band)%tod%scans(i)%pix(i)
     A(pix) = A(pix) + 1.d0
     b(pix) = b(pix) + data(band)%tod%scans(i)%d(pix)
  end do

  ! Solve for map and rms
  data(band)%map = b/A
  data(band)%rms = 1.d0/sqrt(A)
    
end subroutine tod_mapmaker
