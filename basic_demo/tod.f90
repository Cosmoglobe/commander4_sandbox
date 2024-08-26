module tod_mod
use data_mod
use iso_c_binding
implicit none

contains

subroutine tod_init_band(band, nscan)
  integer(c_int64_t), intent(in) :: band, nscan
  
  print *, "tod_init_band called with ", band, nscan

  allocate(data(band)%tod%scans(nscan))
end subroutine tod_init_band

subroutine tod_init_scan(band, scan, d, pix)
  integer(c_int64_t), intent(in) :: band, scan
  real(c_float), dimension(0:), intent(in) :: d
  integer(c_int32_t), dimension(0:), intent(in) :: pix
  integer :: ntod

  ntod = size(d,1)
  if (ntod /= size(pix,1)) then
    print *, "ERROR"
    stop 1
  endif
  print *, "tod_init_scan called with ", band, scan

  allocate(data(band)%tod%scans(scan)%d(0:ntod-1), data(band)%tod%scans(scan)%pix(0:ntod-1))
  data(band)%tod%scans(scan)%d   = d
  data(band)%tod%scans(scan)%pix = pix
end subroutine tod_init_scan

subroutine tod_estimate_sigma0(band, scan, signal)
  integer(c_int64_t), intent(in)  :: band, scan
  real(c_double), dimension(0:), intent(in)  :: signal

  integer(c_int64_t) :: i, pix, ntod
  real(c_double) :: sigmasq

  ntod = size(data(band)%tod%scans(scan)%pix,1)
  sigmasq = 0.d0
  do i = 0, ntod-1
     pix     = data(band)%tod%scans(scan)%pix(i)
     sigmasq = sigmasq + (data(band)%tod%scans(scan)%d(i) -  signal(pix))**2
  end do
  data(band)%tod%scans(scan)%sigma0 = sqrt(sigmasq / (ntod-1))
end subroutine tod_estimate_sigma0

subroutine tod_mapmaker(band, map_ifc, rms_ifc)
  use data_mod
  use iso_c_binding
  implicit none
  integer(c_int64_t),  intent(in) :: band
  real(c_double), intent(inout), dimension(0:) :: map_ifc, rms_ifc


  integer(8) :: i, j, scan, pix
  real(c_double), allocatable, dimension(:) :: A, b

  allocate(A(0:data(band)%npix-1), b(0:data(band)%npix-1))
  A = 0.d0
  b = 0.d0
  ! Collect contributions from all samples

  do scan = 1, size(data(band)%tod%scans)
     do j = 0, size(data(band)%tod%scans(scan)%pix) - 1
       pix    = data(band)%tod%scans(scan)%pix(j)
       A(pix) = A(pix) + 1.d0
       b(pix) = b(pix) + data(band)%tod%scans(scan)%d(j)
     end do
  end do

  ! Solve for map and rms
  data(band)%map = b/A
  data(band)%rms = 1.d0/sqrt(A)
  map_ifc = data(band)%map
  rms_ifc = data(band)%rms
end subroutine tod_mapmaker

end module tod_mod

! foreign language interface routines ... simple forwarders to module members

subroutine tod_init_band_ifc(band, nscan) bind(c, name="tod_init_band_ifc")
  use tod_mod
  implicit none
  integer(c_int64_t), intent(in), value :: band, nscan
  
  call tod_init_band(band, nscan)
end subroutine tod_init_band_ifc

subroutine tod_init_scan_ifc(band, scan, ntod, d, pix) bind(c, name="tod_init_scan_ifc")
  use tod_mod
  implicit none
  integer(c_int64_t), intent(in), value :: band, scan, ntod
  real(c_float),      intent(in) :: d(ntod)
  integer(c_int32_t), intent(in) :: pix(ntod)

  call tod_init_scan(band, scan, d, pix)
end subroutine tod_init_scan_ifc

subroutine tod_estimate_sigma0_ifc(band, scan, signal, npix) bind(c, name="tod_estimate_sigma0_ifc")
  use tod_mod
  implicit none
  integer(c_int64_t), intent(in), value  :: band, scan, npix
  real(c_double),     intent(in)  :: signal(0:npix-1)

  call tod_estimate_sigma0(band, scan, signal)
end subroutine tod_estimate_sigma0_ifc

subroutine tod_mapmaker_ifc(band, map_ifc, rms_ifc, npix) bind(c, name="tod_mapmaker_ifc")
  use tod_mod
  implicit none
  integer(c_int64_t),  intent(in), value :: band, npix
  real(c_double), intent(inout), dimension(0:npix-1) :: map_ifc
  real(c_double), intent(inout), dimension(0:npix-1) :: rms_ifc


  call tod_mapmaker(band, map_ifc, rms_ifc)
end subroutine tod_mapmaker_ifc
