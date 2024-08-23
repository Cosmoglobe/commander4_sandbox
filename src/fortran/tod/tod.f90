module compsep
  use healpix_types
  implicit none

  type :: scan
     integer(i4b) :: ntod
     real(sp),           allocatable, dimension(:)     :: tod            ! Detector values in time domain
     integer(i4b),       allocatable, dimension(:)     :: pix            ! pointer array of pixels 
  end type scan

  type :: tod
     character(len=512) :: label
     character(len=512) :: filelist
     integer(i4b)       :: nscan
     integer(i4b)       :: nside
     type(scan), allocatable, dimension(:) :: scans    ! Array of all scans
  end type comm_tod

contains

end module compsep
