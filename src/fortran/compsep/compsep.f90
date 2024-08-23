module compsep
  use healpix_types
  implicit none
  
  type comm_data_set
     character(len=512)  :: label
     integer(i4b)        :: nside, lmax
     real(dp), allocatable, dimension(:) :: map
     real(dp), allocatable, dimension(:) :: rms
     real(dp), allocatable, dimension(:) :: b_l
  end type comm_data_set
  
  integer(i4b) :: numband
  type(comm_data_set), allocatable, dimension(:) :: data
  
contains

end module compsep
