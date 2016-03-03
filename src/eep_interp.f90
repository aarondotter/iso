program eep_interp

  !local modules
  use iso_eep_support

  implicit none

  character(len=file_path) :: infile(4), outfile
  type(track) :: s(4), t
  real(dp) :: met(4) = [-1.5d0, -1.0d0, 0d0, 0.3d0], new_met = -0.05
  integer :: i

  do i=1,4
     call get_command_argument(i,s(i)% filename)
     call read_iso(s(i))
  enddo

  t% neep = s(1)% neep


end program eep_interp
