program iso_interp
  
  !MESA modules
  use const_def, only: dp
  use utils_lib
  use interp_1d_def
  use interp_1d_lib

  !local modules
  use iso_eep_support

  implicit none

  type(isochrone_set), allocatable :: s(:)
  integer :: ierr, n, i

  ierr = 0
  
  call read_interp_input(ierr)
  if(ierr/=0) stop

  do i=1,n
     call read_isochrone_file(s(i),ierr)
  enddo

  write(*,*) ' ierr = ', ierr

contains

  subroutine read_interp_input(ierr)
    integer, intent(out) :: ierr
    integer :: io, i
    character(len=file_path) :: iso_list
    ierr = 0
    if(command_argument_count() < 1) then
       ierr=-1
       write(0,*) 'iso_interp            '
       write(0,*) '   usage:             '
       write(0,*) ' ./iso_interp [list]  '
       return
    endif

    call get_command_argument(1,iso_list)

    io=alloc_iounit(ierr)
    open(io,file=trim(iso_list),action='read',status='old')
    read(io,*) n
    allocate(s(n))
    do i=1,n
       read(io,'(a)') s(i)% filename
    enddo
    close(io)
    call free_iounit(io)
  end subroutine read_interp_input


end program iso_interp
