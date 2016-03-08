program make_cmd

  !MESA modules
  use utils_lib

  !local modules
  use iso_eep_support
  use iso_color

  implicit none
 
  integer :: ierr
  type(isochrone_set) :: s
  logical :: do_Cstars = .false.
  character(len=file_path) :: BC_table_list = '', cmd_suffix = 'cmd'
  character(len=file_path) :: Cstar_table_list = ''
  real(sp) :: extinction_Av=0.0, extinction_Rv=0.0

  namelist /cmd_controls/ extinction_Av, extinction_Rv, BC_table_list, &
  Cstar_table_list, do_Cstars, cmd_suffix

  call cmd_init(ierr)
  if(ierr==0)then

     call read_isochrone_file(s,ierr)

     call write_cmds_to_file(s)

  endif


contains

  subroutine cmd_init(ierr)
    integer, intent(out) :: ierr
    integer :: io 

    if(command_argument_count()<1)then
       write(*,*) ' make_cmd:   '
       write(*,*) '   usage: ./make_cmd [isochrone file]'
       write(*,*) '     all other options set through cmd_controls in input.nml'
       ierr=-1
       return
    endif

    call get_command_argument(1,s% filename)
    
    io=alloc_iounit(ierr)
    open(io,file='input.nml',action='read',status='old', iostat=ierr)
    read(io,nml=cmd_controls, iostat=ierr)
    close(io)

    s% iso(:)% Av = extinction_Av
    s% iso(:)% Rv = extinction_Rv
    s% cmd_suffix  = cmd_suffix

    call iso_color_init(BC_table_list,do_Cstars,Cstar_table_list,ierr)

  end subroutine cmd_init

end program make_cmd
