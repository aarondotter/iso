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
  logical :: set_fixed_Fe_div_H = .false.
  real(sp) :: extinction_Av=0.0, extinction_Rv=0.0, Fe_div_H = 0.0

  namelist /cmd_controls/ BC_table_list, extinction_Av, extinction_Rv, &
  Cstar_table_list, do_Cstars, cmd_suffix, Fe_div_H, set_fixed_Fe_div_H

  call cmd_init(ierr)
  if(ierr==0) call read_isochrone_file(s,ierr)
  if(ierr==0) call write_cmds_to_file(s)

contains

  subroutine cmd_init(ierr)
    integer, intent(out) :: ierr
    integer :: io , i
    character(len=32) :: phot_string, arg

    if(command_argument_count()<1)then
       write(*,*) ' make_cmd:   '
       write(*,*) '   usage: ./make_cmd [phot string] [isochrone file] [Av]'
       write(*,*) '     [phot string] = UBVRIJHKs, etc.                    '
       write(*,*) '     [isochrone file] = name of isochrone file to transform'
       write(*,*) '     [Av] optional argument; if not set then take value from input.nml'
       write(*,*) '     all other options set through cmd_controls in input.nml'
       ierr=-1
       return
    endif
    
    call get_command_argument(1,phot_string)
    call get_command_argument(2,s% filename)
    
    io=alloc_iounit(ierr)
    open(io,file='input.nml',action='read',status='old', iostat=ierr)
    read(io,nml=cmd_controls, iostat=ierr)
    close(io)
    call free_iounit(io)


    if(command_argument_count()>2)then
       call get_command_argument(3, arg)
       read(arg,*) extinction_Av
    endif

    s% Av = extinction_Av
    s% Rv = extinction_Rv
    do i=1,s% number_of_isochrones
       s% iso(i)% Av = s% Av
       s% iso(i)% Rv = s% Rv
    enddo
    s% cmd_suffix  = phot_string

    call iso_color_init(phot_string,BC_table_list,do_Cstars,Cstar_table_list, &
         set_fixed_Fe_div_H,Fe_div_H,ierr)

  end subroutine cmd_init

end program make_cmd
