program eep_cmd

  !MESA modules
  use utils_lib

  !local modules
  use iso_eep_support
  use iso_eep_color
  
  implicit none

  integer :: ierr = 0
  type(track) :: t
  character(len=file_path) :: BC_table_list = '', cmd_suffix = 'cmd'
  character(len=file_path) :: Cstar_table_list = '', phot_string
  logical :: set_fixed_Fe_div_H = .false., do_Cstars=.false.
  real(sp) :: extinction_Av=0.0, extinction_Rv=3.1, Fe_div_H = 0.0
  
  namelist /cmd_controls/ BC_table_list, extinction_Av, extinction_Rv, &
  Cstar_table_list, do_Cstars, cmd_suffix, Fe_div_H, set_fixed_Fe_div_H

  
  call cmd_init(ierr)

  if(ierr==0) call color_init(phot_string, BC_table_list, do_Cstars, Cstar_table_list, &
       set_fixed_Fe_div_H, Fe_div_H, ierr)
 
  if(ierr==0) call read_eep(t,full_path=.true.,append_eep=.false.)

  t% Av = extinction_Av
  t% cmd_suffix = trim(phot_string)

  if(ierr==0) call write_track_cmd_to_file(t)

contains

  subroutine cmd_init(ierr)
    integer, intent(out) :: ierr
    integer :: c, io
    character(len=32) :: arg
    ierr = 0
    c = command_argument_count()
    if(c<2)then
       write(0,*) ' eep_cmd '
       write(0,*) '    transform a theoretical eep file to a cmd      '
       write(0,*) '    usage: ./eep_cmd [phot string] [eep file] [Av] '
       write(0,*) '       [phot string] = UBVRIplus, etc.             '
       write(0,*) '       [eep file] = name eep file to transform     '
       write(0,*) '       [Av] optional argument, default is zero     '
       ierr=-1
    endif
    call get_command_argument(1,phot_string)
    call get_command_argument(2,t% filename)

    io=alloc_iounit(ierr)
    open(unit=io,file='input.nml', action='read', status='old', iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_track: problem reading input.nml '
       return
    endif
    read(io, nml=cmd_controls, iostat=ierr)
    close(io)
    
    if(c>2) then
       call get_command_argument(3,arg)
       read(arg,*) extinction_Av
    endif
  end subroutine cmd_init
  
end program eep_cmd
