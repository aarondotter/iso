program make_cmd

  !MESA modules
  use utils_lib

  !local modules
  use iso_eep_support
  use iso_eep_color

  implicit none
 
  integer :: ierr
  type(isochrone_set) :: s
  logical, parameter :: do_timing = .false.
  logical :: do_Cstars = .false.
  character(len=file_path) :: BC_table_list = '', cmd_suffix = 'cmd'
  character(len=file_path) :: Cstar_table_list = ''
  logical :: set_fixed_Fe_div_H = .false.
  logical :: include_gravity_darkening = .false.
  real(sp) :: extinction_Av=0.0, extinction_Rv=0.0, Fe_div_H = 0.0
  integer :: count_rate, time(4)

  namelist /cmd_controls/ BC_table_list, extinction_Av, extinction_Rv, &
       Cstar_table_list, do_Cstars, cmd_suffix, Fe_div_H, set_fixed_Fe_div_H, &
       include_gravity_darkening
  

  if(do_timing) call system_clock(time(1),count_rate)
  call cmd_init(ierr)

  s% cmd_suffix = cmd_suffix !override value in bin, if present

  if(do_timing) call system_clock(time(2),count_rate)
  if(ierr==0) call write_cmds_to_file(s)

  if(do_timing) call system_clock(time(3),count_rate)

  if(do_timing) call report_timing
contains
  subroutine report_timing
    real(dp) :: t(4),c
    t=real(time,kind=dp)
    c=real(count_rate,kind=dp)
    t=t/c
    write(0,*) ' Total execution time: ', (t(3)-t(1))
    write(0,*) ' Time to init,read   : ', (t(2)-t(1))
    write(0,*) ' Time to write isos  : ', (t(3)-t(2))
  end subroutine report_timing

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

  if(ierr==0) call read_isochrone_file(s,ierr)
    
    s% Av = extinction_Av
    s% Rv = extinction_Rv
    s% include_gravity_darkening = include_gravity_darkening
    do i=1,s% number_of_isochrones
       s% iso(i)% Av = s% Av
       s% iso(i)% Rv = s% Rv
       s% iso(i)% include_gravity_darkening = s% include_gravity_darkening
    enddo
   
    call color_init(phot_string,BC_table_list,do_Cstars,Cstar_table_list, &
         set_fixed_Fe_div_H,Fe_div_H,ierr)

  end subroutine cmd_init

end program make_cmd
