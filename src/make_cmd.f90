program make_cmd

  use iso_eep_support
  use iso_eep_color

  implicit none

  integer :: ierr
  type(isochrone_set) :: s
  logical, parameter :: do_timing = .false.
  logical :: do_Cstars = .false.
  character(len=file_path) :: BC_table_list = 'bc_table.list', cmd_suffix = ''
  character(len=file_path) :: Cstar_table_list = 'cstar_table.list'
  logical :: set_fixed_Fe_div_H = .false.
  real(sp) :: extinction_Av, extinction_Rv, Fe_div_H
  integer :: count_rate, time(3)
  ierr=0
  
  if(do_timing) call system_clock(time(1),count_rate)

  call cmd_init(ierr)

  s% cmd_suffix = cmd_suffix

  if(do_timing) call system_clock(time(2),count_rate)

  if(ierr==0) call write_cmds_to_file(s)

  if(do_timing) call system_clock(time(3),count_rate)

  if(do_timing) call report_timing

contains

  subroutine report_timing
    real(dp) :: t(3),c
    t=real(time,kind=dp)
    c=real(count_rate,kind=dp)
    t=t/c
    write(0,*) ' Total execution time: ', t(3)-t(1)
    write(0,*) ' Time to init,read   : ', t(2)-t(1)
    write(0,*) ' Time to write isos  : ', t(3)-t(2)
  end subroutine report_timing

  subroutine cmd_init(ierr)
    integer, intent(out) :: ierr
    integer :: i, j, n_arg
    character(len=file_path) :: phot_string, arg, option, result
    
    Fe_div_H = 0.0
    extinction_Av = 0.0
    extinction_Rv = 0.0
    
    if(command_argument_count()<1)then
       write(*,*) ' make_cmd:   '
       write(*,*) '   usage: ./make_cmd phot_string isochrone_file [Av] '
       write(*,*) '     phot_string = UBVRIplus, etc.                         '
       write(*,*) '     isochrone_file = name of isochrone file to transform  '
       write(*,*) '     OPTIONAL -                                            '
       write(*,*) '     [Av] extinction in V band                             '
       ierr=-1
       return
    endif

    call get_command_argument(1,phot_string)
    call get_command_argument(2,s% filename)

    if(cmd_suffix == '') cmd_suffix=phot_string

    do i=1, command_argument_count()
       call get_command_argument(i,arg)
    enddo
    
    !process command arguments
    n_arg = command_argument_count()
    if(n_arg > 2) then
       do i=3,n_arg
          call get_command_argument(i,arg)
          j=index(arg,'=')
          option=arg(1:j-1)
          result=arg(j+1:)
          if(trim(option)=='Av')then
             read(result,*) extinction_Av
          endif
       enddo
    endif
    
    if(ierr==0) call read_isochrone_file(s,ierr)

    if(ierr/=0) write(*,*) 'read_isochrone_file: ierr = ', ierr
    
    s% Av = extinction_Av
    s% Rv = extinction_Rv
    do i=1,s% number_of_isochrones
       s% iso(i)% Av = s% Av
       s% iso(i)% Rv = s% Rv
    enddo

    call color_init(phot_string,BC_table_list,do_Cstars,Cstar_table_list, &
         set_fixed_Fe_div_H,Fe_div_H,ierr)

  end subroutine cmd_init

end program make_cmd
