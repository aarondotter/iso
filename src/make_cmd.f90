program make_cmd

  !MESA modules
  use utils_lib
  use crlibm_lib

  !local modules
  use iso_eep_support
  use iso_eep_color

  implicit none

  integer :: ierr
  type(isochrone_set) :: s
  logical, parameter :: do_timing = .false.
  logical :: do_Cstars = .false., do_spots=.false.
  character(len=file_path) :: BC_table_list = 'bc_table.list', cmd_suffix = ''
  character(len=file_path) :: Cstar_table_list = 'cstar_table.list'
  logical :: set_fixed_Fe_div_H = .false.
  real(sp) :: extinction_Av=0.0, extinction_Rv=0.0, Fe_div_H = 0.0
  integer :: count_rate, time(4)
  real(dp) :: spot_beta=0.0d0, spot_gamma=0.0d0, Teff_spot_full_on=0.0d0, Teff_spot_full_off=0.0d0

!!$  namelist /cmd_controls/ BC_table_list, extinction_Av, extinction_Rv, &
!!$       Cstar_table_list, do_Cstars, cmd_suffix, Fe_div_H, set_fixed_Fe_div_H, &
!!$       Teff_spot_full_on, Teff_spot_full_off

  Teff_spot_full_on = 5.0d3
  Teff_spot_full_off= 7.0d3

  if(do_timing) call system_clock(time(1),count_rate)
  call crlibm_init
  call cmd_init(ierr)

  s% cmd_suffix = cmd_suffix

  if(do_timing) call system_clock(time(2),count_rate)

  if(ierr==0.and.do_spots)then
     call spotify(s,spot_beta,spot_gamma)
  endif

  if(do_timing) call system_clock(time(3),count_rate)
  if(ierr==0) call write_cmds_to_file(s)

  if(do_timing) call system_clock(time(4),count_rate)

  if(do_timing) call report_timing

contains

  subroutine report_timing
    real(dp) :: t(4),c
    t=real(time,kind=dp)
    c=real(count_rate,kind=dp)
    t=t/c
    write(0,*) ' Total execution time: ', t(4)-t(1)
    write(0,*) ' Time to init,read   : ', t(2)-t(1)
    write(0,*) ' Time to add spots   : ', t(3)-t(2)
    write(0,*) ' Time to write isos  : ', t(4)-t(3)
  end subroutine report_timing

  subroutine cmd_init(ierr)
    integer, intent(out) :: ierr
    integer :: i, j, n_arg
    character(len=32) :: phot_string, arg, option, result

    if(command_argument_count()<1)then
       write(*,*) ' make_cmd:   '
       write(*,*) '   usage: ./make_cmd phot_string isochrone_file [Av] [beta] [gamma]'
       write(*,*) '     phot_string = UBVRIplus, etc.                         '
       write(*,*) '     isochrone_file = name of isochrone file to transform  '
       write(*,*) '     OPTIONAL -                                            '
       write(*,*) '     [Av] extinction in V band                             '
       write(*,*) '     [beta] governs fraction of surface covered by spots   '
       write(*,*) '     [gamma] governs ratio of T_spot to T_effective        '
       ierr=-1
       return
    endif

    call get_command_argument(1,phot_string)
    call get_command_argument(2,s% filename)

    if(cmd_suffix == '') cmd_suffix=phot_string
    
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
          elseif(trim(option)=='beta')then
             read(result,*) spot_beta
          elseif(trim(option)=='gamma')then
             read(result,*) spot_gamma 
          endif
       enddo
    endif

    do_spots = (spot_beta > 0.0d0).and.(spot_gamma < 1.0d0)
    
    if(ierr==0) call read_isochrone_file(s,ierr)

    s% Av = extinction_Av
    s% Rv = extinction_Rv
    do i=1,s% number_of_isochrones
       s% iso(i)% Av = s% Av
       s% iso(i)% Rv = s% Rv
    enddo

    call color_init(phot_string,BC_table_list,do_Cstars,Cstar_table_list, &
         set_fixed_Fe_div_H,Fe_div_H,ierr)

  end subroutine cmd_init

  subroutine spotify(s,beta0,gamma)
    !based on Chabrier et al. (2007) with minor extension*
    !beta = fraction of surface covered by spots
    !*gamma = temperature ratio, Tspot/Teff (Chabrier+2007 set gamma=0)
    !Teff and radius are changed without modifying luminosity
    type(isochrone_set), intent(inout) :: s
    real(dp), intent(in) :: beta0, gamma
    real(dp) :: a, Tfac, alfa, beta, R, Teff, g, log_Teff
    real(dp) :: log_Teff_spot_full_on, log_Teff_spot_full_off, gnew, Rnew, Tnew
    integer :: i, j, kT=0, kR=0, kg=0

    !first locate log_Teff and log_L
    do i=1,s% iso(1)% ncol
       if(trim(adjustl(s% iso(1)% cols(i)% name)) == 'log_Teff') then
          kT=i
       else if(trim(adjustl(s% iso(1)% cols(i)% name))== 'log_R') then
          kR=i
       else if(trim(adjustl(s% iso(1)% cols(i)% name))=='log_g') then
          kg=i
       endif
    enddo

    log_Teff_spot_full_on = log10_cr(Teff_spot_full_on)
    log_Teff_spot_full_off = log10_cr(Teff_spot_full_off)
    
    do i=1, s% number_of_isochrones
       do j=1, s% iso(i)% neep
          log_Teff = s% iso(i)% data(kT,j)
          if(log_Teff > log_Teff_spot_full_off) then
             a = 0.0d0
          elseif(log_Teff < log_Teff_spot_full_on) then
             a = 1.0d0
          else
             a = (log_Teff_spot_full_off - log_Teff)/(log_Teff_spot_full_off-log_Teff_spot_full_on)
          endif

          Tfac = 0.5d0*(1d0 - cospi_cr(a))

          !combine all multiplicative factors into one
          beta = Tfac*beta0 !this is the fraction of the spots
          alfa = 1d0 - beta !this is the fraction of the unspotted star

          g = exp10_cr(s% iso(i)% data(kg,j))
          R = exp10_cr(s% iso(i)% data(kR,j))
          Teff = exp10_cr(log_Teff)
          Tnew = pow_cr(alfa*powi_cr(Teff,4) + beta*powi_cr(Teff*gamma,4),2.5d-1)
          Rnew = R*powi_cr(Teff/Tnew,2)
          gnew = g*powi_cr(R/Rnew,2)
          s% iso(i)% data(kR,j) = log10_cr(Rnew)
          s% iso(i)% data(kT,j) = log10_cr(Tnew)
          s% iso(i)% data(kg,j) = log10_cr(gnew)
       enddo
    enddo
  end subroutine spotify

end program make_cmd
