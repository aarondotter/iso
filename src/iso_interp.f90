program iso_interp

  !MESA modules
  use const_def, only: dp
  use utils_lib
  use interp_1d_def
  use interp_1d_lib

  !local modules
  use iso_eep_support
  use iso_color

  implicit none

  type(isochrone_set), allocatable :: s(:)
  type(isochrone_set) :: t
  real(dp), allocatable :: Z_div_H(:)
  real(dp) :: new_Z_div_H
  integer :: ierr = 0, n, i
  logical :: do_colors

  call read_interp_input(ierr)
  if(ierr/=0) stop ' iso_interp: failed reading input list'

  do i=1,n
     call read_isochrone_file(s(i),ierr)
     if(ierr/=0)  stop ' iso_interp: failed in read_isochrone_file'
  enddo

  call consistency_check(s,ierr)
  if(ierr/=0) stop ' iso_interp: failed consistency check(s)'

  !if all checks pass, then proceed with interpolation
  call interpolate(n,Z_div_H,new_Z_div_H,s,t,ierr)

  !write the new isochrone to file
  if(ierr==0) then
     call write_isochrones_to_file(t)
     if(do_colors) call write_cmds_to_file(t)
  endif

contains


  subroutine read_interp_input(ierr)
    integer, intent(out) :: ierr
    integer :: io, i
    character(len=file_path) :: iso_list, arg
    character(len=file_path) :: bc_list, cstar_list
    ierr = 0
    if(command_argument_count() < 3) then
       ierr=-1
       write(0,*) 'iso_interp            '
       write(0,*) '   usage:             '
       write(0,*) ' ./iso_interp [list] [Z/H] [output] [Av]'
       write(0,*) '  Av is optional;                       '
       write(0,*) '   will generate a .cmd file if included'
       return
    else if(command_argument_count() > 3) then
       do_colors=.true.
       t% cmd_suffix = 'cmd'
       bc_list = 'bc_table.list'
       cstar_list = 'cstar_table.list'
       call iso_color_init(bc_list,.false.,cstar_list,ierr)
       if(ierr/=0) then 
          write(0,'(a)') ' problem reading BC_table_list = ', trim(bc_list)
          do_colors=.false.
       endif
    endif

    call get_command_argument(1,iso_list)
    call get_command_argument(2,arg)
    read(arg,*) new_Z_div_H
    call get_command_argument(3,t% filename)
    if(do_colors)then
       call get_command_argument(4, arg)
       read(arg,*) t% Av
       t% Rv = 3.1
    endif

    io=alloc_iounit(ierr)
    open(io,file=trim(iso_list),action='read',status='old')
    read(io,*) n
    allocate(s(n),Z_div_H(n))
    do i=1,n
       read(io,'(f5.2,1x,a)') Z_div_H(i), s(i)% filename
    enddo
    close(io)
    call free_iounit(io)
  end subroutine read_interp_input


  subroutine consistency_check(s,ierr)
    type(isochrone_set), intent(in) :: s(:)
    integer, intent(out) :: ierr
    integer :: n, i

    ierr = 0

    n=size(s)

    if(any(s(:)% number_of_isochrones /= s(1)% number_of_isochrones))then
       write(0,*) ' iso_interp : failed consistency check - inconsistent number of isochrones '
       ierr=-1
       do i=1,n
          write(0,*) trim(s(i)% filename), s(i)% version_number
       enddo
    endif

    if(any(s(:)% version_number /= s(1)% version_number))then
       write(0,*) ' iso_interp : failed consistency check - inconsistent version number '
       ierr=-1
       do i=1,n
          write(0,*) trim(s(i)% filename), s(i)% number_of_isochrones
       enddo
    endif

    !more checks here as needed

  end subroutine consistency_check


  subroutine interpolate(n,Z,newZ,s,t,ierr)
    integer, intent(in) :: n
    real(dp), intent(in) :: Z(n), newZ
    type(isochrone_set), intent(in) :: s(n)
    type(isochrone_set), intent(inout) :: t
    integer, intent(out) :: ierr
    integer :: loc=0, lo=0, hi=0, order=0

    ierr = 0

    do loc=1,n-1
       if(newZ >= Z_div_H(loc) .and. newZ < Z_div_H(loc+1)) exit
    enddo

    if(loc==0)then
       ierr=-1
       write(0,*) ' iso_interp: no extrapolation! '
       return
    endif

    !try for cubic interpolation but do something else if not
    lo=max(loc-1,1)
    hi=min(loc+2,n)

    order = hi - lo + 1  ! either 4, 3, or 2

    t% number_of_isochrones =s(lo)% number_of_isochrones
    t% version_number = s(lo)% version_number
    allocate(t% iso(t% number_of_isochrones))
    t% iso(:)% age_scale = s(lo)% iso(:)% age_scale
    t% iso(:)% has_phase = s(lo)% iso(:)% has_phase
    t% iso(:)% age       = s(lo)% iso(:)% age
    t% iso(:)% ncol      = s(lo)% iso(:)% ncol
    do i = 1, t% number_of_isochrones
       allocate(t% iso(i)% cols(t% iso(i)% ncol))
       t% iso(i)% cols(:)% name = s(lo)% iso(i)% cols(:)% name
       t% iso(i)% cols(:)% type = s(lo)% iso(i)% cols(:)% type
       t% iso(i)% cols(:)% loc = s(lo)% iso(i)% cols(:)% loc
    enddo
    call do_interp(order,Z(lo:hi),newZ,s(lo:hi),t)
  end subroutine interpolate


  subroutine do_interp(n,Z,newZ,s,t)
    integer, intent(in) :: n
    real(dp), intent(in) :: Z(n), newZ
    type(isochrone_set), intent(in) :: s(n)
    type(isochrone_set), intent(inout) :: t
    integer :: i, j, k, neep, eeps_lo(n), eeps_hi(n)
    integer :: eep, eep_lo, eep_hi, ioff(n), ncol
    real(dp) :: f(n)

    call coeff(Z,f,newZ,n)

    do i=1,n
       write(*,*) trim(s(i)% filename), Z(i), f(i)
    enddo

    do i=1,t% number_of_isochrones
       do j=1,n
          eeps_lo(j) = s(j)% iso(i)% eep(1)
          neep=s(j)% iso(i)% neep
          eeps_hi(j) = s(j)% iso(i)% eep(neep)
       enddo
       eep_lo = maxval(eeps_lo)
       eep_hi = minval(eeps_hi)
       t% iso(i)% neep = eep_hi - eep_lo + 1
       neep = t% iso(i)% neep
       ncol = t% iso(i)% ncol
       allocate(t% iso(i)% eep(neep), t% iso(i)% data(ncol,neep))
       if(t% iso(i)% has_phase) allocate(t% iso(i)% phase(neep))      
       t% iso(i)% data = 0d0

       !the EEP offsets are for alignment
       do k=1,n
          ioff(k) = eep_lo - s(k)% iso(i)% eep(1)  
       enddo

       !loop over the EEPs
       do eep=1,neep
          t% iso(i)% eep(eep) = eep_lo + eep - 1
          t% iso(i)% phase(eep) = s(1)% iso(i)% phase(eep+ioff(1))
          do k=1,n
             t% iso(i)% data(:,eep) = t% iso(i)% data(:,eep) + f(k)*s(k)% iso(i)% data(:,eep+ioff(k))
          enddo
       enddo
    enddo
  end subroutine do_interp


  subroutine coeff(a,b,x,n)
    ! {a} are the tabulated values for use in interpolation
    ! {b} are coefficients of the interpolating polynomial
    !  x  is the abscissa to be interpolated
    !  n  is the number of points to be used, interpolating polynomial
    !     has order n-1 
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n), x
    real(dp), intent(out) :: b(n)
    integer :: i,j
    b=1d0
    do i=1,n
       do j=1,n
          if(j/=i) b(i)=b(i)*(x-a(j))/(a(i)-a(j))
       enddo
    enddo
  end subroutine coeff


end program iso_interp
