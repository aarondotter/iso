program iso_interp_met

  !MESA modules
  use utils_lib

  !local modules
  use iso_eep_support
  use iso_interp_support

  implicit none

  type(isochrone_set), allocatable :: s(:)
  type(isochrone_set) :: t
  real(dp), allocatable :: Z_div_H(:)
  real(dp) :: new_Z_div_H
  integer :: ierr, n, i, i_Minit
  logical, parameter :: force_linear = .true.
  logical, parameter :: debug = .false.
  logical, parameter :: do_PAV = .true.
  
  call read_interp_input(ierr)
  if(ierr/=0) stop ' iso_interp_met: failed reading input list'

  do i=1,n
     call read_isochrone_file(s(i),ierr)
     if(ierr/=0)  stop ' iso_interp_met: failed in read_isochrone_file'
  enddo

  !need this for PAV
  i_Minit=locate_Minit(s(1)% iso(1))

  call consistency_check(s,ierr)
  if(ierr/=0) stop ' iso_interp_met: failed consistency check(s)'

  !if all checks pass, then proceed with interpolation
  call interpolate_Z(n,Z_div_H,new_Z_div_H,s,t,ierr)

  !write the new isochrone to file
  if(ierr==0) then
     call write_isochrones_to_file(t)
  endif

contains

  function locate_Minit(iso) result(i)
    type(isochrone), intent(in) :: iso
    integer :: i, j
    do j=1,iso% ncol
       if(trim(adjustl(iso% cols(j)% name))=='initial_mass')then
          i=j
          return
       endif
    enddo
    i=0
  end function locate_Minit

  subroutine read_interp_input(ierr)
    integer, intent(out) :: ierr
    integer :: io, i
    character(len=file_path) :: iso_list, arg
    ierr = 0

    if(command_argument_count() < 3) then
       ierr=-1
       write(0,*) 'iso_interp_met        '
       write(0,*) '   usage:             '
       write(0,*) ' ./iso_interp_met [list] [Z/H] [output]'
       return
    endif

    call get_command_argument(1,iso_list)
    call get_command_argument(2,arg)
    read(arg,*) new_Z_div_H
    call get_command_argument(3,t% filename)

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


  subroutine interpolate_Z(n,Z,newZ,s,t,ierr)
    integer, intent(in) :: n
    real(dp), intent(in) :: Z(n), newZ
    type(isochrone_set), intent(in) :: s(n)
    type(isochrone_set), intent(inout) :: t
    integer, intent(out) :: ierr
    integer :: i, loc, lo, hi, order
    loc=0; order=0; lo=0; hi=0

    ierr = 0

    do i=1,n-1
       if(newZ >= Z(i) .and. newZ < Z(i+1)) then
          loc=i
          exit
       endif
    enddo

    if(loc==0)then
       ierr=-1
       write(0,*) ' iso_interp_met: no extrapolation! '
       return
    endif

    !this results in either cubic (order=4) or linear (order=2)
    if(force_linear)then
       lo=loc
       hi=loc+1
    elseif(loc==n-1 .or. loc==1)then
       lo=loc
       hi=loc+1
    else
       lo=loc-1
       hi=loc+2
    endif

    order = hi - lo + 1  ! either 4, 3, or 2

    t% Fe_div_H = new_Z_div_H
    t% number_of_isochrones =s(lo)% number_of_isochrones
    t% MESA_revision_number = s(lo)% MESA_revision_number
    t% version_string = s(lo)% version_string
    !initial_Y, initial_Z, 
    t% v_div_vcrit = s(lo)% v_div_vcrit
    t% alpha_div_Fe = s(lo)% alpha_div_Fe
    allocate(t% iso(t% number_of_isochrones))
    t% iso(:)% age_scale = s(lo)% iso(:)% age_scale
    t% iso(:)% has_phase = s(lo)% iso(:)% has_phase
    t% iso(:)% age       = s(lo)% iso(:)% age
    t% iso(:)% ncol      = s(lo)% iso(:)% ncol
    t% iso(:)% alpha_div_Fe = t% alpha_div_Fe
    t% iso(:)% v_div_vcrit = t% v_div_vcrit
    t% iso(:)% Fe_div_H = t% Fe_div_H
    do i = 1, t% number_of_isochrones
       allocate(t% iso(i)% cols(t% iso(i)% ncol))
       t% iso(i)% cols(:)% name = s(lo)% iso(i)% cols(:)% name
       t% iso(i)% cols(:)% type = s(lo)% iso(i)% cols(:)% type
       t% iso(i)% cols(:)% loc = s(lo)% iso(i)% cols(:)% loc
    enddo

    call do_interp_all(order,Z(lo:hi),newZ,s(lo:hi),t)

    !PAV?
    if(do_PAV)then
       do i=1,t% number_of_isochrones
          call PAV(t% iso(i)% data(i_Minit,:))
       enddo
    endif
  end subroutine interpolate_Z

  subroutine do_interp_all(n,Z,newZ,s,t)
    integer, intent(in) :: n
    real(dp), intent(in) :: Z(n), newZ
    type(isochrone_set), intent(in) :: s(n)
    type(isochrone_set), intent(inout) :: t
    integer :: i, j, k, neep, eeps_lo(n), eeps_hi(n)
    integer :: eep, eep_lo, eep_hi, ioff(n), ncol
    logical :: eep_good(n)
    logical, pointer :: good(:)=>NULL()
    integer, pointer :: eep_index(:)=>NULL()
    real(dp) :: Yinit(1),Zinit(1),vv(1),ww(1),xx(1),yy(1)

    Zinit = 0d0
    Yinit = 0d0

    if(n==2)then
       xx(1) = s(1)% initial_Z
       yy(1) = s(2)% initial_Z
       Zinit = linear(1,Z,newZ,xx,yy)

       xx(1) = s(1)% initial_Y
       yy(1) = s(2)% initial_Y
       Yinit = linear(1,Z,newZ,xx,yy)       
    elseif(n==4)then      
       vv(1) = s(1)% initial_Z
       ww(1) = s(2)% initial_Z
       xx(1) = s(3)% initial_Z
       yy(1) = s(4)% initial_Z
       Zinit = cubic_pm(1,Z,newZ,vv,ww,xx,yy)

       vv(1) = s(1)% initial_Y
       ww(1) = s(2)% initial_Y
       xx(1) = s(3)% initial_Y
       yy(1) = s(4)% initial_Y
       Yinit = cubic_pm(1,Z,newZ,vv,ww,xx,yy)
    endif

    t% initial_Z = Zinit(1)
    t% initial_Y = Yinit(1)
    
    t% iso(:)% initial_Z = t% initial_Z
    t% iso(:)% initial_Y = t% initial_Y

    do i=1,n
       write(*,*) trim(s(i)% filename), Z(i)
    enddo

    do i=1,t% number_of_isochrones
       do j=1,n
          eeps_lo(j) = s(j)% iso(i)% eep(1)
          neep=s(j)% iso(i)% neep
          eeps_hi(j) = s(j)% iso(i)% eep(neep)
       enddo
       eep_lo = maxval(eeps_lo)
       eep_hi = minval(eeps_hi)

       allocate(good(eep_hi),eep_index(eep_hi))
       good=.false.
       eep_index=0
       neep=0
       eep_loop: do eep=eep_lo, eep_hi
          eep_good=.false.
          j_loop: do j=1,n
             k_loop: do k=1,s(j)% iso(i)% neep
                if(s(j)% iso(i)% eep(k)==eep)then
                   eep_good(j)=.true.
                   exit k_loop
                endif
             enddo k_loop
          enddo j_loop
          if(all(eep_good,1)) then
             neep = neep + 1
             good(eep) = .true.
             eep_index(eep) = neep
          elseif(debug)then
             write(*,*) ' bad EEP = ', eep
          endif
       enddo eep_loop

       t% iso(i)% neep = neep
       ncol = t% iso(i)% ncol
       allocate(t% iso(i)% eep(neep), t% iso(i)% data(ncol,neep))
       t% iso(i)% eep  = 0
       t% iso(i)% data = 0d0

       if(t% iso(i)% has_phase) then
          allocate(t% iso(i)% phase(neep))      
          t% iso(i)% phase = 0d0
       endif

       !loop over the EEPs
!$omp parallel do private(eep,ioff)
       do eep=eep_lo, eep_hi
          if(good(eep)) then
             t% iso(i)% eep(eep_index(eep)) = eep
             ioff=0
             do j=1,n
                do k=1,s(j)% iso(i)% neep
                   if(s(j)% iso(i)% eep(k) == eep)then
                      ioff(j)=k
                      exit
                   endif
                enddo
             enddo
             if(t%iso(i)% has_phase) t% iso(i)% phase(eep_index(eep)) = s(1)% iso(i)% phase(ioff(1))
             if(n==2)then 
                t% iso(i)% data(:,eep_index(eep)) = linear( ncol, Z, newZ, &
                     s(1)% iso(i)% data(:,ioff(1)), &
                     s(2)% iso(i)% data(:,ioff(2)) )
             else if(n==4) then
                t% iso(i)% data(:,eep_index(eep)) = cubic_pm( ncol, Z, newZ, &
                     s(1)% iso(i)% data(:,ioff(1)), &
                     s(2)% iso(i)% data(:,ioff(2)), &
                     s(3)% iso(i)% data(:,ioff(3)), &
                     s(4)% iso(i)% data(:,ioff(4)) )
             endif
          endif
       enddo
!$omp end parallel do
       deallocate(good)
       nullify(good)
    enddo
  end subroutine do_interp_all

end program iso_interp_met
