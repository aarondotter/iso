program charlies_angel

  use iso_eep_support

  implicit none
  
  type(track) :: m
  type(isochrone_set) :: s
  integer :: ierr
  character(len=file_path) :: output

  if(command_argument_count()<3)then
     write(0,*) "*** Charlie's Angel ***"
     write(0,*) '   ./charlies_angel [eep file] [iso file] [output]'
     write(0,*) '   appends points to [iso file] from [eep file] at each age'
     write(0,*) '   [output] is the new isochrone file name'
     stop '   badness: not enough command arguments provided'
  endif
  
  call get_command_argument(1,m% filename)
  call read_eep(m,.true.)

  call get_command_argument(2,s% filename)
  call read_isochrone_file(s,ierr)
  if(ierr/=0) stop 'failed reading isochrone file'
  
  call get_command_argument(3,output)

  call interpolate_eep_onto_iso(m,s,output) 
  
contains
  
  subroutine interpolate_eep_onto_iso(t,s,output)
    type(track), intent(in) :: t
    type(isochrone_set), intent(in) :: s
    character(len=file_path) :: output
    integer :: i, ierr
    type(isochrone_set) :: q
    real(dp), allocatable :: res(:)

    q% filename = output

    q% MESA_revision_number = s% MESA_revision_number
    q% version_string = s% version_string
    q% number_of_isochrones = s% number_of_isochrones
    q% Fe_div_H = s% Fe_div_H
    q% alpha_div_Fe = s% alpha_div_Fe
    q% v_div_vcrit =  s% v_div_vcrit
    q% initial_Z = s% initial_Z
    q% initial_Y = s% initial_Y

    allocate(q% iso(q% number_of_isochrones))
    allocate(res(s% iso(1)% ncol)) !not safe!
    
    q% iso(:)% has_phase = s% iso(:)% has_phase
    q% iso(:)% age_scale = s% iso(:)% age_scale
    q% iso(:)% age       = s% iso(:)% age
    q% iso(:)% initial_Z = q% initial_Z
    q% iso(:)% initial_Y = q% initial_Y
    q% iso(:)% Fe_div_H = q% Fe_div_H
    q% iso(:)% alpha_div_Fe = q% alpha_div_Fe
    q% iso(:)% v_div_vcrit = q% v_div_vcrit

    do i=1,s% number_of_isochrones
       res = 0d0
       q% iso(i)% ncol = s% iso(i)% ncol
       allocate(q% iso(i)% cols(q% iso(i)% ncol))
       q% iso(i)% cols = s% iso(i)% cols
       q% iso(i)% age = s% iso(i)% age
       call interpolate_age(t,pow10(q% iso(i)% age),res,ierr)
       if(ierr/=0) stop 99 
       
       q% iso(i)% neep = s% iso(i)% neep + 1

       allocate(q% iso(i)% eep(q% iso(i)% neep))
       allocate(q% iso(i)% data(q% iso(i)% ncol, q% iso(i)% neep))
       if(q% iso(i)% has_phase) allocate( q% iso(i)% phase(q% iso(i)% neep))
       
       q% iso(i)% eep(1) = s% iso(i)% eep(1) -1 
       q% iso(i)% eep(2:q% iso(i)% neep) = s% iso(i)% eep

       q% iso(i)% data(:,1) = res
       q% iso(i)% data(:,2:q% iso(i)% neep) = s% iso(i)% data

    enddo

    call write_isochrones_to_file(q)

  end subroutine interpolate_eep_onto_iso

    
  subroutine interpolate_age(t,age,res,ierr) 
    type(track), intent(in) :: t
    real(dp), intent(in) :: age
    real(dp), intent(out) :: res(:)
    integer, intent(out) :: ierr
    integer :: i, i_age, lo, hi
    real(dp) :: alfa, beta

    res = 0d0; lo=0; hi=0; ierr=0

    if(size(res) /= t% ncol) then
       ierr=-1
       return
    endif
    
    i_age = 0

    do i=1,t% ncol
       if(trim(adjustl(t% cols(i)% name))=='star_age') then
          i_age = i
          exit
       endif
    enddo

    if(i_age == 0) then
       ierr=-1
       return
    endif

    do i=1,t% ntrack-1
       if(age >= t% tr(i_age,i) .and. age < t% tr(i_age,i+1))then
          lo=i
          hi=i+1
          exit
       endif
    enddo

    if(lo==0)then
       ierr=-1
       return
    endif

    !interpolate
    beta = (age - t% tr(i_age,lo))/(t% tr(i_age,hi)-t% tr(i_age,lo))
    alfa = 1d0 - beta
    !switch age to mass for isochrone
    res = alfa*t% tr(:,lo) + beta*t% tr(:,hi)
    res(i_age) = t% initial_mass
  end subroutine interpolate_age
    
end program charlies_angel
