program iso_interp_age

  !MESA modules
  use const_def, only: dp
  use utils_lib

  !local modules
  use iso_eep_support
  use iso_interp_support
  use iso_color

  implicit none

  type(isochrone_set) :: old, new
  integer, parameter :: single=0, range=1, list=2
  integer :: ierr
  logical, parameter :: debug = .false.
  logical, parameter :: force_linear=.true.
  
  call read_interp_age_input(ierr)
  if(ierr/=0) stop ' iso_interp_age: failed reading input list'

  call interpolate_age(old,new,ierr)

  !write the new isochrone to file
  if(ierr==0) call write_isochrones_to_file(new)

contains

  subroutine read_interp_age_input(ierr)
    integer, intent(out) :: ierr
    integer :: age_scale, age_type, number_of_ages, i
    character(len=file_path) :: input, output, arg
    real(dp) :: age_range_lo, age_range_hi, delta
    real(dp), allocatable :: age(:)
    ierr = 0

    if(command_argument_count() < 4) then
       ierr=-1
       write(0,*) 'iso_interp_age        '
       write(0,*) '   usage:             '
       write(0,*) ' ./iso_interp_age [input] [output] [age_scale] [age_type] [age(s)]'
       write(0,*) '     [input] = input file name   '
       write(0,*) '     [output] = output file name (can be same as input)'
       write(0,*) '     [age_scale] = integer, 0 for linear or 1 for log10'
       write(0,*) '     [age_type]  = integer, 0 for single, 1 for range, 2 for list'
       write(0,*)
       write(0,*) '     for age_type = single, enter the age'
       write(0,*)
       write(0,*) '     for age_type = range, enter the number, minimum and maximum age'


       write(0,*) '     [age(s)] = real               '
       return
    endif

    call get_command_argument(1,input)
    call get_command_argument(2,output)
    call get_command_argument(3,arg)
    read(arg,*) age_scale
    call get_command_argument(4,arg)
    read(arg,*) age_type
    if(age_type==single)then
       number_of_ages = 1
       allocate(age(1))
       call get_command_argument(5,arg)
       read(arg,*) age(1)
    else
       call get_command_argument(5,arg)
       read(arg,*) number_of_ages

       allocate(age(number_of_ages))
       if(age_type==range) then
          call get_command_argument(6,arg)
          read(arg,*) age_range_lo
          call get_command_argument(7,arg)
          read(arg,*) age_range_hi
          delta = (age_range_hi - age_range_lo)/real(number_of_ages-1,kind=dp)
          do i=1,number_of_ages
             age(i) = age_range_lo + delta*real(i-1,kind=dp)
          enddo
       else if(age_type==list)then
          do i=1,number_of_ages
             call get_command_argument(5+i,arg)
             read(arg,*) age(i)
          enddo
       endif
    endif

    !now set up the isochrone types...
    old% filename = trim(adjustl(input))
    new% filename = trim(adjustl(output))

    call read_isochrone_file(old,ierr)
    if(ierr/=0) return    

    new% version_number = old% version_number
    new% number_of_isochrones = number_of_ages

    allocate(new% iso(new% number_of_isochrones))
    new% iso(:)% has_phase = old% iso(1)% has_phase


    do i=1,new% number_of_isochrones
       new% iso(i)% age_scale = age_scale
       if(new% iso(i)% age_scale == age_scale_linear)then
          new% iso(i)% age = log10(age(i))
       elseif(new% iso(i)% age_scale == age_scale_log10)then
          new% iso(i)% age = age(i)
       endif
    enddo

    if(.false.)then
       do i=1,new% number_of_isochrones
          write(*,*) new% iso(i)% age
       enddo
       
       do i=1,old% number_of_isochrones
          write(*,*) old% iso(i)% age
       enddo
       
       ierr=-1
    endif

  end subroutine read_interp_age_input

  subroutine interpolate_age(old,new,ierr)
    type(isochrone_set), intent(in) :: old
    type(isochrone_set), intent(inout) :: new
    integer, intent(out) :: ierr
    integer :: i, j, n, loc, lo, hi, order
    loc=0; order=0; lo=0; hi=0

    ierr = 0

    !we have two sets of ages, new and old

    !first make sure they are compatible
    if(minval(new% iso(:)% age) < minval(old% iso(:)% age))then
       write(0,*) ' iso_interp_age error: new age minimum is less than old age minimum '
       ierr=-1
    else if(maxval(new% iso(:)% age) > maxval(old% iso(:)% age))then
       write(0,*) ' iso_interp_age error: new age maximum is greater than old age maximum '
       ierr=-1
    endif

    if(ierr/=0) return

    n = old% number_of_isochrones

    do i = 1, new% number_of_isochrones

       do j = 1, n-1
          if(new% iso(i)% age > old% iso(j)% age .and. new% iso(i)% age <= old% iso(j+1)% age)then
             loc=j
             exit
          endif
       enddo

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
       
       order = hi - lo + 1  ! either 4 or 2 

       new% iso(i)% has_phase = old% iso(lo)% has_phase
       new% iso(i)% ncol = old% iso(lo)% ncol
       allocate(new% iso(i)% cols(new% iso(i)% ncol))
       do j = 1, new% number_of_isochrones
          new% iso(i)% cols(:)% name = old% iso(lo)% cols(:)% name
          new% iso(i)% cols(:)% type = old% iso(lo)% cols(:)% type
          new% iso(i)% cols(:)% loc  = old% iso(lo)% cols(:)% loc
       enddo
       call do_interp_one(order,old% iso(lo:hi),new% iso(i))
    enddo

  end subroutine interpolate_age


  subroutine do_interp_one(n,s,t)
    integer, intent(in) :: n
    type(isochrone), intent(in) :: s(n)
    type(isochrone), intent(inout) :: t
    integer :: j, k, neep, eeps_lo(n), eeps_hi(n)
    integer :: eep, eep_lo, eep_hi, ioff(n), ncol
    real(dp) :: ages(n), new_age
    logical :: eep_good(n)
    logical, allocatable :: good(:)
    integer, allocatable :: eep_index(:)
    new_age = t% age

    do j=1,n
       eeps_lo(j) = s(j)% eep(1)
       neep=s(j)% neep
       eeps_hi(j) = s(j)% eep(neep)
       ages(j) = s(j)% age
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
          k_loop: do k=1,s(j)% neep
             if(s(j)% eep(k)==eep)then
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

    t% neep = neep
    ncol = t% ncol
    allocate(t% eep(neep), t% data(ncol,neep))
    t% eep  = 0
    t% data = 0d0

    if(t% has_phase) then
       allocate(t% phase(neep))      
       t% phase = 0d0
    endif

    !loop over the EEPs
    !$omp parallel do private(eep)
    do eep=eep_lo,eep_hi
       if(good(eep))then
          t% eep(eep_index(eep)) = eep
          ioff=0
          do j=1,n
             do k=1,s(j)% neep
                if(s(j)% eep(k) == eep)then
                   ioff(j)=k
                   exit
                endif
             enddo
          enddo
          if(t% has_phase) t% phase(eep_index(eep)) = s(1)% phase(ioff(1))
          if(n==2)then 
             t% data(:,eep_index(eep)) = linear(ncol, ages, new_age, &
                  s(1)% data(:,ioff(1)), &
                  s(2)% data(:,ioff(2)))
          else if(n==4) then 
             t% data(:,eep_index(eep)) = cubic_pm( ncol, ages, new_age, &
                  s(1)% data(:,ioff(1)), &
                  s(2)% data(:,ioff(2)), &
                  s(3)% data(:,ioff(3)), &
                  s(4)% data(:,ioff(4)))
          endif
       endif
    enddo
    !$omp end parallel do

  end subroutine do_interp_one

end program iso_interp_age
