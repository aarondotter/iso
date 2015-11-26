program make_isochrone

  !MESA modules
  use const_def, only: dp, sp
  use utils_lib
  use interp_1d_def
  use interp_1d_lib

  !local modules
  use iso_eep_support
  use iso_color

  implicit none

  character(len=file_path) :: input_file, history_columns_list
  integer :: i, ierr, io, j, niso, ntrk, ngood, first, prev, i_Minit
  type(track), allocatable :: s(:), t(:), q
  type(isochrone_set) :: set

  logical :: use_double_eep
  integer, parameter :: piecewise_monotonic = 4
  integer, parameter :: age_scale_linear = 0
  integer, parameter :: age_scale_log10  = 1
  integer :: age_scale

  !default some namelist parameters
  logical :: iso_debug = .false.
  logical :: do_tracks = .false.
  logical :: do_isochrones = .true.
  logical :: do_smooth = .true.
  logical :: do_PAV = .true.
  logical :: do_colors = .false.
  logical :: do_Cstars = .false.
  character(len=file_path) :: BC_table_list='', Cstar_table_list=''
  character(len=file_path) :: cmd_suffix = '.cmd'
  character(len=file_path) :: track_filename = 'out.trk'
  real(sp) :: extinction_Av = 0.0, extinction_Rv = 0.0
  real(dp) :: track_initial_mass = 1d0

  namelist /iso_controls/ iso_debug, do_tracks, do_isochrones, do_smooth, &
       do_PAV, do_colors, do_Cstars, BC_table_list, Cstar_table_list, cmd_suffix, &
       extinction_Av, extinction_Rv, track_initial_mass, track_filename
  

  ierr=0
  if(command_argument_count()<1) then
     write(*,*) '   make_iso                   '
     write(*,*) '   usage: ./make_iso <input>  '
     stop       '   no command line argument   '
  endif

  !begin
  call read_iso_input(ierr)
  if(ierr/=0) write(*,*) '  read_iso_input: ierr = ', ierr

  !read eep files to fill s()
  first=ntrk
  prev=first
  ngood=0

  do i=1,ntrk
     call read_eep(s(i))
     if(s(i)% ignore) cycle
     if(i<first) then 
        first=i             !this marks the first valid eep file
        prev=first
     elseif(prev>first)then
        prev=first
     endif
     if(iso_debug) write(*,'(a50,f8.2,99i8)') &
          trim(s(i)% filename), s(i)% initial_mass, s(i)% eep
     !check for monotonic mass, consistent phase info and version number
     if(ngood > 1)then
        if( s(i)% initial_mass < s(prev)% initial_mass ) &
             stop ' make_iso: masses out of order'
        if( s(i)% has_phase.neqv.s(prev)% has_phase ) &
             stop ' make_iso: inconsistent phase info in tracks'
        if( s(i)% version_number /= s(prev)% version_number )&
             stop ' make_iso: inconsistent version number in tracks'
     endif
     ngood=ngood+1
  enddo

  allocate(t(ngood))
  j=0
  do i=1,ntrk
     if(.not.s(i)% ignore) then
        j=j+1
        t(j)=s(i)
     endif
  enddo
  deallocate(s)

  !above checks pass => these are safe assignments
  set% iso(:)% has_phase = t(1)% has_phase
  set% version_number = t(1)% version_number

  !interpolate a new track . . .
  if(do_tracks)then
     allocate(q)
     q% initial_mass = track_initial_mass
     q% filename = track_filename
     write(*,*) ' call interpolate_track'
     call interpolate_track(t,q)
     call write_track(q)
     deallocate(q)
  endif

  !create isochrones 
  if(do_isochrones)then
     do i=1,niso
        call do_isochrone_for_age(t,set% iso(i))
     enddo
     call write_isochrones_to_file(set)
     if(do_colors) call write_cmds_to_file(set)
  endif

  !all done.
  deallocate(t,cols)

contains

  subroutine do_isochrone_for_age(s,iso)
    type(track), intent(in) :: s(:)
    type(isochrone), intent(inout) :: iso
    integer :: eep, hi, ierr, index, interp_method, j, k, l, lo, max_eep, n, pass
    integer :: max_neep_low, max_neep_high
    character(len=col_width) :: mass_age_string = 'mass from age'
    real(dp) :: age, mass, min_age, max_age
    real(dp), pointer :: ages(:)=>NULL(), masses(:)=>NULL()
    real(dp), allocatable :: result1(:,:), result2(:,:), mass_tmp(:)
    logical, allocatable :: skip(:,:)
    integer, allocatable :: valid(:), count(:)
    real(dp), parameter :: age_delta = 0.75d0

    ierr = 0

    !initialize some quantities
    n = size(s) ! is the number of tracks
    max_eep = maxval(s(:)% ntrack) !is the largest number of EEPs in any track
    allocate(skip(n,max_eep),count(max_eep))

    !set method and options for interpolation
    ! - for interp_m3: average, quartic, or super_bee
    !interp_method = average
    ! - for interp_pm: piecewise_monotonic
    interp_method = piecewise_monotonic

    if(age_scale==age_scale_linear)then
       age = log10(iso% age)
    elseif(age_scale==age_scale_log10)then
       age = iso% age
    endif
    mass=0d0

    !this is temporary storage for the isochrone data:
    !result1 stores the data for all EEPs, valid tells
    !which ones are good and will be returned via the iso 
    !derived type
    allocate(result1(ncol,max_eep),result2(ncol,max_eep),valid(max_eep))
    result1 = 0d0
    result2 = 0d0
    valid = 0
    skip = .false.
    count = 0

    !determine the largest number of EEPs in tracks of different types
    max_neep_low = 0
    max_neep_high = 0
    do k=1,n
       if(s(k)% star_type < star_high_mass)then
          max_neep_low = max(max_neep_low, s(k)% neep)
       else ! high-mass star
          max_neep_high = max(max_neep_high, s(k)% neep)
       endif
    enddo

    !now check each track to make sure it is complete for its typ
    do k=1,n
       if(s(k)% star_type == star_high_mass .and. s(k)% neep < max_neep_high) then
          skip(k,:) = .true.
       else if(s(k)% initial_mass >= 0.6d0 .and. s(k)% star_type == star_low_mass .and. s(k)% neep < max_neep_low) then
          skip(k,:) = .true.
       endif
    enddo

    eep_loop1: do eep=1,max_eep

       !determine tracks for which the ith eep is defined
       !the skip logical array determines whether or not a given
       !track will be included in the ensuing interpolation steps.
       !count keeps track of how many tracks will be used. if
       !fewer than 2 tracks satisfy the condition, skip the EEP
       do k=1,n
          max_age=age + age_delta
          min_age=age - age_delta

          if(s(k)% eep(1) > eep .or. s(k)% eep(s(k)% neep) < eep ) then
             skip(k,eep) = .true.
          else if(.not.use_double_eep .and. log10(s(k)% tr(i_age,eep)) > max_age ) then
             skip(k,eep) = .true.
          else if(.not.use_double_eep .and. log10(s(k)% tr(i_age,eep)) < min_age ) then
             skip(k,eep) = .true.
          endif
       enddo

       !this loop attempts to pick out non-monotonic points

       if(.true.) then !top-down

          if(.not.use_double_eep)then
             do k=n,2,-1
                if(skip(k,eep)) cycle
                do l=k-1,1
                   if( skip(l,eep)) cycle
                   if( s(k)% tr(i_age,eep) > s(l)% tr(i_age,eep) ) skip(l,eep) = .true.
                enddo
             enddo
          endif

       else !bottom-up

          if(.not.use_double_eep)then
             do k=1,n-1
                if(skip(k,eep)) cycle
                do l=k+1,n
                   if( skip(l,eep)) cycle
                   if( s(k)% tr(i_age,eep) < s(l)% tr(i_age,eep) ) skip(l,eep) = .true.
                enddo
             enddo
          endif
          
       endif

       do k=1,n
          if(.not.skip(k,eep)) count(eep)=count(eep)+1
       enddo



       if(iso_debug) write(*,*) '  EEP, count, n = ', eep, count(eep), n
       if(count(eep) < 2)then
          if(iso_debug) write(*,*) 'not enough eeps to interpolate'
          cycle eep_loop1
       endif

       !allocate ages and masses that pass the above test
       !I use pointers here because they can be allocated
       !and deallocated to the proper size for each EEP
       if(associated(ages)) then
          deallocate(ages)
          nullify(ages)
       endif
       if(associated(masses)) then
          deallocate(masses)
          nullify(masses)
       endif
       allocate(ages(count(eep)),masses(count(eep)))

       !this step fills the masses and ages arrays that are
       !used as the basis for interpolation
       j=0
       do k=1,n
          if(.not.skip(k,eep))then
             j=j+1
             ages(j) = log10(s(k)% tr(i_age,eep))
             masses(j) = s(k)% initial_mass
          endif
       enddo

       if(do_smooth.and.all(masses>0.5,dim=1)) call smooth(masses,ages)

       !check to see if the input age is found within the
       !current set of ages. if not, skip to the next EEP.
       j = binary_search( count(eep), ages, 1, age)
       if( j < 1 .or. j > count(eep)-1 ) cycle eep_loop1

       !check to see if masses and ages are monotonic
       !if not, then interpolation will fail
       if(.not.monotonic(masses)) then
          write(*,*) ' masses not monotonic in do_isochrone_for_age: ', age
          stop 99
       endif

       if(iso_debug) then 
          write(*,*) ' j, count = ', j, count(eep)
          write(*,*) skip(:,eep)
          do k=1,count(eep)
             write(*,*) masses(k), ages(k), age
          enddo
       endif

       !this block checks between two tracks at the current EEP to see if the 
       !age lies in between the two. if it does, then it outputs a linear mass 
       !interpolation. it is checking to see if this happens more than once.
       !j gives the number of times that this happens for a given EEP.
       if(use_double_eep)then
          pass=0 
          k_loop: do k=1,count(eep)-1

             if(.not.skip(k,eep))then

                if( ages(k) > ages(k+1) .and. ages(k) >= age .and. ages(k+1) < age )then
                   pass = pass + 1
                   if(pass==1)then
                      call monotonic_mass_range(ages,k,lo,hi)
                      if(iso_debug) write(*,*) mass, masses(lo), masses(hi)

                      mass=interp_x_from_y(ages(lo:hi),masses(lo:hi),age,mass_age_string,ierr)
                      if(ierr/=0)then
                         write(0,*) ' age interpolation failed for eep, pass = ', eep, pass
                         cycle eep_loop1
                      endif

                      if(iso_debug) write(*,'(2i5,f9.5,3f14.9)') eep, pass, age, mass, masses(lo), masses(hi)

                      result1(i_Minit,eep) = mass
                      valid(eep)=1
                   else if(pass==2)then
                      call monotonic_mass_range(ages,k,lo,hi)
                      if(iso_debug) write(*,*) mass, masses(lo), masses(hi)

                      mass=interp_x_from_y(ages(lo:hi),masses(lo:hi),age,mass_age_string,ierr)
                      if(ierr/=0)then
                         write(0,*) ' age interpolation failed for eep, pass = ', eep, pass
                         cycle eep_loop1
                      endif

                      if(iso_debug) write(*,'(2i5,f9.5,3f14.9)') eep, pass, age, mass, masses(lo), masses(hi)

                      result2(i_Minit,eep) = mass
                      valid(eep)=2
                   endif

                endif

             endif

          enddo k_loop

       else ! single EEP case; this is the original method.
          ! interpolate in age to find the EEP's initial mass
          index = i_Minit     ! special case for iso_intepolate

          mass = iso_interpolate( eep, interp_method, n, &
               s, skip(:,eep), count(eep), index, ages, age, mass_age_string, ierr)
          if(ierr/=0)then
             write(0,*) ' interpolation failed in age->mass'
             cycle eep_loop1
          endif
          result1(i_Minit,eep) = mass
          valid(eep)=1
       endif

    enddo eep_loop1


    !PAV forces monotonicity in the masses
    if(do_PAV .and. .not.use_double_eep)then
       j=0
       allocate(mass_tmp(max_eep))
       do eep=1,max_eep
          if(valid(eep)>0) then 
             j = j+1
             mass_tmp(j) = result1(i_Minit,eep)
          endif
       enddo

       call PAV(mass_tmp(1:j))

       j=0
       do eep=1,max_eep
          if(valid(eep)>0)then
             j=j+1
             result1(i_Minit,eep) = mass_tmp(j)
          endif
       enddo
       deallocate(mass_tmp)
    endif


    eep_loop2: do eep=1,max_eep

       if(count(eep)<2.or.valid(eep)<1) cycle eep_loop2

       !allocate ages and masses that pass the above test
       !I use pointers here because they can be allocated
       !and deallocated to the proper size for each EEP
       if(associated(ages)) then
          deallocate(ages)
          nullify(ages)
       endif
       if(associated(masses)) then
          deallocate(masses)
          nullify(masses)
       endif
       allocate(ages(count(eep)),masses(count(eep)))

       !this step fills the masses and ages arrays that are
       !used as the basis for interpolation
       j=0
       do k=1,n
          if(.not.skip(k,eep))then
             j=j+1
             ages(j) = log10(s(k)% tr(i_age,eep))
             masses(j) = s(k)% initial_mass
          endif
       enddo

       if(do_smooth.and.all(masses>0.5,dim=1)) call smooth(masses,ages)

       !this block checks between two tracks at the current EEP to see if the 
       !age lies in between the two. if it does, then it outputs a linear mass 
       !interpolation. it is checking to see if this happens more than once.
       !j gives the number of times that this happens for a given EEP.
       if(use_double_eep)then
          pass=0 
          j_loop: do j=1,count(eep)-1

             if(.not.skip(j,eep))then

                if( ages(j) > ages(j+1) .and. ages(j) >= age .and. ages(j+1) < age )then
                   pass = pass + 1
                   if(pass==1)then

                      do index = 2, ncol
                         mass = result1(i_Minit,eep)
                         result1(index,eep) = iso_interpolate(eep, interp_method, n, &
                              s, skip(:,eep), count(eep), index, masses, mass, cols(index), ierr)
                         if(ierr/=0) then
                            write(0,*) ' mass interpolation failed for index = ', &
                                 trim(cols(index))
                            valid(eep)=0
                            cycle eep_loop2
                         endif
                      enddo
                      valid(eep)=1

                   else if(pass==2)then

                      do index = 2, ncol
                         mass = result2(i_Minit,eep)
                         result2(index,eep) = iso_interpolate(eep, interp_method, n, &
                              s, skip(:,eep), count(eep), index, masses, mass, cols(index), ierr)
                         if(ierr/=0) then
                            write(0,*) ' mass interpolation failed for index = ', &
                                 trim(cols(index))
                            valid(eep)=0
                            cycle eep_loop2
                         endif
                      enddo
                      valid(eep)=2
                   endif

                endif

             endif

          enddo j_loop


       else ! single EEP case; this is the original method.

          do index = 2, ncol
             mass = result1(i_Minit,eep)
             result1(index,eep) = iso_interpolate(eep, interp_method, n, &
                  s, skip(:,eep), count(eep), index, masses, mass, cols(index), ierr)
             if(ierr/=0) then
                write(0,*) ' mass interpolation failed for index = ', trim(cols(index))
                valid(eep)=0
                cycle eep_loop2
             endif
          enddo
          valid(eep) = 1

       endif
    enddo eep_loop2

    deallocate(masses,ages)


    !now result1 and valid are full for all EEPs,
    !we can pass the data to the iso derived type
    iso% ncol = ncol
    iso% neep = sum(valid)
    allocate(iso% cols(iso% ncol))
    allocate(iso% data(iso% ncol, iso% neep), iso% eep(iso% neep))
    if(iso% has_phase) allocate(iso% phase(iso% neep))

    iso% cols = cols

    j=0
    do eep=1,max_eep
       if(valid(eep)>0) then
          j=j+1
          iso% eep(j) = eep
          iso% data(:,j) = result1(:,eep)
          if(iso% has_phase)then
             do k=1,n-1
                if(iso% data(i_Minit,j) < s(k)% initial_mass) exit
             enddo
             iso% phase(j) = s(k)% phase(min(s(k)% ntrack,eep))
          endif
       endif

       if(valid(eep)==2) then
          j=j+1
          iso% eep(j) = eep
          iso% data(:,j) = result2(:,eep)
          if(iso% has_phase) then
             do k=1,n-1
                if(iso% data(i_Minit,j) < s(k)% initial_mass) exit
             enddo
             iso% phase(j) = s(k)% phase(min(s(k)% ntrack,eep))
          endif
       endif

    enddo

    !all done
    deallocate(result1,result2,skip,valid)

  end subroutine do_isochrone_for_age

  subroutine monotonic_mass_range(ages,k,lo,hi)
    !this subroutine determines the upper and lower limits that
    !should be used in the EEP interpolation for the case of 
    !double EEPs; default result is lo=1, hi=size(ages)
    real(dp), intent(in) :: ages(:)
    integer, intent(in) :: k
    integer, intent(out) :: lo, hi
    integer :: j,n
    n=size(ages)
    !find lower limit:
    lo=1
    lower_loop: do j=k,1,-1
       if(ages(j+1) < ages(j)) then
          lo=j
       else
          exit lower_loop
       endif
    enddo lower_loop
    !find upper limit
    hi=n
    upper_loop: do j=lo+1,n
       if(ages(j-1) > ages(j)) then
          hi = j
       else
          exit upper_loop
       endif
    enddo upper_loop
  end subroutine monotonic_mass_range


  real(dp) function interp_x_from_y(x,y,x_in,label,ierr)
    real(dp), intent(in) :: x(:), y(:), x_in
    integer, parameter :: n_new = 1
    real(dp) :: x0(n_new), y0(n_new)
    real(dp), pointer :: work1(:)
    character(len=col_width), intent(in) :: label
    integer, intent(out) :: ierr
    integer :: n_old
    ierr=0
    if(size(x)/=size(y))then
       ierr=-1
       interp_x_from_y = 0d0
       return
    endif
    n_old = size(x)
    x0(n_new) = x_in
    allocate(work1(n_old*pm_work_size))
    call interpolate_vector_pm( n_old, x, n_new, x0, y, y0, work1, label, ierr )
    deallocate(work1)
    interp_x_from_y = y0(n_new)
  end function interp_x_from_y

  real(dp) function iso_interpolate(eep, method, n, s, skip, count, index, x_array, x, label, ierr)
    integer, intent(in) :: eep, method, n
    type(track), intent(in) :: s(n)
    logical, intent(in) :: skip(n)
    integer, intent(in) :: count, index
    real(dp), intent(in) :: x_array(count), x
    character(len=col_width), intent(in) :: label
    integer, intent(out) :: ierr
    integer :: j,k
    integer, parameter :: nwork = max(mp_work_size,pm_work_size)
    real(dp) :: y 
    real(dp), target :: f_ary(4*count), work_ary(count*nwork)
    real(dp), pointer :: f1(:)=>NULL(), f(:,:)=>NULL(), work(:)=>NULL()

    ierr = 0
    iso_interpolate = 0d0

    !set up pointers for interpolation
    f1 => f_ary
    f(1:4,1:count) => f1(1:4*count)
    work => work_ary

    !check again that enough tracks are defined to interpolate
    if(count < 2) then
       ierr=-1
       return
    endif

    !fill the interpolant f(1,:)
    j=0
    do k=1,n
       if(.not.skip(k))then
          j=j+1
          if(index == i_Minit)then
             f(1,j) = s(k)% initial_mass
          else
             f(1,j) = s(k)% tr(index,eep)
          endif
       endif
    enddo

    !perform the interpolation, y~f(x), using the input method
    if(method == piecewise_monotonic)then
       call interp_pm(x_array,count,f1,nwork,work,label,ierr)
    else
       call interp_m3(x_array,count,f1,method,nwork,work,label,ierr)
    endif

    call interp_value(x_array,count,f1,x,y,ierr)      

    if(is_bad_num(x) .or. is_bad_num(y))then
       write(0,*) ' eep = ', eep
       write(0,*) ' x = ', x
       write(0,*) ' y = ', y
       do k=1,n
          write(0,*) x_array(k), f(1,k), skip(k)
       enddo
       ierr=-1
    endif
    iso_interpolate = y

  end function iso_interpolate


  !writes a series of n isochrones to filename
  subroutine write_isochrones_to_file(set)
    type(isochrone_set), intent(in) :: set
    integer :: i, ierr, io,n
    io=alloc_iounit(ierr)
    n=size(set% iso)
    write(0,*) ' isochrone output file = ', trim(set% filename)
    open(io,file=trim(set% filename),action='write',status='unknown',iostat=ierr)
    write(io,'(a25,i5)') '# number of isochrones = ', n
    write(io,'(a25,i5)') '# MESA version number  = ', set% version_number
    do i=1,n
       call write_isochrone_to_file(io,set% iso(i))
       if(i<n) write(io,*)
       if(i<n) write(io,*)
    enddo
    close(io)
    call free_iounit(io)
  end subroutine write_isochrones_to_file

  !writes one age isochrone to the open io unit
  subroutine write_isochrone_to_file(io,iso)
    integer, intent(in) :: io
    type(isochrone), intent(in) :: iso
    if(iso% has_phase)then
       call write_isochrone_to_file_phase(io,iso)
    else
       call write_isochrone_to_file_orig(io,iso)
    endif
  end subroutine write_isochrone_to_file

  subroutine write_isochrone_to_file_orig(io,iso)
    integer, intent(in) :: io
    type(isochrone), intent(in) :: iso
    integer :: my_ncol
    my_ncol = iso% ncol + 2 !add two for eep and age
    write(io,'(a25,2i5)') '# number of EEPs, cols = ', iso% neep, my_ncol
    write(io,'(a1,i4,299i32)') '#    ', (i,i=1,my_ncol)
    if(age_scale==age_scale_log10)then
       write(io,'(a5,299a32)') '# EEP', 'log10_isochrone_age_yr', adjustr(iso% cols)
    elseif(age_scale==age_scale_linear)then
       write(io,'(a5,299a32)') '# EEP', 'isochrone_age_yr', adjustr(iso% cols)
    endif
    do i=1,iso% neep
       write(io,'(i5,299(1pes32.16e3))') iso% eep(i), iso% age, iso% data(:,i)
    enddo
  end subroutine write_isochrone_to_file_orig

  subroutine write_isochrone_to_file_phase(io,iso)
    integer, intent(in) :: io
    type(isochrone), intent(in) :: iso
    integer :: my_ncol
    my_ncol = iso% ncol + 3 !add three for eep, phase, and age
    write(io,'(a25,2i5)') '# number of EEPs, cols = ', iso% neep, my_ncol
    write(io,'(a1,i4,299i32)') '#    ', (i,i=1,my_ncol)
    if(age_scale==age_scale_log10)then
       write(io,'(a5,299a32)') '# EEP', 'log10_isochrone_age_yr', adjustr(iso% cols)
    elseif(age_scale==age_scale_linear)then
       write(io,'(a5,299a32)') '# EEP', 'isochrone_age_yr', adjustr(iso% cols)
    endif
    do i=1,iso% neep
       write(io,'(i5,299(1pes32.16e3))') iso% eep(i), iso% age, &
            iso% data(:,i), real(iso% phase(i),kind=dp)
    enddo
  end subroutine write_isochrone_to_file_phase

  subroutine write_cmds_to_file(set)
    type(isochrone_set), intent(inout) :: set
    character(len=256) :: output
    integer :: i, io, n, ierr
    if(set% cmd_suffix/='') output = trim(set% filename) // '.' // trim(set% cmd_suffix)
    n=size(set% iso)
    write(0,*) ' cmd output file = ', trim(output)
    io = alloc_iounit(ierr)
    if(ierr/=0) return
    open(io,file=trim(output),action='write',status='unknown',iostat=ierr)
    if(ierr/=0) return
    write(io,'(a25,i5)')    '# number of isochrones = ', n
    write(io,'(a25,i5)')    '# MESA version number  = ', set% version_number
    write(io,'(a25,2f6.3)') '# values of Av and Rv  = ', set% Av
    do i=1,n
       set% iso(i)% Av = set% Av
       call write_cmd_to_file(io, set% iso(i))
       if(i<n) write(io,*)
       if(i<n) write(io,*)
    enddo
    close(io)
    call free_iounit(io)
  end subroutine write_cmds_to_file

  subroutine write_cmd_to_file(io,iso)
    integer, intent(in) :: io
    type(isochrone), intent(inout) :: iso
    integer :: i, iT, ig, iL, ierr
    ierr=0; iT=0; ig=0; iL=0

    do i=1, iso% ncol
       if(trim(iso% cols(i)) == 'log_Teff') then
          iT=i
       else if(trim(iso% cols(i)) == 'log_g') then
          ig=i
       else if(trim(iso% cols(i))== 'log_L') then
          iL=i
       endif
    enddo

    call get_mags(iso,do_Cstars,ierr)
    if(ierr/=0) write(0,*) ' problem in get_mags '

    write(io,'(a25,2i5)') '# number of EEPs, cols = ', iso% neep, iso% nfil + 6
    write(io,'(a1,i4,5i32,299(17x,i3))') '#    ', (i,i=1,iso% nfil+6)

    if(age_scale==age_scale_linear)then
       write(io,'(a5,5a32,299a20)') '# EEP', 'isochrone_age_yr', 'initial_mass', 'log_Teff', &
            'log_g', 'log_L', adjustr(iso% labels)
    elseif(age_scale==age_scale_log10)then
       write(io,'(a5,5a32,299a20)') '# EEP', 'log10_isochrone_age_yr', 'initial_mass', 'log_Teff', &
            'log_g', 'log_L', adjustr(iso% labels)
    endif

    do i = 1,iso% neep
       write(io,'(i5,5(1pes32.16e3),299(0pf20.6))') iso% eep(i), iso% age, iso% data(i_Minit,i), &
            iso% data(iT,i), iso% data(ig,i), iso% data(iL,i), iso% mags(:,i)
    enddo
  end subroutine write_cmd_to_file

  !takes a set of EEP-defined tracks and interpolates a new
  !track for the desired initial mass
  subroutine interpolate_track(a,b)
    type(track), intent(in) :: a(:)
    type(track), intent(inout) :: b
    real(dp) :: f(3), dm, x(4), y(4)
    real(dp), pointer :: initial_mass(:)=>NULL() !(n)
    integer :: i, j, k, m, mlo, mhi, n

    n = size(a)

    allocate(initial_mass(n))
    initial_mass = a(:)% initial_mass
    m = binary_search(n, initial_mass, 1, b% initial_mass)
    write(*,*) b% initial_mass
    write(*,*) initial_mass(m:m+1)

    mlo = min(max(1,m-1),n-3)
    mhi = max(min(m+2,n),4)
    m = mlo+1

    write(*,*) '   mlo, m, mhi = ', mlo, m, mhi
    write(*,*) initial_mass(mlo:mhi)
    write(*,*) a(mlo:mhi)% neep

    b% neep = minval(a(mlo:mhi)% neep)

    write(*,*) b% neep

    b% ntrack = a(m)% eep(b% neep)

    write(*,*) b% ntrack

    b% ncol = a(1)% ncol
    b% has_phase = a(1)% has_phase
    allocate(b% cols(b% ncol))
    b% cols = a(1)% cols
    b% version_number = a(1)% version_number
    allocate(b% eep(b% neep))
    allocate(b% tr(b% ncol, b% ntrack))
    allocate(b% dist(b% ntrack))
    if(b% has_phase) allocate(b% phase(b% ntrack))
    b% eep = a(m)% eep(1:b% neep)
    if(a(m)% has_phase) b% phase = a(m)% phase(1:b% ntrack)
    b% tr = 0d0
    b% dist = 0d0

    dm = b% initial_mass - initial_mass(m)
    x = initial_mass(mlo:mhi)
    do i=1,b% ntrack
       do j=1,b% ncol
          forall(k=1:4) y(k) = a(mlo-1+k)% tr(j,i)
          call interp_4pt_pm(x, y, f)
          b% tr(j,i) = y(2) + dm*(f(1) + dm*(f(2) + dm*f(3)))
       enddo
    enddo

    write(*,*) ' ierr = ', ierr

    deallocate(initial_mass)
  end subroutine interpolate_track

  subroutine read_iso_input(ierr)
    use iso_color, only: iso_color_init
    integer, intent(out) :: ierr
    character(len=col_width) :: col_name
    character(len=10) :: list_type, age_type
    character(len=6) :: eep_style
    integer :: i
    real(dp) :: age_low, age_high, age_step
    ierr=0
    ntrk=0
    io=alloc_iounit(ierr)
    open(io, file='input.nml', action='read', status='old', iostat=ierr)
    read(io, nml=iso_controls, iostat=ierr)
    close(io)

    if(ierr/=0) then
       write(0,'(a)') 'make_iso: problem reading iso_controls namelist'
       return
    endif

    if(do_colors)then
       set% Av = extinction_Av
       set% Rv = extinction_Rv
       set% cmd_suffix = cmd_suffix
       call iso_color_init(BC_table_list,do_Cstars,Cstar_table_list,ierr)
       if(ierr/=0) write(0,'(a)') ' problem reading BC_table_list = ', trim(BC_table_list)
    endif

    call get_command_argument(1,input_file)

    !read info about into tracks
    open(unit=io,file=trim(input_file),status='old',action='read',iostat=ierr)
    if(ierr/=0)then
       write(0,*) ' make_iso: problem reading ', trim(input_file)
       return
    endif
    read(io,*) !skip comment line
    read(io,'(a)') history_dir
    read(io,'(a)') eep_dir
    read(io,'(a)') iso_dir
    read(io,*) !skip comment line
    read(io,'(a)') history_columns_list
    read(io,*) !skip comment line
    read(io,*) ntrk
    allocate(s(ntrk))
    do i=1,ntrk
       read(io,'(a)',iostat=ierr) s(i)% filename 
       s(i)% filename = trim(s(i)% filename)
       if(ierr/=0) exit
    enddo
    !read info about output isochrones
    read(io,*) !skip this line
    read(io,*) set% filename
    set% filename = trim(iso_dir) // '/' // trim(set% filename)
    read(io,'(a)') list_type
    read(io,'(a)') age_type
    if(trim(age_type)=='linear')then
       age_scale = age_scale_linear
    elseif(trim(age_type)=='log10')then
       age_scale = age_scale_log10
    else
       stop ' make_iso: age scale must be given as "linear" or "log10"'
    endif
    read(io,*) niso
    allocate(set% iso(niso))
    if(trim(list_type)=='min_max') then
       read(io,*) age_low
       read(io,*) age_high
       if(age_high < age_low) stop '  make_iso: max age < min age'
       !assign ages
       if(niso > 1) then
          age_step = (age_high - age_low)/real(niso-1,kind=dp)
       else
          age_step = 0d0
       endif
       do i=1,niso
          set% iso(i)% age = age_low + age_step*real(i-1,kind=dp)
       enddo
    else if(trim(list_type)=='list') then
       do i=1,niso
          read(io,*) set% iso(i)% age
       enddo
    else
       stop ' make_iso: ages must be given as "list" or "min_max"'
    endif

    read(io,'(a6)',iostat=ierr) eep_style
    if(ierr==0 .and. eep_style=='double') then 
       use_double_eep=.true.
    else
       use_double_eep=.false.
    endif
    close(io)
    call free_iounit(io)

    ! this section reads the column names to use from history_columns.list
    ! and locates those that need to be identified for isochrone construction
    call process_history_columns(history_columns_list,cols,ierr)
    ncol = size(cols)

    col_name = 'star_age'; i_age = locate_column(col_name,cols)
    col_name = 'star_mass'; i_mass= locate_column(col_name,cols)

    ! hack: replace star_age column with initial_mass in isochrones
    i_Minit = i_age
    cols(i_Minit)   = 'initial_mass'

  end subroutine read_iso_input

  logical function monotonic(array)
    real(dp), intent(in) :: array(:)
    logical :: ascending = .false.
    integer :: i, n
    monotonic = .false.
    n = size(array)
    if(n<=2)then
       write(*,*) ' Warning, monotonic: array of length <= 2'
    else
       ascending = array(1) <= array(n)
       if(ascending)then
          do i=1,n-1
             if(array(i) > array(i+1)) then 
                if(iso_debug) write(*,*) array(i), ' > ', array(i+1)
                return
             endif
          enddo
       else !descending
          do i=1,n-1
             if(array(i) < array(i+1)) then
                if(iso_debug) write(*,*) array(i), ' > ', array(i+1)
                return
             endif
          enddo
       endif
    endif
    monotonic = .true.
  end function monotonic

  subroutine check_monotonic(array,ierr)
    real(dp), intent(in) :: array(:)
    integer, intent(out) :: ierr
    integer :: i
    ierr=0
    if(.not.monotonic(array))then
       ierr=-1
       write(*,*) ' array not monotonic '
       do i=1,size(array)
          write(*,*) i, array(i)
       enddo
    endif
  end subroutine check_monotonic

  subroutine smooth(x,y)
    real(dp), intent(in) :: x(:)
    real(dp), intent(inout) :: y(:)
    integer :: i,n
    n=size(y)
    if(n<8) return
    y(2)=npoint(x(1:3),y(1:3),x(2))
    y(n-1)=npoint(x(n-2:n),y(n-2:n),x(n-1))
    do i=3,n-2
       y(i)=npoint(x(i-2:i+2),y(i-2:i+2),x(i))
    enddo
  end subroutine smooth

  function npoint(x,y,x0) result(y0)
    real(dp), intent(in)  :: x(:), y(:), x0
    real(dp) :: y0
    real(dp), allocatable :: w(:)
    real(dp), parameter :: eps=1d-1
    allocate(w(size(x)))
    w = 1d0/(eps+abs(x-x0))
    y0 = sum(y*w)/sum(w)      
    deallocate(w)
  end function npoint
  
  subroutine PAV(y) !pool-adjancent violators algorithm applied in place
    !based on python implementation at https://gist.github.com/fabianp/3081831
    real(dp), intent(inout) :: y(:)
    integer :: i, n, start, last, m
    real(dp), allocatable :: d(:)
    integer, allocatable :: lvls(:,:)
    n=size(y)
    allocate(d(n-1),lvls(n,2))
    do i=1,n
       lvls(i,:)=i
    enddo
    do while(.true.)
       d=y(2:n)-y(1:n-1)
       if(all(d>=0)) exit !test for monotonicity
       i=locate(d) !finds the first point in d that is < 0
       start = lvls(i,1)
       last = lvls(i+1,2)
       m = last - start + 1
       y(start:last) = sum(y(start:last))/real(m)
       lvls(start:last,1)=start
       lvls(start:last,2)=last
    enddo
  end subroutine PAV  

  integer function locate(y)
    real(dp), intent(in) :: y(:)
    integer :: i, n
    n=size(y)
    do i=1,n
       if(y(i)<0.0)then
          locate=i
          return
       endif
    enddo
    locate=0
  end function locate

end program make_isochrone
