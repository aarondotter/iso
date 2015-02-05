program make_isochrone

  !MESA modules
  use const_def, only: dp
  use utils_lib
  use interp_1d_def
  use interp_1d_lib

  !local modules
  use iso_eep_support
  use iso_color

  implicit none

  character(len=file_path) :: input_file, iso_file, history_columns_list
  integer :: i, ierr, io, j, niso, ntrk, ngood, first, prev, i_Minit
  type(track), allocatable :: s(:), t(:), q
  type(isochrone_set) :: set

  logical :: use_double_eep, do_colors
  integer, parameter :: piecewise_monotonic = 4
  logical, parameter :: iso_debug = .false.
  logical, parameter :: do_tracks = .false.
  logical, parameter :: do_isochrones = .true.
  logical, parameter :: do_smooth = .true.
  logical, parameter :: skip_non_monotonic = .false.

  ierr=0
  if(command_argument_count()<1) then
     write(*,*) '   make_iso                   '
     write(*,*) '   usage: ./make_iso <input>  '
     stop       '   no command line argument   '
  endif

  !begin
  call read_iso_input(ierr)
  if(ierr/=0) write(*,*) '  read_iso_input: ierr = ', ierr

  call read_color_input(ierr)
  do_colors = (ierr == 0)

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
     q% initial_mass = 1.2d0
     q% filename = trim(history_dir) // '/out.trk'
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
     call write_isochrones_to_file(iso_file,set)
     if(do_colors) call write_cmds_to_file(iso_file,set)
  endif

  !all done.
  deallocate(t,cols)

contains

  subroutine do_isochrone_for_age(s,iso)
    type(track), intent(in) :: s(:)
    type(isochrone), intent(inout) :: iso
    integer :: count, eep, hi, ierr, index, j, k, lo, max_eep, n
    integer :: interp_method, jlo, jhi, pass
    integer, parameter :: jinc = 6
    character(len=col_width) :: mass_age_string = 'mass from age'
    real(dp) :: age, mass
    real(dp), pointer :: ages(:)=>NULL(), masses(:)=>NULL()
    real(dp), allocatable :: result1(:,:), result2(:,:)
    logical, allocatable :: skip(:)
    integer, allocatable :: valid(:)

    ierr = 0
    n = size(s) ! is the number of tracks
    allocate(skip(n))

    !set method and options for interpolation
    !for interp_m3: average, quartic, or super_bee
    !interp_method = average
    !for interp_pm: piecewise_monotonic           
    interp_method = piecewise_monotonic

    !initialize some quantities
    age = iso% age
    max_eep = maxval(s(:)% ntrack)
    mass=0d0

    !this is temporary storage for the isochrone data:
    !result1 stores the data for all EEPs, valid tells
    !which ones are good and will be returned via the iso 
    !derived type
    allocate(result1(ncol,max_eep),result2(ncol,max_eep),valid(max_eep))
    result1 = 0d0
    result2 = 0d0
    valid = 0


    eep_loop: do eep=1,max_eep

       !determine tracks for which the ith eep is defined
       !the skip logical array determines whether or not a given
       !track will be included in the ensuing interpolation steps.
       !count keeps track of how many tracks will be used. if
       !fewer than 2 tracks satisfy the condition, skip the EEP
       skip = .false.
       count = 0
       hi = max_eep
       lo = eep
       do k=1,n
          if(s(k)% eep(1) > lo .or. s(k)% eep(s(k)% neep) < lo ) then
             skip(k) = .true.
          else
             count = count + 1
          endif
       enddo
       if(iso_debug) write(*,*) '  eep, count, n = ', eep, count, n
       if(count < 2)then
          if(iso_debug) write(*,*) 'not enough eeps to interpolate'
          cycle
       endif

       !allocate ages and masses that pass the above test
       !I use pointers here because they can be allocated
       !and deallocated to the proper size at each EEP
       if(associated(ages)) then
          if(iso_debug) write(0,*) "ages was associated for eep ", eep
          deallocate(ages)
          nullify(ages)
       endif
       if(associated(masses)) then
          if(iso_debug) write(0,*) "masses was associated for eep ", eep
          deallocate(masses)
          nullify(masses)
       endif
       allocate(ages(count),masses(count))

       !this step fills the masses and ages arrays that are
       !used as the basis for interpolation
       j=0
       do k=1,n
          if(.not.skip(k))then
             j=j+1
             ages(j) = log10(s(k)% tr(i_age,eep))
             masses(j) = s(k)% initial_mass
          endif
       enddo

       !do k=1,count
       !   write(112,*) eep, ages(k), log10(masses(k))
       !enddo


       if(do_smooth.and.all(masses>0.5,dim=1)) call smooth(masses,ages)

       !do k=1,count
       !   write(111,*)  eep, ages(k), log10(masses(k))
       !enddo

       !check to see if the input age is found within the
       !current set of ages. if not, skip to the next EEP.
       j = binary_search( count, ages, 1, age)
       if( j < 1 .or. j > count-1 ) then
          deallocate(masses,ages)
          cycle
       endif
       jlo = max(1,j-jinc)
       jhi = min(count,j+jinc)   

       if(iso_debug) write(*,'(i4,99f7.2)') eep, masses

       !check to see if masses and ages are monotonic
       !if not, then interpolation will fail
       if(.not.monotonic(masses)) then
          write(*,*) ' masses not monotonic in do_isochrone_for_age: ', age
          stop 99
       endif

       if(iso_debug .and. .not.monotonic(ages(jlo:jhi))) then
          write(*,'(a8,i5)') '  eep = ', eep
          do k=jlo,jhi
             write(*,'(a20,i5,f9.3,f12.8)') 'i, mass, age = ', k, masses(k), ages(k)
          enddo
          write(*,*) ' ages not monotonic in do_isochrone_for_age: ', age
       endif

       if(iso_debug) then 
          write(*,*) ' j, count = ', j, count
          write(*,*) skip
          do k=1,count
             write(*,*) masses(k), ages(k), age
          enddo
       endif

       !this block checks between two tracks at the current EEP to see if the 
       !age lies in between the two. if it does, then it outputs a linear mass 
       !interpolation. it is checking to see if this happens more than once.
       !j gives the number of times that this happens for a given EEP.
       if(use_double_eep)then
          pass=0 
          k_loop: do k=1,count-1

             if(.not.skip(k))then

                if( ages(k) > ages(k+1) .and. ages(k) >= age .and. ages(k+1) < age )then
                   pass = pass + 1
                   if(pass==1)then
                      call monotonic_mass_range(ages,k,jlo,jhi)
                      if(iso_debug) write(*,*) mass, masses(jlo), masses(jhi)

                      mass=interp_x_from_y(ages(jlo:jhi),masses(jlo:jhi),age,mass_age_string,ierr)
                      if(ierr/=0)then
                         write(0,*) ' age interpolation failed for eep, pass = ', eep, pass
                         cycle eep_loop
                      endif

                      write(*,'(2i5,f9.5,3f14.9)') eep, pass, age, mass, masses(jlo), masses(jhi)

                      result1(i_Minit,eep) = mass
                      do index = 2, ncol
                         result1(index,eep) = iso_interpolate(eep, interp_method, n, &
                              s, skip, count, index, masses, mass, cols(index), ierr)
                         if(ierr/=0) then
                            write(0,*) ' mass interpolation failed for index = ', &
                                 trim(cols(index))
                            cycle eep_loop
                         endif
                      enddo
                      valid(eep)=valid(eep)+1

                   else if(pass==2)then
                      call monotonic_mass_range(ages,k,jlo,jhi)
                      if(iso_debug) write(*,*) mass, masses(jlo), masses(jhi)

                      mass=interp_x_from_y(ages(jlo:jhi),masses(jlo:jhi),age,mass_age_string,ierr)
                      if(ierr/=0)then
                         write(0,*) ' age interpolation failed for eep, pass = ', eep, pass
                         cycle eep_loop
                      endif

                      write(*,'(2i5,f9.5,3f14.9)') eep, pass, age, mass, masses(jlo), masses(jhi)

                      result2(i_Minit,eep) = mass
                      do index = 2, ncol
                         result2(index,eep) = iso_interpolate(eep, interp_method, n, &
                              s, skip, count, index, masses, mass, cols(index), ierr)
                         if(ierr/=0) then
                            write(0,*) ' mass interpolation failed for index = ', &
                                 trim(cols(index))
                            cycle eep_loop
                         endif
                      enddo
                      valid(eep)=valid(eep)+1
                   endif

                endif

             endif

          enddo k_loop


       else ! single EEP case; this is the original method.
          ! interpolate in age to find the EEP's initial mass
          index = i_Minit     ! special case for iso_intepolate

          mass = iso_interpolate( eep, interp_method, n, &
               s, skip, count, index, ages, age, mass_age_string, ierr)
          if(ierr/=0)then
             write(0,*) ' interpolation failed in age->mass'
             cycle eep_loop
          endif
          result1(i_Minit,eep) = mass

          do index = 2, ncol
             result1(index,eep) = iso_interpolate(eep, interp_method, n, &
                  s, skip, count, index, masses, mass, cols(index), ierr)
             if(ierr/=0) then
                write(0,*) ' mass interpolation failed for index = ', trim(cols(index))
                cycle eep_loop
             endif
          enddo

          !completed all interpolations, so this is a valid EEP
          !if validi > 0 then it is a valid EEP
          valid(eep) = valid(eep)+1

       endif

       ! clean up for each iteration
       deallocate(ages,masses)
       nullify(ages,masses)
    enddo eep_loop

    if(skip_non_monotonic .and. .not.use_double_eep)then !check for non-monotonic EEPs
       do eep=2,max_eep
          if(valid(eep)>0)then
             do j=eep-1,2,-1
                if(valid(j)>0) exit
             enddo
             if( result1(i_Minit,eep) < result1(i_Minit,j) ) valid(eep)=0
          endif
       enddo
    endif

    !now result1 and valid are full for all EEPs,
    !we can pass the data to the iso derived type
    iso% ncol = ncol
    iso% neep = sum(valid)
    allocate(iso% cols(iso% ncol))
    !do eep=1,max_eep
    !   iso% neep = iso% neep + valid(eep)
    !enddo
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

  subroutine monotonic_mass_range(ages,k,jlo,jhi)
    !this subroutine determines the upper and lower limits that
    !should be used in the EEP interpolation for the case of 
    !double EEPs; default result is jlo=1, jhi=size(ages)
    real(dp), intent(in) :: ages(:)
    integer, intent(in) :: k
    integer, intent(out) :: jlo, jhi
    integer :: j,n
    n=size(ages)
    !find lower limit:
    jlo=1
    lower_loop: do j=k,1,-1
       if(ages(j+1) < ages(j)) then
          jlo=j
       else
          exit lower_loop
       endif
    enddo lower_loop
    !find upper limit
    jhi=n
    upper_loop: do j=jlo+1,n
       if(ages(j-1) > ages(j)) then
          jhi = j
       else
          exit upper_loop
       endif
    enddo upper_loop
  end subroutine monotonic_mass_range


  real(dp) function interp_x_from_y(x,y,x_in,label,ierr)
    !     subroutine interpolate_vector_pm( &
    !     n_old, x_old, n_new, x_new, v_old, v_new, work1, str, ierr)
    !     use interp_1d_def, only: pm_work_size

    !     real(dp), intent(in) :: x_old(:) !(n_old)
    !     real(dp), intent(in) :: v_old(:) !(n_old)
    !     real(dp), intent(in) :: x_new(:) !(n_new)
    !     real(dp), intent(out) :: v_new(:) ! (n_new)
    !     real(dp), intent(inout), pointer :: work1(:) ! =(n_old, pm_work_size)
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
    real(dp) :: y     !f(4,count), work(count,nwork)
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
       write(*,*) ' eep = ', eep
       do k=1,n
          write(*,*) x_array(k), f(1,k), skip(k)
       enddo
       write(*,*) x
       write(*,*) y
       write(*,*) 
       write(*,*)
       stop
    endif

    iso_interpolate = y

  end function iso_interpolate


  !writes a series of n isochrones to filename
  subroutine write_isochrones_to_file(filename,set)
    character(len=file_path), intent(in) :: filename
    type(isochrone_set), intent(in) :: set
    integer :: i, ierr, io,n
    io=alloc_iounit(ierr)
    n=size(set% iso)
    write(0,*) ' isochrone output file = ', trim(filename)
    open(io,file=trim(filename),action='write',iostat=ierr)
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
    write(io,'(a5,299a32)') '# EEP', 'log_age', adjustr(iso% cols)
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
    write(io,'(a5,299a32)') '# EEP', 'log_age', adjustr(iso% cols), 'phase'
    do i=1,iso% neep
       write(io,'(i5,299(1pes32.16e3))') iso% eep(i), iso% age, &
            iso% data(:,i), real(iso% phase(i),kind=dp)
    enddo
  end subroutine write_isochrone_to_file_phase

  subroutine write_cmds_to_file(filename,s)
    character(len=256), intent(in) :: filename
    type(isochrone_set), intent(inout) :: s
    character(len=256) :: output
    integer :: i, io, n, ierr
    if(color_suffix/='') output = trim(filename) // '.' // trim(color_suffix)
    n=size(s% iso)
    write(0,*) ' cmd output file = ', trim(output)
    io = alloc_iounit(ierr)
    if(ierr/=0) return
    open(io,file=trim(output),action='write',iostat=ierr)
    if(ierr/=0) return
    write(io,'(a25,i5)')    '# number of isochrones = ', n
    write(io,'(a25,i5)')    '# MESA version number  = ', s% version_number
    write(io,'(a25,2f6.3)') '# values of Av and Rv  = ', bc% Av(1), bc% Rv(1)
    do i=1,n
       call write_cmd_to_file(io,s% iso(i))
       if(i<n) write(io,*)
       if(i<n) write(io,*)
    enddo
    close(io)
    call free_iounit(io)
  end subroutine write_cmds_to_file

  subroutine write_cmd_to_file(io,iso)
    integer, intent(in) :: io
    type(isochrone), intent(inout) :: iso
    integer :: i, iT, ig, iL
    iT=0; ig=0; iL=0

    do i=1, iso% ncol
       if(trim(iso% cols(i)) == 'log_Teff') then
          iT=i
       else if(trim(iso% cols(i)) == 'log_g') then
          ig=i
       else if(trim(iso% cols(i))== 'log_L') then
          iL=i
       endif
    enddo
    call get_mags(iso,iT,ig,iL)

    write(io,'(a25,2i5)') '# number of EEPs, cols = ', iso% neep, iso% nfil + 5
    write(io,'(a1,i4,4i32,299i12)') '#    ', (i,i=1,iso% nfil+5)
    write(io,'(a5,4a32,299a12)') '# EEP', 'log_age', 'log_Teff', 'log_g', 'log_L', &
         adjustr(iso% labels)
    do i = 1,iso% neep
       write(io,'(i5,4(1pes32.16e3),299(0pf12.6))') iso% eep(i), iso% age,  &
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
    integer, intent(out) :: ierr
    character(len=col_width) :: col_name
    character(len=10) :: list_type
    character(len=6) :: eep_style
    integer :: i
    real(dp) :: age_low, age_high, age_step
    ierr=0
    call get_command_argument(1,input_file)
    io=alloc_iounit(ierr)
    !read info about into tracks
    open(unit=io,file=trim(input_file))
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
    read(io,*) iso_file
    iso_file = trim(iso_dir) // '/' // trim(iso_file)
    read(io,'(a)') list_type
    read(io,*) niso
    allocate(set% iso(niso))
    if(trim(list_type)=='min_max') then
       read(io,*) age_low
       read(io,*) age_high
       if(age_high < age_low) stop '  make_iso: max age < min age'
       !assign ages
       if(niso > 1) then
          age_step = (age_high - age_low)/dble(niso-1)
       else
          age_step = 0d0
       endif
       do i=1,niso
          set% iso(i)% age = age_low + age_step*dble(i-1)
       enddo
    else if(trim(list_type)=='list') then
       do i=1,niso
          read(io,*) set% iso(i)% age
       enddo
    else
       stop ' make_iso: ages must be given as "list" or "min_max"'
    endif

    read(io,'(a6)') eep_style
    if(eep_style=='double') then 
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

end program make_isochrone
