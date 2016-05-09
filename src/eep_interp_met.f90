program eep_interp_feh

  !MESA modules
  use interp_1d_def
  use interp_1d_lib

  !local modules
  use iso_eep_support
  use iso_eep_color
  use iso_interp_support

  implicit none

  type(track), allocatable :: s1(:), s2(:), t(:)
  character(len=file_path), allocatable :: input_grid(:)
  character(len=file_path) :: input_file, phot_string
  real(dp), allocatable :: Fe_div_H(:)
  real(dp) :: new_Fe_div_H
  integer :: i, ierr, num_grids, num_tracks_s1, num_tracks_s2, num_tracks_t, lo, hi
  logical, parameter :: debug=.true.
  logical :: output_to_eep_dir = .false., do_Cstars = .false., do_Z_interp
  logical :: set_fixed_Fe_div_H = .false., do_CMDs = .false.
  character(len=file_path) :: BC_table_list = '', cmd_suffix = 'cmd'
  character(len=file_path) :: Cstar_table_list = ''
  real(sp) :: extinction_Av=0.0, extinction_Rv=0.0, fixed_Fe_div_H=0.0

  namelist /cmd_controls/ BC_table_list, extinction_Av, extinction_Rv, &
       Cstar_table_list, do_Cstars, cmd_suffix, fixed_Fe_div_H, set_fixed_Fe_div_H

  namelist /track_controls/ output_to_eep_dir

  !set this for increased speed
  make_bin_eeps=.true.
  
  call process_command_args

  call read_grid_lists
  if(ierr/=0) stop 'make_track: failed in read_grid_lists'

  call read_input
  if(ierr/=0) stop 'make_track: failed in read_input'

  t(:)% Fe_div_H = new_Fe_div_H

  do i=1,num_tracks_t
     if(do_Z_interp)then
        call interpolate_mass_Z(s1,s2,t(i),ierr)
     else
        call interpolate_mass(s1,t(i),ierr)
     endif
     if(ierr/=0) then
        write(0,*) 'make_track: interpolation failed for ', trim(t(i)% filename)
        write(0,*) '            no output written '
        cycle
     endif

     if(output_to_eep_dir) t(i)% filename = trim(eep_dir) // '/' // trim(t(i)% filename)
     call write_track(t(i))
     if(do_CMDs)then
        t(i)% Av = extinction_Av
        t(i)% cmd_suffix = cmd_suffix
        call write_track_cmd_to_file(t(i))
     endif
  enddo

contains

  subroutine interpolate_mass_Z(a1,a2,b,ierr)
    type(track), intent(in) :: a1(:), a2(:)
    type(track), intent(inout) :: b
    integer, intent(out) :: ierr
    type(track) :: t1, t2
    integer :: i
    real(dp) :: alfa, beta
    ierr=0

    t1% initial_mass = b% initial_mass
    t2% initial_mass = b% initial_mass

    call interpolate_mass(a1,t1,ierr)
    if(ierr/=0) return
    call interpolate_mass(a2,t2,ierr)
    if(ierr/=0) return

    alfa = (t2% Fe_div_H - b% Fe_div_H)/(t2% Fe_div_H - t1% Fe_div_H)
    beta = 1d0 - alfa
    write(*,*) ' alfa, beta = ', alfa, beta

    b% star_type = t1% star_type
    b% alpha_div_Fe = t1% alpha_div_Fe
    b% v_div_vcrit = t1% v_div_vcrit
    b% version_string = t1% version_string
    b% MESA_revision_number = t1% MESA_revision_number
    b% initial_mass = t1% initial_mass
    b% has_phase = t1% has_phase
    b% ncol = t1% ncol
    allocate(b% cols(b% ncol))
    b% cols = t1% cols
    call set_initial_Y_and_Z_for_eep(b)

    b% neep = min(t1% neep, t2% neep)
    allocate(b% eep(b% neep))
    b% eep = t1% eep(1:b% neep)
    b% ntrack = min(t1% ntrack, t2% ntrack)
    allocate(b% tr(b% ncol, b% ntrack))

    write(*,*) ' t1% ncol = ', t1% ncol
    write(*,*) ' t2% ncol = ', t2% ncol
    write(*,*) ' b% ncol = ', b% ncol

    do i=1,b% ntrack
       b% tr(:,i) = alfa*t1% tr(:,i) + beta*t2% tr(:,i)      
    enddo
    if(b% has_phase)then
       allocate(b% phase(b% ntrack))
       b% phase = t1% phase(1:b% ntrack)
    endif
  end subroutine interpolate_mass_Z

  !takes a set of EEP-defined tracks and interpolates a new
  !track for the desired initial mass
  subroutine interpolate_mass(a,b,ierr)
    type(track), intent(in) :: a(:)
    type(track), intent(inout) :: b
    integer, intent(out) :: ierr
    real(dp) :: f(3), dx, x(4), y(4)
    real(dp), pointer :: initial_mass(:)=>NULL() !(n)
    integer :: i, j, k, m, mlo, mhi, n

    n = size(a)

    nullify(initial_mass)
    allocate(initial_mass(n))
    initial_mass = a(:)% initial_mass
    m = binary_search(n, initial_mass, 1, b% initial_mass)
    if(debug)then
       write(*,*) b% initial_mass
       write(*,*) initial_mass(m:m+1)
    endif

    mlo = min(max(1,m-1),n-3)
    mhi = max(min(m+2,n),4)

    if(debug)then
       write(*,*) '   mlo, m, mhi = ', mlo, m, mhi
       write(*,*) initial_mass(mlo:mhi)
       write(*,*) a(mlo:mhi)% neep
    endif

    k = minloc(a(mlo:mhi)% neep,dim=1) + mlo - 1

    b% neep = a(k)% neep
    b% ntrack = a(k)% ntrack
    b% star_type = a(m)% star_type
    b% version_string = a(m)% version_string
    b% initial_Z = a(m)% initial_Z
    b% initial_Y = a(m)% initial_Y
    b% Fe_div_H = a(m)% Fe_div_H
    b% alpha_div_Fe = a(m)% alpha_div_Fe
    b% v_div_vcrit = a(m)% v_div_vcrit
    b% ncol = a(m)% ncol
    b% has_phase = a(m)% has_phase
    allocate(b% cols(b% ncol))
    b% cols = a(m)% cols
    b% MESA_revision_number = a(m)% MESA_revision_number
    allocate(b% eep(b% neep))
    allocate(b% tr(b% ncol, b% ntrack))
    if(b% has_phase) allocate(b% phase(b% ntrack))
    b% eep = a(m)% eep(1:b% neep)
    if(a(m)% has_phase) b% phase = a(m)% phase
    b% tr = 0d0

    x = initial_mass(mlo:mhi)
    dx = b% initial_mass - x(2)

    do i=1,b% ntrack
       do j=1,b% ncol
          do k=1,4
             y(k) = a(mlo-1+k)% tr(j,i)
          enddo
          call interp_4pt_pm(x, y, f)
          b% tr(j,i) = y(2) + dx*(f(1) + dx*(f(2) + dx*f(3)))
       enddo
    enddo

    deallocate(initial_mass)
    if(debug) write(*,*) ' ierr = ', ierr
  end subroutine interpolate_mass


  subroutine read_input
    integer :: io, i, j, k
    character(len=file_path) :: data_line

    io=alloc_iounit(ierr)

    open(unit=io,file='input.nml', action='read', status='old', iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_track: problem reading input.nml '
       return
    endif
    read(io, nml=track_controls, iostat=ierr)
    rewind(io) !i've always wanted to use rewind!
    read(io, nml=cmd_controls, iostat=ierr)
    close(io)

    open(unit=io,file=trim(input_file),status='old',action='read',iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_track: problem reading ', trim(input_file)
       return
    endif

    read(io,*) !skip comment
    read(io,*) new_Fe_div_H
    read(io,*) !skip comment
    read(io,*) num_tracks_t

    allocate(t(num_tracks_t))

    read(io,*) !skip comment
    do i=1,num_tracks_t
       read(io,'(a)') data_line
       j=index(data_line, ' ')
       k=len_trim(data_line)
       read(data_line(1:j-1),*) t(i)% initial_mass
       read(data_line(j+1:k),'(a)') t(i)% filename
    enddo
    close(io)
    call free_iounit(io)

    lo=0; hi=0
    do_Z_interp=.false.

    do i=1,num_grids
       if(abs(Fe_div_H(i) - new_Fe_div_H) < 1d-6)then
          lo = i
          hi = lo
          exit
       endif
    enddo

    if(lo==0)then
       do i=1,num_grids-1
          if(Fe_div_H(i) < new_Fe_div_H .and. Fe_div_H(i+1) > new_Fe_div_H)then
             lo = i
             hi = i+1
             exit
          endif
       enddo
    endif

    if(lo==0)then
       write(0,*) ' [Fe/H] not within the grid! '
       return
    endif

    do_Z_interp = lo < hi

    if(debug)then
       write(0,*) ' new [Fe/H] = ', new_Fe_div_H
       write(0,*) ' grid [Fe/H] = ', Fe_div_H(lo)
       if(do_Z_interp) write(0,*) ' grid [Fe/H] = ', Fe_div_H(hi)
    endif

    call read_eep_file(s1,lo,num_tracks_s1)
    if(do_Z_interp) call read_eep_file(s2,hi,num_tracks_s2)

    if(debug)then
       write(0,*) 'num_tracks_s1 = ', num_tracks_s1
       if(do_Z_interp) write(0,*) 'num_tracks_s2 = ', num_tracks_s2
    endif

  end subroutine read_input

  subroutine read_eep_file(x,q,num)
    type(track), allocatable, intent(inout) :: x(:)
    integer, intent(in) :: q
    integer, intent(out) :: num
    integer :: i, io
    io=alloc_iounit(ierr)
    open(unit=io,file=trim(input_grid(q)),status='old',action='read',iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_track: problem reading ', trim(input_grid(q))
       return
    endif
    read(io,*) !skip1
    read(io,*) !skip2
    read(io,*) !skip3
    read(io,*) !skip4
    read(io,*) !skip5
    read(io,*) !skip6
    read(io,'(a)') eep_dir
    read(io,*) !skip8
    read(io,*) !skip9
    read(io,*) !skip10
    read(io,*) !skip11
    read(io,*) num

    allocate(x(num))
    do i=1,num
       read(io,'(a)',iostat=ierr) x(i)% filename
       if(ierr/=0) exit
       call read_eep(x(i))
    enddo
    close(io)

    call free_iounit(io)

  end subroutine read_eep_file


  subroutine read_grid_lists
    integer :: i, io
    real(dp) :: junk
    ierr=0
    io=alloc_iounit(ierr)
    open(io,file='input_grids.list',action='read',status='old',iostat=ierr)
    read(io,*) num_grids
    allocate(input_grid(num_grids),Fe_div_H(num_grids))
    if(debug) write(0,*) ' number of grids = ', num_grids
    do i=1,num_grids
       read(io,'(a)') input_grid(i)
       if(debug) write(0,*) ' grid file: ', trim(input_grid(i))
    enddo
    close(io)

    do i=1,num_grids
       open(io,file=trim(input_grid(i)),action='read',status='old',iostat=ierr)
       read(io,*)
       read(io,*)
       read(io,*)
       read(io,*) junk, junk, Fe_div_H(i), junk, junk
       close(io)
       if(debug) write(0,*) trim(input_grid(i)), Fe_div_H(i)
    enddo

    call free_iounit(io)

    if(debug) write(0,*) ' ierr = ', ierr
  end subroutine read_grid_lists

  subroutine process_command_args
    character(len=32) :: arg
    if(command_argument_count() < 1) then
       write(0,*) ' eep_interp_feh '
       write(0,*) ' usage: ./eep_interp_feh [input] [phot_string] [Av]'
       stop 'no command line argument'
    else
       call get_command_argument(1,input_file)
       if(command_argument_count()>1)then
          call get_command_argument(2,phot_string)
          call color_init(phot_string, BC_table_list, &
               do_Cstars, Cstar_table_list, &
               set_fixed_Fe_div_H, fixed_Fe_div_H, ierr)
          do_CMDs = ierr==0
       endif
    endif

    if(do_CMDs .and. command_argument_count()>2) then
       call get_command_argument(3,arg)
       read(arg,*) extinction_Av
    else
       extinction_Av = 0.0
    endif

  end subroutine process_command_args

end program eep_interp_feh
