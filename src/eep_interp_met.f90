program eep_interp_feh

  !MESA modules
  use interp_1d_def
  use interp_1d_lib

  !local modules
  use iso_eep_support
  use iso_eep_color
  use interp_support

  implicit none

  type(track), allocatable :: s1(:), s2(:), t(:)
  character(len=file_path), allocatable :: input_grid(:)
  character(len=file_path) :: input_file, phot_string
  real(dp), allocatable :: Fe_div_H(:)
  real(dp) :: new_Fe_div_H
  integer :: i, ierr, num_grids, num_tracks_s1, num_tracks_s2, num_tracks_t, lo, hi
  logical, parameter :: debug=.false.
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

  !set this to use linear interpolation in mass
  force_linear=.true.
  
  call process_command_args

  call read_grid_lists
  if(ierr/=0) stop 'make_track: failed in read_grid_lists'

  call read_input
  if(ierr/=0) stop 'make_track: failed in read_input'


  if(do_CMDs) then
     call color_init(phot_string, BC_table_list, &
          do_Cstars, Cstar_table_list, &
          set_fixed_Fe_div_H, fixed_Fe_div_H, ierr)
     if(ierr/=0) then
        write(0,*) 'failed to initialize BC tables'
        do_CMDs = .false.
     endif
  endif

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
          do_CMDs = .true.
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
