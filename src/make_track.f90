program make_track

  !MESA modules
  use interp_1d_def
  use interp_1d_lib

  !local modules
  use iso_eep_support
  use iso_eep_color
  use interp_support

  implicit none

  character(len=file_path) :: input_file
  character(len=32) :: phot_string, arg
  type(track), allocatable :: s(:), t(:) !existing set
  integer :: i, ierr=0, num_tracks_t=0, num_tracks_s=0
  logical, parameter :: debug=.false.
  logical :: output_to_eep_dir = .false.
  logical :: set_fixed_Fe_div_H = .false., do_CMDs = .false.
  character(len=file_path) :: BC_table_list = '', cmd_suffix = 'cmd'
  real(sp) :: extinction_Av=0.0, extinction_Rv=0.0, Fe_div_H = 0.0

  namelist /cmd_controls/ BC_table_list, extinction_Av, extinction_Rv, &
  cmd_suffix, Fe_div_H, set_fixed_Fe_div_H

  namelist /track_controls/ output_to_eep_dir

  force_linear=.true.

  !check command line arguments
  if(command_argument_count()<1) then
     write(*,*) '   make_track                 '
     write(*,*) '   usage: ./make_track [input] [phot_string] [Av]'
     write(*,*) '   input = input file         '
     write(*,*) '   (optional) phot_string = UBVRIJHKsKp, etc.'
     write(*,*) '   (optional) Av = extinction, 0 <= Av <= 6, defaults to zero'
     stop       '   no command line argument   '
  endif

  call get_command_argument(1,input_file)

  !read input file
  call read_input
  if(ierr/=0) stop 'make_track: failed in read_input'

  if(command_argument_count()>1) then
     call get_command_argument(2,phot_string)
     call color_init(phot_string, BC_table_list, set_fixed_Fe_div_H, Fe_div_H, ierr)
     do_CMDs = ierr==0
  endif

  if(do_CMDs .and. command_argument_count()>2) then
     call get_command_argument(3,arg)
     read(arg,*) extinction_Av
  else
     extinction_Av = 0.0
  endif

  !read in existing tracks
  do i=1,num_tracks_s
     call read_eep(s(i))
     if(debug) write(*,'(a50,f8.2,99i8)') trim(s(i)% filename), s(i)% initial_mass, s(i)% eep
  enddo

  do i=1,num_tracks_t
     call interpolate_mass(s,t(i),ierr)
     if(ierr/=0) then
        write(0,*) 'make_track: interpolation failed for ', trim(t(i)% filename)
        write(0,*) '            no output written '
        cycle
     endif
     if(output_to_eep_dir) t(i)% filename = trim(eep_dir) // '/' // trim(t(i)% filename)
     call write_track(t(i))
     if(do_CMDs) then
        t(i)% Av = extinction_Av
        t(i)% cmd_suffix = cmd_suffix
        call write_track_cmd_to_file(t(i))
     endif
  enddo

contains

  subroutine read_input
    integer :: io, i, j, k
    character(len=file_path) :: eep_file, data_line

    open(newunit=io,file='input.nml', action='read', status='old', iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_track: problem reading input.nml '
       return
    endif
    read(io, nml=track_controls, iostat=ierr)
    rewind(io)
    read(io, nml=cmd_controls, iostat=ierr)
    close(io)

    open(newunit=io,file=trim(input_file),status='old',action='read',iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_track: problem reading ', trim(input_file)
       return
    endif
    read(io,*) !skip comment
    read(io,'(a)') eep_file
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

    open(newunit=io,file=trim(eep_file),status='old',action='read',iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_track: problem reading ', trim(eep_file)
       return
    endif
    read(io,*) !skip
    read(io,*) !skip
    read(io,*) !skip
    read(io,*) !skip
    read(io,*) !skip
    read(io,*) !skip
    read(io,'(a)') eep_dir
    read(io,*) !skip
    read(io,*) !skip
    read(io,*) !skip
    read(io,*) !skip
    read(io,*) num_tracks_s

    allocate(s(num_tracks_s))
    do i=1,num_tracks_s
       read(io,'(a)',iostat=ierr) s(i)% filename
       if(ierr/=0) exit
    enddo
    close(io)

  end subroutine read_input

end program make_track
