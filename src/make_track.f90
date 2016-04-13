program make_track

  !MESA modules
  use interp_1d_def
  use interp_1d_lib

  !local modules
  use iso_eep_support

  implicit none

  character(len=file_path) :: input_file
  type(track), allocatable :: s(:), t(:) !existing set
  integer :: i, ierr=0, num_tracks_t=0, num_tracks_s=0
  logical, parameter :: debug=.false.
  logical :: output_to_eep_dir = .false.
  namelist /track_controls/ output_to_eep_dir

  !check command line arguments
  if(command_argument_count()<1) then
     write(*,*) '   make_track                 '
     write(*,*) '   usage: ./make_track <input>'
     stop       '   no command line argument   '
  else
     call get_command_argument(1,input_file)
  endif

  !read input file
  call read_input
  if(ierr/=0) stop 'make_track: failed in read_input'

  !read in existing tracks
  do i=1,num_tracks_s
     call read_eep(s(i))
     if(debug) write(*,'(a50,f8.2,99i8)') trim(s(i)% filename), s(i)% initial_mass, s(i)% eep
  enddo

  do i=1,num_tracks_t
     call interpolate_track(s,t(i),ierr)
     if(ierr/=0) then
        write(0,*) 'make_track: interpolation failed for ', trim(t(i)% filename)
        write(0,*) '            no output written '
        cycle
     endif
     if(output_to_eep_dir) t(i)% filename = trim(eep_dir) // '/' // trim(t(i)% filename)
     call write_track(t(i))
  enddo

contains

  !takes a set of EEP-defined tracks and interpolates a new
  !track for the desired initial mass
  subroutine interpolate_track(a,b,ierr)
    type(track), intent(in) :: a(:)
    type(track), intent(inout) :: b
    integer, intent(out) :: ierr
    real(dp) :: f(3), dx, x(4), y(4)
    real(dp), pointer :: initial_mass(:)=>NULL() !(n)
    integer :: i, j, k, m, mlo, mhi, n

    n = size(a)

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
  end subroutine interpolate_track

  subroutine read_input
    integer :: io
    character(len=file_path) :: eep_file

    io=alloc_iounit(ierr)

    open(unit=io,file='input.nml', action='read', status='old', iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_track: problem reading input.nml '
       return
    endif
    read(io, nml=track_controls, iostat=ierr)
    close(io)

    open(unit=io,file=trim(input_file),status='old',action='read',iostat=ierr)
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
       read(io,*) t(i)% initial_mass, t(i)% filename
    enddo
    close(io)

    open(unit=io,file=trim(eep_file),status='old',action='read',iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_track: problem reading ', trim(eep_file)
       return
    endif
    read(io,*) !skip first line
    read(io,*) !skip
    read(io,'(a)') eep_dir
    read(io,*) !skip
    read(io,*) !skip comment line
    read(io,*) !skip
    read(io,*) !skip comment line
    read(io,*) num_tracks_s
    allocate(s(num_tracks_s))
    do i=1,num_tracks_s
       read(io,'(a)',iostat=ierr) s(i)% filename
       if(ierr/=0) exit
    enddo
    close(io)

    call free_iounit(io)

  end subroutine read_input

end program make_track
