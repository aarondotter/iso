program make_eeps

  !local modules
  use iso_eep_support
  use eep
  use phases

  implicit none

  integer :: i, ierr, io, num
  character(len=8) :: version_string='unknown'
  character(len=file_path) :: input_file, history_columns_list
  character(len=file_path), allocatable :: history_files(:)
  type(track), pointer :: t=>NULL(), s=>NULL()
  logical :: do_phases = .true.
  real(dp) :: initial_Y, initial_Z, Fe_div_H, v_div_vcrit, alpha_div_Fe

  namelist /eep_controls/ do_phases, center_gamma_limit, &
       center_carbon_limit, log_center_T_limit, &
       high_mass_limit, very_low_mass_limit, weight_center_rho_T_by_Xc, &
       Teff_scale, logL_scale, age_scale, Tc_scale, Rhoc_scale, &
       make_bin_tracks

  ierr=0

  if(command_argument_count()<1) then
     write(*,*) '   make_eeps                  '
     write(*,*) '   usage: ./make_eeps <input> '
     stop       '   no command line argument   '
  endif

  !read input file, set up columns, eeps, format specs
  call read_input(ierr)
  if(ierr/=0) stop 'make_eeps: failed in read_input'

  ! allocate tracks for history files, read them in and convert to eep
  do i=1,num
     call alloc_track(history_files(i),t)
     call read_history_file(t,ierr)
     write(*,*) trim(t% filename), t% neep, t% MESA_revision_number
     !now set header info
     t% initial_Y = initial_Y
     t% initial_Z = initial_Z
     t% Fe_div_H  = Fe_div_H
     t% alpha_div_Fe = alpha_div_Fe
     t% v_div_vcrit = v_div_vcrit
     t% version_string = version_string
     if(ierr/=0) then
        write(0,*) 'make_eep: problem reading!'
        cycle
     endif
     call primary_eep(t)
     write(*,'(99i8)') t% eep
     if( all(t% eep == 0) ) then
        write(*,*) ' PROBLEM WITH TRACK: NO EEPS DEFINED '
     else
        call alloc_track(t% filename,s)
        call secondary_eep(t,s)
        if(do_phases) call set_track_phase(s)
        s% filename = trim(eep_dir) // '/' // trim(s% filename) // '.eep'
        call write_track(s)
        deallocate(s)
        nullify(s)
     endif
     deallocate(t)
     nullify(t)
  enddo

contains

  subroutine read_input(ierr)
    integer, intent(out) :: ierr

    ierr=0
    io=alloc_iounit(ierr)

    open(unit=io,file='input.nml', action='read', status='old', iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_eeps: problem reading input.nml '
       return
    endif
    read(io, nml=eep_controls, iostat=ierr)
    close(io)

    version_string = adjustr(version_string)

    call get_command_argument(1,input_file)

    open(unit=io,file=trim(input_file),status='old',action='read',iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_eeps: problem reading ', trim(input_file)
       return
    endif
    read(io,*) !skip comment
    read(io,'(a8)') version_string
    read(io,*) !skip comment
    read(io,*) initial_Y, initial_Z, Fe_div_H, alpha_div_Fe, v_div_vcrit
    read(io,*) !skip comment
    read(io,'(a)') history_dir
    read(io,'(a)') eep_dir
    read(io,'(a)') iso_dir
    read(io,*) !skip comment
    read(io,'(a)') history_columns_list
    read(io,*) !skip comment
    read(io,*) num
    allocate(history_files(num))
    do i=1,num
       read(io,'(a)',iostat=ierr) history_files(i)
       if(ierr/=0) exit
    enddo
    close(io)

    !set number of secondary EEPs between each primary EEP
    call set_eep_interval(ierr)
    if(ierr/=0) then
       write(0,*) ' make_eeps: problem reading input.eep'
       write(0,*) '            setting default EEPs     '
    endif
    !read history file format specs

    open(unit=io,file='input.format',status='old',action='read',iostat=ierr)
    if(ierr/=0)then
       write(0,*) ' make_eeps: problem reading input.format'
       return
    endif
    read(io,*) 
    read(io,*) head
    read(io,*) main
    read(io,*) xtra
    close(io)
    call free_iounit(io)
    !set up columns to be used
    call setup_columns(history_columns_list,ierr)
  end subroutine read_input

end program make_eeps
