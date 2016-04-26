module iso_eep_support

  !MESA modules
  use const_def, only: dp, sp
  use utils_lib, only: alloc_iounit, free_iounit, has_bad_real

  implicit none

  logical, parameter :: verbose = .false.
  logical, parameter :: old_core_mass_names=.false.
  integer, parameter :: col_width = 32, file_path = 256

  !for isochrones
  integer, parameter :: age_scale_linear = 0
  integer, parameter :: age_scale_log10  = 1

  logical, parameter ::check_initial_mass = .true.
  real(dp) :: mass_eps = 1d-6

  real(dp), parameter :: ln10=log(1.0d1)
  real(sp), parameter :: ln10_sp = log(10.0)

  character(len=file_path) :: history_dir, eep_dir, iso_dir  

  !stellar types for handling primary eeps
  integer, parameter :: unknown           =  1 !for initialization only
  integer, parameter :: sub_stellar       =  2 !no fusion = brown dwarf
  integer, parameter :: star_low_mass     =  3 !ends as a WD
  integer, parameter :: star_high_mass    =  4 !does not end as a WD

  character(len=10) :: star_label(4)=['   unknown', 'substellar', '  low-mass', ' high-mass']

  !eep_controls quantities
  logical :: make_bin_tracks=.true. !faster for repeated isochrone construction
  logical :: make_bin_isos  =.true. 

  ! central limits for high- / intermediate-mass stars, set these from input eep_controls nml
  real(dp) :: center_gamma_limit=1d2 
  real(dp) :: center_carbon_limit=1d-4
  real(dp) :: log_center_T_limit=9d0
  real(dp) :: high_mass_limit = 1d1 !Msun
  real(dp) :: very_low_mass_limit = 0.5d0 !Msun

  ! default column format specs
  integer :: head !=29
  integer :: main !=28
  integer :: xtra !=0

  ! quantities from history file that need to be identified
  integer :: i_age, i_mass, i_logLH, i_logLHe, i_logTe, i_logL
  integer :: i_logg, i_Tc, i_Rhoc, i_Xc, i_Yc, i_he_core, i_co_core
  integer :: i_Cc, i_gamma, i_surfH

  ! for use when constructing EEP distance
  logical :: weight_center_rho_T_by_Xc
  real(dp) :: Teff_scale=2d0
  real(dp) :: logL_scale=0.125d0
  real(dp) :: age_scale=0.05d0
  real(dp) :: Rhoc_scale=1d0
  real(dp) :: Tc_scale=1d0

  !for columns
  integer, parameter :: max_col = 180
  integer :: ncol
  integer, parameter :: column_int=0
  integer, parameter :: column_dbl=1
  type column
     character(len=col_width) :: name
     integer :: type, loc
  end type column
  type(column), allocatable :: cols(:) !(ncol)

  !EEP arrays
  integer, parameter :: primary = 10 ! number of primary EEPs
  ! as set by primary_eep
  integer :: eep_interval(primary-1) ! number of secondary EEPs
  ! between the primaries

  !holds an evolutionary track, use an array of these for multiple tracks
  type track
     character(len=file_path) :: filename, cmd_suffix
     character(len=8) :: version_string
     type(column), allocatable :: cols(:)
     logical :: has_phase = .false., ignore=.false.
     integer :: ncol, ntrack, neep, MESA_revision_number
     integer :: star_type = unknown
     integer :: nfil !number of filters
     character(len=20), allocatable :: labels(:) !(nfil) for mags
     real(sp), allocatable :: mags(:,:) !(nfil,neep)
     integer, allocatable :: eep(:)
     real(dp) :: initial_mass, initial_Y, Fe_div_H, initial_Z, v_div_vcrit, alpha_div_Fe
     real(dp), allocatable :: tr(:,:), dist(:), phase(:)
     !these are used internally as an intermediate step
     real(dp), allocatable :: eep_tr(:,:), eep_dist(:) !(ncol,neep), (neep)
     real(sp) :: Av=0.0, Rv=3.1
  end type track

  !holds one isochrone
  type isochrone
     integer :: neep !number of eeps
     integer :: ncol !number of columns in history file
     integer :: nfil !number of filters for mags
     integer :: age_scale ! either linear or log10
     type(column), allocatable :: cols(:) !for history columns
     character(len=20), allocatable :: labels(:) !for mags
     logical :: has_phase = .false.
     integer, allocatable :: eep(:)
     real(dp) :: age, Fe_div_H, initial_Y, initial_Z, v_div_vcrit, alpha_div_Fe ! log(age in yrs)
     real(dp), allocatable :: phase(:) !neep
     real(dp), allocatable :: data(:,:) !(ncol,neep)
     real(sp), allocatable :: mags(:,:) !(num filters, neep)
     real(sp) :: Av, Rv
  end type isochrone

  !holds a set of isochrones
  type isochrone_set
     integer :: MESA_revision_number, number_of_isochrones
     integer, allocatable :: num_valid_eeps(:)
     type(isochrone), allocatable :: iso(:)
     real(sp) :: Av, Rv
     real(dp) :: initial_Y, initial_Z, Fe_div_H, v_div_vcrit, alpha_div_Fe
     character(len=file_path) :: cmd_suffix, filename
     character(len=8) :: version_string
  end type isochrone_set

contains

  include 'num_binary_search.inc'

  elemental function pow10_sg(x) result(y)
    real(sp), intent(in) :: x
    real(sp) :: y
    y = exp(ln10_sp*x)
  end function pow10_sg

  elemental function pow10(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = exp(ln10*x)
  end function pow10

  subroutine process_history_columns(history_columns_list,ierr)
    character(len=file_path), intent(in) :: history_columns_list
    integer, intent(out) :: ierr
    integer :: i, io, ncols(2), nchar, column_length, pass
    character(len=file_path) :: line, column_name
    logical :: is_int
    ierr=0
    io=alloc_iounit(ierr)
    open(io,file=trim(history_columns_list),action='read',status='old',iostat=ierr)
    if(ierr/=0) then
       write(*,*) 'failed to open history columns list: ', trim(history_columns_list)
       call free_iounit(io)
       return
    endif
    ncols=0
    do pass=1,2
       if(pass==2) allocate(cols(ncols(1)))
       inner_loop: do while(.true.)
          is_int = .false.
          read(io,'(a)',iostat=ierr) line
          if(ierr/=0) exit inner_loop

          !remove any nasty tabs
          do while(index(line,char(9))>0)
             write(*,*) ' found a tab in line: ', trim(line)
             i=index(line,char(9))
             line(i:i)=' '
          enddo

          nchar=len_trim(line)
          if(nchar==0) cycle inner_loop ! ignore blank line

          line=adjustl(line)
          i=index(line,'!')-1
          if(i<0) i=len_trim(line)

          if(i==0) then       !comment line
             if(verbose) write(*,*) ' comment: ', trim(line)
             cycle inner_loop
          else if(index(line(1:i),'number')>0) then
             if(verbose) write(*,*) '****** ', trim(line)
             if(verbose) write(*,*) '****** index of number > 0 => integer'
             is_int = .true.
          else if(index(line(1:i),'num_')==1)then
             if(verbose) write(*,*) '****** ', trim(line)
             if(verbose) write(*,*) '****** index of num_ == 1 => integer'
             is_int = .true.
          endif

          column_name = line
          ncols(pass)=ncols(pass)+1
          if(i==0) then
             column_length=len_trim(column_name)
          else
             column_length=len_trim(column_name(1:i))
          endif
          do i=1,column_length
             if(column_name(i:i)==' ') column_name(i:i)='_'
          enddo
          if(verbose) write(*,'(2i5,a32,i5)') pass, ncols(pass),trim(column_name(1:column_length)), column_length
          if(pass==2) then
             cols(ncols(pass))% name = trim(column_name(1:column_length))
             if(is_int) then
                cols(ncols(pass))% type = column_int
             else
                cols(ncols(pass))% type = column_dbl
             endif
             cols(ncols(pass))% loc = ncols(pass)
          endif
       end do inner_loop
       if(pass==1) rewind(io)
       if(pass==2) close(io)
    end do
    if(ncols(1)==ncols(2)) then
       ierr=0
       ncol=ncols(1)
    endif
    call free_iounit(io)

    if(verbose) write(*,*) 'process_history_columns: ncol = ', ncol

  end subroutine process_history_columns

  integer function locate_column(col_name)
    character(len=col_width), intent(in) :: col_name
    integer :: i
    locate_column = -1
    do i=1,size(cols)
       if(trim(cols(i)% name)==trim(col_name)) then 
          locate_column = i
          return
       endif
    enddo
    if(locate_column<0) then
       write(0,*) ' locate_column, could not find column: ', trim(col_name)
    endif
  end function locate_column

  subroutine alloc_track(filename,t)
    character(len=file_path), intent(in) :: filename
    type(track), pointer :: t
    allocate(t)
    t% ncol = ncol
    t% neep = primary
    t% filename = trim(filename)
    allocate(t% eep(t% neep)) 
  end subroutine alloc_track

  subroutine write_track(x)
    type(track), intent(in) :: x
    integer :: io, ierr, j
    character(len=8) :: have_phase
    io=alloc_iounit(ierr)
    write(*,*) '    ', trim(x% filename)
    open(io,file=trim(x% filename),action='write',status='unknown')
    if(x% has_phase)then
       have_phase = 'YES'
    else
       have_phase = 'NO'
    endif
    have_phase = adjustr(have_phase)
    write(io,'(a25,a8)') '# MIST version number  = ', x% version_string
    write(io,'(a25,i8)') '# MESA revision number = ', x% MESA_revision_number
!                      123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a88)') '#  Yinit        Zinit   [Fe/H]   [a/Fe]  v/vcrit                                        '
    write(io,'(a2,f6.4,1p1e13.5,0p3f9.2)') '# ', x% initial_Y, x% initial_Z, x% Fe_div_H, x% alpha_div_Fe, x% v_div_vcrit
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a1,1x,a16,4a8,2x,a10)') '#','initial_mass', 'N_pts', 'N_EEP', 'N_col', 'phase', 'type'
    write(io,'(a1,1x,1p1e16.10,3i8,a8,2x,a10)') '#', x% initial_mass, x% ntrack, x% neep, x% ncol, have_phase, &
         star_label(x% star_type)
    write(io,'(a8,20i8)') '# EEPs: ', x% eep
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'

    if(x% has_phase)then
       write(io,'(299(27x,i5))') (j,j=1,x% ncol+1)
       write(io,'(299a32)') adjustr(x% cols(:)% name), 'phase'
       do j=x% eep(1),x% ntrack
          write(io,'(299(1pes32.16e3))') x% tr(:,j), x% phase(j)
       enddo
    else
       write(io,'(299(27x,i5))') (j,j=1,x% ncol)
       write(io,'(299a32)') adjustr(x% cols(:)% name)
       do j=x% eep(1),x% ntrack
          write(io,'(299(1pes32.16e3))') x% tr(:,j)
       enddo
    endif
    close(io)
    call free_iounit(io)
  end subroutine write_track

  !writes a series of n isochrones to filename
  subroutine write_isochrones_to_file(set)
    type(isochrone_set), intent(in) :: set
    integer :: i, ierr, io, n
    io=alloc_iounit(ierr)
    n=set% number_of_isochrones
    write(0,*) ' isochrone output file = ', trim(set% filename)
    open(io,file=trim(set% filename),action='write',status='unknown',iostat=ierr)
    write(io,'(a25,a8)') '# MIST version number  = ', set% version_string
    write(io,'(a25,i8)') '# MESA revision number = ', set% MESA_revision_number
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a88)') '#  Yinit        Zinit   [Fe/H]   [a/Fe]  v/vcrit                                        '
    write(io,'(a2,f6.4,1p1e13.5,0p3f9.2)') '# ', set% initial_Y, set% initial_Z, set% Fe_div_H, set% alpha_div_Fe, &
         set% v_div_vcrit
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a25,i5)') '# number of isochrones = ', n
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
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
    integer :: i, my_ncol
    my_ncol = iso% ncol + 2 !add two for eep and age
    write(io,'(a25,2i5)') '# number of EEPs, cols = ', iso% neep, my_ncol
    write(io,'(a1,i4,299i32)') '#    ', (i,i=1,my_ncol)
    if(iso% age_scale==age_scale_log10)then
       write(io,'(a5,299a32)') '# EEP', 'log10_isochrone_age_yr', adjustr(iso% cols(:)% name)
    elseif(iso% age_scale==age_scale_linear)then
       write(io,'(a5,299a32)') '# EEP', 'isochrone_age_yr', adjustr(iso% cols(:)% name)
    endif
    do i=1,iso% neep
       write(io,'(i5,299(1pes32.16e3))') iso% eep(i), iso% age, iso% data(:,i)
    enddo
  end subroutine write_isochrone_to_file_orig

  subroutine write_isochrone_to_file_phase(io,iso)
    integer, intent(in) :: io
    type(isochrone), intent(in) :: iso
    integer :: i, my_ncol
    my_ncol = iso% ncol + 3 !add three for eep, phase, and age
    write(io,'(a25,2i5)') '# number of EEPs, cols = ', iso% neep, my_ncol
    write(io,'(a1,i4,299i32)') '#    ', (i,i=1,my_ncol)
    if(iso% age_scale==age_scale_log10)then
       write(io,'(a5,299a32)') '# EEP', 'log10_isochrone_age_yr', adjustr(iso% cols(:)% name), 'phase'
    elseif(iso% age_scale==age_scale_linear)then
       write(io,'(a5,299a32)') '# EEP', 'isochrone_age_yr', adjustr(iso% cols(:)% name), 'phase'
    endif
    do i=1,iso% neep
       write(io,'(i5,299(1pes32.16e3))') iso% eep(i), iso% age, &
            iso% data(:,i), real(iso% phase(i),kind=dp)
    enddo
  end subroutine write_isochrone_to_file_phase


  subroutine read_eep(x,full_path)
    type(track), intent(inout) :: x
    logical, optional :: full_path
    logical :: use_full_path
    integer :: ierr, io, j
    character(len=8) :: phase_info
    character(len=file_path) :: eepfile
    character(len=10) :: type_label
    io=alloc_iounit(ierr)

    if(present(full_path))then
       use_full_path = full_path
    else
       use_full_path = .false.
    endif

    if(use_full_path)then
       eepfile = trim(x% filename)
    else
       eepfile = trim(eep_dir) // '/' // trim(x% filename) // '.eep'
    endif

    open(io,file=trim(eepfile),status='old',action='read',iostat=ierr)

    !check if the file was opened successfully; if not, then fail
    if(ierr/=0) then
       x% ignore=.true.
       write(*,*) '  PROBLEM OPENING EEP FILE: ', trim(eepfile)
       close(io)
       call free_iounit(io)
       return
    endif

    read(io,'(25x,a8)') x% version_string
    read(io,'(25x,i8)') x% MESA_revision_number
    read(io,*) !comment line
    read(io,*) !comment line
    read(io,'(2x,f6.4,1p1e13.5,0p3f9.2)') x% initial_Y, x% initial_Z, x% Fe_div_H, x% alpha_div_Fe, x% v_div_vcrit
    read(io,*) !comment line
    read(io,*) !comment line
    read(io,'(2x,1p1e16.10,3i8,a8,2x,a10)') x% initial_mass, x% ntrack, x% neep, x% ncol, phase_info, type_label

    call set_star_type_from_label(type_label,x)

    allocate(x% tr(x% ncol, x% ntrack), x% eep(x% neep), x% cols(x% ncol))
    if(index(phase_info,'YES')/=0) then
       x% has_phase = .true.
       allocate(x% phase(x% ntrack))
    else
       x% has_phase = .false.
    endif
    read(io,'(8x,299i8)') x% eep
    read(io,*) ! comment line
    read(io,*) ! column numbers
    if(x% has_phase) then
       read(io,'(299a32)') x% cols(:)% name ! column names
       do j=x% eep(1), x% ntrack
          read(io,'(299(1pes32.16e3))')  x% tr(:,j), x% phase(j)
       enddo
    else
       read(io,'(299a32)') x% cols
       do j=x% eep(1), x% ntrack
          read(io,'(299(1pes32.16e3))') x% tr(:,j)
       enddo
    endif
    close(io)
    call free_iounit(io)
  end subroutine read_eep

  subroutine read_history_file(t,ierr)
    type(track), intent(inout) :: t
    integer, intent(out) :: ierr
    character(len=8192) :: line
    character(len=file_path) :: binfile
    integer :: i, ilo, ihi, io, j, imass, iversion
    integer, allocatable :: output(:) !ncol
    logical :: binfile_exists

    ierr = 0
    if(verbose)then
       write(*,*)  '    main = ', main
       write(*,*)  '    head = ', head
       write(*,*)  '    xtra = ', xtra
    endif

    ! using unformatted binary files makes the process of creating
    ! EEP files much faster. so check to see if the .bin exists and,
    ! if it does, read it and be done. otherwise read the .data
    ! file and write a new .bin at the end.
    ! time goes from t=2min to t<3sec for 94 tracks. tight!
    if(make_bin_tracks)then
       binfile=trim(history_dir) // '/' // trim(t% filename) // '.bin'
       inquire(file=binfile,exist=binfile_exists)

       if(binfile_exists)then
          call read_history_bin(t)
          call distance_along_track(t)
          return
       endif
    endif

    ! if the binfile does not exist, then we read the history files and write new
    ! .bins.  slow.
    io=alloc_iounit(ierr)
    open(unit=io,file=trim(trim(history_dir) // '/' // t% filename),status='old',action='read')
    !read first 3 lines of header
    !currently don't use all of this info, but could...
    imass=0
    iversion=0
    do j=1,3
       read(io,'(a)') line
       do i=1,7
          ilo =   1 + head*(i-1) + xtra*(i-1)
          ihi = ilo + head-1
          if(j==2)then
             if(adjustl(adjustr(line(ilo:ihi)))=='version_number') iversion=i
             if(adjustl(adjustr(line(ilo:ihi)))=='initial_mass')   imass=i
          else if(j==3)then
             if(i==iversion) read(line(ilo:ihi),*) t% MESA_revision_number
             if(i==imass) read(line(ilo:ihi),*) t% initial_mass
          endif
       enddo
    enddo

    read(io,*) !blank line

    if(verbose) write(*,*) trim(t% filename), t% initial_Mass, t% MESA_revision_number

    !read first two lines of main section
    read(io,*)
    read(io,'(a)') line
    call dict(line,output)

    !figure out how many data lines
    j=1
    do while(.true.)
       read(io,*,iostat=ierr)
       if(ierr/=0) exit
       j=j+1
    enddo
    t% ntrack = j-1

    t% ncol = ncol
    allocate(t% tr(t% ncol, t% ntrack),t% cols(t% ncol))
    allocate(t% eep_tr(t% ncol, primary), t% eep_dist(primary))

    t% cols = cols

    !ignore file header, already read it once
    rewind(io)
    do i=1,6
       read(io,*)
    enddo

    !read track data
    do j=1,t% ntrack
       read(io,'(a)',iostat=ierr) line
       if(ierr/=0) exit
       do i=1,t% ncol
          if(output(i)/=0)then
             ilo =   1 + main*(output(i)-1)+xtra*(output(i)-1)
             ihi = ilo + main-1
             read(line(ilo:ihi),*) t% tr(i,j)
          endif
       enddo
    enddo
    close(io)

    call free_iounit(io)

    !compute distance along track
    allocate(t% dist(t% ntrack))
    call distance_along_track(t)
    call set_star_type_from_history(t)

    !finally check if initial mass is correct and, if not, replace it
    if(check_initial_mass) then
       if(abs(t% initial_mass - t% tr(i_mass,1)) > mass_eps) t% initial_mass = t% tr(i_mass,1)
    endif

    if(make_bin_tracks) call write_history_bin(t)

  end subroutine read_history_file

  subroutine read_history_bin(t)
    type(track), intent(inout) :: t
    integer :: io, ierr
    character(len=file_path) :: binfile
    io=alloc_iounit(ierr)
    binfile = trim(history_dir) // '/' // trim(t% filename) // '.bin'
    open(io,file=trim(binfile),form='unformatted',status='old',action='read')
    read(io) t% filename
    read(io) t% ncol, t% ntrack, t% neep, t% MESA_revision_number, t% star_type
    read(io) t% initial_mass
    allocate(t% tr(t% ncol, t% ntrack),t% cols(t% ncol),t% dist(t% ntrack))
    allocate(t% eep_tr(t% ncol, t% neep), t% eep_dist(t% neep))
    read(io) t% cols
    read(io) t% tr
    read(io) t% dist
    close(io) 
    t% eep_tr = 0d0
    t% eep_dist = 0d0
    call free_iounit(io)

  end subroutine read_history_bin

  subroutine write_history_bin(t)
    type(track), intent(in) :: t
    integer :: io, ierr
    character(len=file_path) :: binfile
    io=alloc_iounit(ierr)
    binfile = trim(history_dir) // '/' // trim(t% filename) // '.bin'
    open(io,file=trim(binfile),form='unformatted',action='write',status='unknown')
    write(io) t% filename
    write(io) t% ncol, t% ntrack, t% neep, t% MESA_revision_number, t% star_type
    write(io) t% initial_mass
    write(io) t% cols
    write(io) t% tr
    write(io) t% dist
    close(io) 
    call free_iounit(io)
  end subroutine write_history_bin


  subroutine read_isochrone_file(s,ierr)
    type(isochrone_set), intent(inout) :: s
    integer, intent(out) :: ierr
    integer :: io, i, n
    ierr=0
    io=alloc_iounit(ierr)
    open(io,file=trim(s% filename), action='read', status='old')
    read(io,'(25x,a8)') s% version_string
    read(io,'(25x,i8)') s% MESA_revision_number
    read(io,*) !comment line
    read(io,*) !comment line
    read(io,'(2x,f6.4,1p1e13.5,0p3f9.2)') s% initial_Y, s% initial_Z, s% Fe_div_H, s% alpha_div_Fe, s% v_div_vcrit
    read(io,*) !comment line
    read(io,'(25x,i5)') s% number_of_isochrones
    read(io,*) !comment line

    !make room
    n=s% number_of_isochrones
    allocate(s% iso(n))
    !read the data
    do i=1,n
       call read_one_isochrone_from_file(io,s% iso(i))
       s% iso(i)% Fe_div_H = s% Fe_div_H
       s% iso(i)% initial_Y = s% initial_Y
       s% iso(i)% initial_Z = s% initial_Z
       if(i<n) read(io,*)
       if(i<n) read(io,*)
    enddo    
    close(io)
    call free_iounit(io)
  end subroutine read_isochrone_file

  
  subroutine read_one_isochrone_from_file(io,iso)
    integer, intent(in) :: io
    type(isochrone), intent(out) :: iso
    integer :: i, my_ncol
    type(column), allocatable :: cols(:)
    read(io,'(25x,2i5)') iso% neep, my_ncol
    read(io,*) !skip column numbers
    allocate(cols(my_ncol))
    read(io,'(2x,a3,299a32)') cols(:)% name
    
    if(index(cols(2)% name, 'log10')>0)then
       iso% age_scale = age_scale_log10
    else
       iso% age_scale = age_scale_linear
    endif

    iso% has_phase = index(cols(my_ncol)% name, 'phase') > 0 

    if(iso% has_phase)then
       iso% ncol = my_ncol - 3
    else
       iso% ncol = my_ncol - 2
    endif
    
    allocate(iso% cols(iso% ncol))
    iso% cols(:)% name = cols(3:iso% ncol+2)% name

    if(iso% has_phase) allocate(iso% phase(iso% neep))

    allocate(iso% eep(iso% neep), iso% data(iso% ncol, iso% neep))
    do i=1,iso% neep
       if(iso% has_phase)then
          read(io,'(i5,299(1pes32.16e3))') iso% eep(i), iso% age, iso% data(:,i), iso% phase(i)
       else
          read(io,'(i5,299(1pes32.16e3))') iso% eep(i), iso% age, iso% data(:,i)
       endif
    enddo
    
  end subroutine read_one_isochrone_from_file

  subroutine distance_along_track(t)
    type(track), intent(inout) :: t
    real(dp) :: tmp_dist, weight, max_center_h1
    integer :: j

    if(weight_center_rho_T_by_Xc)then
       max_center_h1 = maxval(t% tr(i_Xc,:))
       if(max_center_h1 <= 0d0) max_center_h1 = 1d0
    else
       max_center_h1 = 1d0
       weight = 1d0
    endif

    t% dist(1) = 0d0
    if(t% ntrack > 3)then
       do j = 2, t% ntrack
          
          if(weight_center_rho_T_by_Xc)then
             weight = max(0d0, t% tr(i_Xc,j)/max_center_h1)
          endif
          
          !build up the distance between EEPs piece by piece
          tmp_dist =            Teff_scale*sqdiff(t% tr(i_logTe,j) , t% tr(i_logTe,j-1))
          tmp_dist = tmp_dist + logL_scale*sqdiff(t% tr(i_logL, j) , t% tr(i_logL, j-1))
          tmp_dist = tmp_dist + weight * Rhoc_scale * sqdiff(t% tr(i_Rhoc, j) , t% tr(i_Rhoc, j-1))
          tmp_dist = tmp_dist + weight * Tc_scale*  sqdiff(t% tr(i_Tc,   j) , t% tr(i_Tc,   j-1))
          tmp_dist = tmp_dist + age_scale* sqdiff(log10(t% tr(i_age,j)) , log10(t% tr(i_age,j-1))) 

          t% dist(j) = t% dist(j-1) + sqrt(tmp_dist)
       enddo
    endif
  end subroutine distance_along_track

  elemental function sqdiff(x0,x1) result(y) !square of x, y=x*x
    real(dp), intent(in) :: x0, x1
    real(dp) :: y, dx
    dx=x0-x1
    y = dx*dx
  end function sqdiff

  subroutine dict(input,output)
    character(len=*), intent(in) :: input
    integer, allocatable, intent(out) :: output(:)
    integer :: i,ihi,ilo,j
    logical :: have_col(ncol)
    have_col = .false.
    allocate(output(ncol))
    output = 0
    if (verbose) write(*,*) 'number of columns: ', ncol
    do j=1,ncol
       iloop: do i=1,max_col
          ilo =   1 + main*(i-1) + xtra*(i-1)
          ihi = ilo + main-1
          if(adjustl(adjustr(input(ilo:ihi)))==trim(cols(j)% name))then
             output(j)=i
             have_col(j) = .true.
             exit iloop
          endif
       enddo iloop
    enddo
    do j=1,ncol
       if(.not.have_col(j)) write(*,*) 'do not have ', trim(cols(j)% name)
    enddo
  end subroutine dict

  subroutine set_eep_interval(ierr)
    integer, intent(out) :: ierr
    integer :: i, j, io, my_eep_interval, my_num_eep
    io=alloc_iounit(ierr)
    open(unit=io,file='input.eep',status='old',action='read',iostat=ierr)
    if(ierr/=0) then
       call set_default_eep_interval
       return
    endif
    read(io,*) !comments
    read(io,'(i3)') my_num_eep
    if(my_num_eep /= (primary-1)) then
       write(0,*) 'incorrect number of eep intervals in input.eep'
       call set_default_eep_interval
       return
    endif
    do i=1,primary-1
       read(io,'(i3,i5)') j, my_eep_interval
       eep_interval(j) = my_eep_interval
    enddo
    close(io)
  end subroutine set_eep_interval

  subroutine set_default_eep_interval
    !determines total number and relative density of EEPs
    !along an evolutionary track. Total number of EEPs
    !will be sum(eep_interval) + primary.
    eep_interval = 50
  end subroutine set_default_eep_interval

  subroutine setup_columns(history_columns_list,ierr)
    !reads a history_columns.list file to determine what columns to write
    !to the .eep files; it identifies those important columns that are 
    !required to identify eeps in the evolutionary tracks
    character(len=file_path) :: history_columns_list
    integer, intent(out) :: ierr
    character(len=col_width) :: col_name
    call process_history_columns(history_columns_list,ierr)
    if(ierr/=0) then
       write(*,*) 'failed in process_history_columns'
       return
    endif
    ncol = size(cols) 
    if(verbose) write(*,*) ' number of history columns = ', ncol
    col_name = 'star_age'; i_age = locate_column(col_name)
    col_name = 'star_mass'; i_mass= locate_column(col_name)
    col_name='log_LH'; i_logLH=locate_column(col_name)
    col_name='log_LHe'; i_logLHe=locate_column(col_name)
    col_name='log_Teff'; i_logTe=locate_column(col_name)
    col_name='log_L'; i_logL=locate_column(col_name)
    col_name='log_g'; i_logg=locate_column(col_name)
    col_name='log_center_T'; i_Tc=locate_column(col_name)
    col_name='log_center_Rho'; i_Rhoc=locate_column(col_name)
    col_name='center_h1'; i_Xc=locate_column(col_name)
    col_name='center_he4'; i_Yc=locate_column(col_name)
    col_name='center_c12'; i_Cc=locate_column(col_name)
    col_name='center_gamma'; i_gamma=locate_column(col_name)
    col_name='surface_h1'; i_surfH=locate_column(col_name)
    if(old_core_mass_names)then
       col_name='h1_boundary_mass'; i_he_core = locate_column(col_name)
       col_name='he4_boundary_mass'; i_co_core = locate_column(col_name)
    else
       col_name='he_core_mass'; i_he_core = locate_column(col_name)
       col_name='c_core_mass'; i_co_core = locate_column(col_name)
    endif
    if(verbose)then
       write(*,*) ' star_age column = ', i_age
       write(*,*) ' star_mass column= ', i_mass      
    endif
  end subroutine setup_columns

  subroutine set_star_type_from_label(label,t)
    character(len=10), intent(in) :: label
    type(track), intent(inout) :: t
    integer :: n,i
    n=size(star_label)
    do i=1,n
       if(label==star_label(i)) t% star_type = i
    enddo
  end subroutine set_star_type_from_label

  subroutine set_star_type_from_history(t)
    type(track), intent(inout) :: t
    integer :: n

    n=t% ntrack

    !simple test for substellar is that central H is unchanged
    if( t% initial_mass < 0.1d0 .and. maxval( abs(t% tr(i_Xc,:) - t% tr(i_Xc,1)) ) < 0.1d0)then
       t% star_type = sub_stellar
       return
    endif

    !only reach center_gamma_limit if the star evolves to a WD
    if( t% tr(i_gamma,n) > center_gamma_limit) then
       t% star_type = star_low_mass
       return
    endif

    !simple test for high-mass stars is that central C is depleted
    if(maxval(t% tr(i_Cc,:)) > 0.4d0 .and. t% tr(i_Cc,n) < center_carbon_limit)then
       t% star_type = star_high_mass
       return
    endif

    !alternative test for high-mass stars is that they reach a
    !central temperature threshhold
    if(t% tr(i_Tc,n) > log_center_T_limit)then
       t% star_type = star_high_mass
    else
       t% star_type = star_low_mass
    endif

    !last gasp test for high-mass stars is the initial mass...
    if(t% initial_mass >= high_mass_limit) then
       t% star_type = star_high_mass
    else
       t% star_type = star_low_mass
    endif
  end subroutine set_star_type_from_history

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
   
end module iso_eep_support
