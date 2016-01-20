module iso_eep_support

  !MESA modules
  use const_def, only: dp, sp
  use utils_lib, only: alloc_iounit, free_iounit

  implicit none

  logical, parameter :: verbose = .false.
  logical, parameter :: old_core_mass_names=.false.
  integer, parameter :: col_width = 32, file_path = 256


  logical, parameter ::check_initial_mass = .true.
  real(dp) :: mass_eps = 1d-6

  character(len=file_path) :: history_dir, eep_dir, iso_dir  

  !stellar types for handling primary eeps
  integer, parameter :: unknown           =  1 !for initialization only
  integer, parameter :: sub_stellar       =  2 !no fusion = brown dwarf
  integer, parameter :: star_low_mass     =  3 !ends as a WD
  integer, parameter :: star_high_mass    =  4 !does not end as a WD

  character(len=10) :: star_label(4)=['   unknown', 'substellar', '  low-mass', ' high-mass']

  ! central gamma limit for high- / intermediate-mass stars
  real(dp) :: center_gamma_limit=1d2 
  real(dp) :: center_carbon_limit=1d-4
  real(dp) :: log_center_T_limit=9d0
  real(dp) :: min_for_high_mass_star=1d1 !Msun

  ! format specs
  integer :: head !=29
  integer :: main !=28
  integer :: xtra !=0

  ! quantities from history file that need to be identified
  integer :: i_age, i_mass, i_logLH, i_logLHe, i_logTe, i_logL
  integer :: i_logg, i_Tc, i_Rhoc, i_Xc, i_Yc, i_he_core, i_co_core
  integer :: i_Cc, i_gamma, i_surfH

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
     character(len=file_path) :: filename
     type(column), allocatable :: cols(:)
     logical :: has_phase = .false., ignore=.false.
     integer :: ncol, ntrack, neep, version_number
     integer :: star_type = unknown
     integer, allocatable :: eep(:)
     real(dp) :: initial_mass
     real(dp), allocatable :: tr(:,:), dist(:), phase(:)
     !these are used internally as an intermediate step
     real(dp), allocatable :: eep_tr(:,:), eep_dist(:) !(ncol,neep), (neep)
  end type track

  !holds one isochrone
  type isochrone
     integer :: neep !number of eeps
     integer :: ncol !number of columns in history file
     integer :: nfil !number of filters for mags
     type(column), allocatable :: cols(:) !for history columns
     character(len=20), allocatable :: labels(:) !for mags
     logical :: has_phase = .false.
     integer, allocatable :: eep(:)
     real(dp) :: age ! log(age in yrs)
     real(dp), allocatable :: phase(:) !neep
     real(dp), allocatable :: data(:,:) !(ncol,neep)
     real(sp), allocatable :: mags(:,:) !(num filters, neep)
     real(sp) :: Av, Rv
  end type isochrone

  !holds a set of isochrones
  type isochrone_set
     integer :: version_number
     integer, allocatable :: num_valid_eeps(:)
     type(isochrone), allocatable :: iso(:)
     real(sp) :: Av, Rv
     character(len=file_path) :: cmd_suffix, filename
  end type isochrone_set

contains

  include 'num_binary_search.inc'

  elemental function pow10(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    real(dp), parameter :: ln10=2.3025850929940459d0
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
    if(x% has_phase)then
       call write_track_phase(x)
    else
       call write_track_orig(x)
    endif
  end subroutine write_track

  subroutine write_track_orig(x)
    type(track), intent(in) :: x
    integer :: io,ierr,j
    io=alloc_iounit(ierr)
    write(*,*) '    ', trim(x% filename)
    open(io,file=trim(x% filename),action='write',status='unknown')
    write(io,'(a20,5a8,2x,a10)') 'initial_mass', 'N_pts', 'N_EEP', 'N_col', 'version', 'phase', 'type'
    write(io,'(1p1e20.10,4i8,a8,2x,a10)') x% initial_mass, x% ntrack, x% neep, x% ncol, &
         x% version_number, 'NO', star_label(x% star_type)
    write(io,'(a10,20i8)') '   EEPs:  ', x% eep
    write(io,'(299(27x,i5))') (j,j=1,x% ncol)
    write(io,'(299a32)') adjustr(x% cols(:)% name)
    do j=x% eep(1),x% ntrack
       write(io,'(299(1pes32.16e3))') x% tr(:,j)
    enddo
    close(io)
    call free_iounit(io)
  end subroutine write_track_orig

  subroutine write_track_phase(x)
    type(track), intent(in) :: x
    integer :: io,ierr,j
    io=alloc_iounit(ierr)
    write(*,*) '    ', trim(x% filename)
    open(io,file=trim(x% filename),action='write',status='unknown')
    write(io,'(a20,5a8,2x,a10)') 'initial_mass', 'N_pts', 'N_EEP', 'N_col', 'version', 'phase', 'type'
    write(io,'(1p1e20.10,4i8,a8,2x,a10)') x% initial_mass, x% ntrack, x% neep, x% ncol, & 
         x% version_number, 'YES', star_label(x% star_type)
    write(io,'(a10,20i8)') '   EEPs:  ', x% eep
    write(io,'(299(27x,i5))') (j,j=1,x% ncol+1)
    write(io,'(299a32)') adjustr(x% cols(:)% name), 'phase'
    do j=x% eep(1),x% ntrack
       write(io,'(299(1pes32.16e3))') x% tr(:,j), x% phase(j)
    enddo
    close(io)
    call free_iounit(io)
  end subroutine write_track_phase

  subroutine read_eep(x)
    type(track), intent(inout) :: x
    integer :: ierr, io, j
    character(len=8) :: phase_info
    character(len=file_path) :: eepfile
    character(len=10) :: type_label
    io=alloc_iounit(ierr)
    eepfile = trim(eep_dir) // '/' // trim(x% filename) // '.eep'
    open(io,file=trim(eepfile),status='old',action='read',iostat=ierr)

    !check if the file was opened successfully; if not, then fail
    if(ierr/=0) then
       x% ignore=.true.
       write(*,*) '  PROBLEM OPENING EEP FILE: ', trim(eepfile)
       close(io)
       call free_iounit(io)
       return
    endif

    read(io,*)
    read(io,'(1p1e20.10,4i8,a8,2x,a10)') x% initial_mass, x% ntrack, x% neep, x% ncol, x% version_number, phase_info, type_label

    call set_star_type_from_label(type_label,x)

    allocate(x% tr(x% ncol, x% ntrack), x% eep(x% neep), x% cols(x% ncol))
    if(index(phase_info,'YES')/=0) then
       x% has_phase = .true.
       allocate(x% phase(x% ntrack))
    else
       x% has_phase = .false.
    endif
    read(io,'(10x,299i8)') x% eep
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
    binfile=trim(history_dir) // '/' // trim(t% filename) // '.bin'
    inquire(file=binfile,exist=binfile_exists)

    if(binfile_exists)then
       call read_history_bin(t)
       call distance_along_track(t)
       return
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
             if(i==iversion) read(line(ilo:ihi),*) t% version_number
             if(i==imass) read(line(ilo:ihi),*) t% initial_mass
          endif
       enddo
    enddo

    read(io,*) !blank line

    if(verbose) write(*,*) trim(t% filename), t% initial_Mass, t% version_number

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

    call write_history_bin(t)

  end subroutine read_history_file

  subroutine read_history_bin(t)
    type(track), intent(inout) :: t
    integer :: io, ierr
    character(len=file_path) :: binfile
    io=alloc_iounit(ierr)
    binfile = trim(history_dir) // '/' // trim(t% filename) // '.bin'
    open(io,file=trim(binfile),form='unformatted',status='old',action='read')
    read(io) t% filename
    read(io) t% ncol, t% ntrack, t% neep, t% version_number, t% star_type
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
    write(io) t% ncol, t% ntrack, t% neep, t% version_number, t% star_type
    write(io) t% initial_mass
    write(io) t% cols
    write(io) t% tr
    write(io) t% dist
    close(io) 
    call free_iounit(io)
  end subroutine write_history_bin

  subroutine distance_along_track(t)
    type(track), intent(inout) :: t
    real(dp), parameter :: Teff_scale=2d0
    real(dp), parameter :: logL_scale=0.125d0
    real(dp), parameter :: age_scale=0.05d0
    real(dp), parameter :: Rhoc_scale=0.01d0
    real(dp), parameter :: Tc_scale=0.01d0
    integer :: j

    t% dist(1) = 0d0
    if(t% ntrack > 3)then
       do j = 2, t% ntrack
          t% dist(j) = t% dist(j-1) + sqrt( &
                 Teff_scale*sqdiff(t% tr(i_logTe,j) , t% tr(i_logTe,j-1))  &
               + logL_scale*sqdiff(t% tr(i_logL, j) , t% tr(i_logL, j-1))  &
               + Rhoc_scale*sqdiff(t% tr(i_Rhoc, j) , t% tr(i_Rhoc, j-1))  &
               + Tc_scale*  sqdiff(t% tr(i_Tc,   j) , t% tr(i_Tc,   j-1))  &
               + age_scale* sqdiff(log10(t% tr(i_age,j)) , log10(t% tr(i_age,j-1)))  &
               )
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
    if( maxval( abs(t% tr(i_Xc,:) - t% tr(i_Xc,1)) ) < 0.1d0)then
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
    if(t% initial_mass >= min_for_high_mass_star) then
       t% star_type = star_high_mass
    else
       t% star_type = star_low_mass
    endif
  end subroutine set_star_type_from_history
    
end module iso_eep_support
