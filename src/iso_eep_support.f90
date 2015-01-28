      module iso_eep_support

      !MESA modules
      use const_def, only: dp, sp
      use utils_lib, only: alloc_iounit, free_iounit
      
      implicit none

      integer, parameter :: col_width = 32, file_path = 256
      logical, parameter :: verbose = .false., old_core_mass_names=.false.

      logical, parameter ::check_initial_mass = .true.
      real(dp) :: mass_eps = 1d-6

      character(len=file_path) :: history_dir, eep_dir, iso_dir

      ! maximum number of columns in history file
      integer, parameter :: max_col = 200

      ! format specs
      integer :: head !=29
      integer :: main !=28
      integer :: xtra !=0

      ! quantities from history file that need to be identified
      integer :: i_age, i_mass, i_logLH, i_logLHe, i_logTe, i_logL
      integer :: i_logg, i_Tc, i_Rhoc, i_Xc, i_Yc, i_he_core, i_co_core

      !for columns
      integer :: ncol
      character(len=col_width), pointer :: cols(:) !(ncol)

      !EEP arrays
      integer, parameter :: primary = 11 ! number of primary EEPs
                                         ! as set by primary_eep
      integer :: eep_interval(primary-1) ! number of secondary EEPs
                                         ! between the primaries

      !holds an evolutionary track
      type track
         character(len=file_path) :: filename
         character(len=col_width), pointer :: cols(:)
         logical :: has_phase = .false., ignore=.false.
         integer :: ncol, ntrack, neep, version_number
         integer, allocatable :: eep(:)
         real(dp) :: initial_mass
         real(dp), allocatable :: tr(:,:), dist(:), phase(:)
      end type track

      !holds one isochrone
      type isochrone
      integer :: neep !number of eeps
      integer :: ncol !number of columns in history file
      integer :: nfil !number of filters for mags
      character(len=col_width), pointer :: cols(:)
      logical :: has_phase = .false.
      integer, allocatable :: eep(:)
      real(dp) :: age ! log(age in yrs)
      real(dp), allocatable :: phase(:) !neep
      real(dp), allocatable :: data(:,:) !(ncol,neep)
      real(sp), allocatable :: mags(:,:) !(num filters, neep)
      end type isochrone

      !holds a set of isochrones
      type isochrone_set
      integer :: version_number
      integer, allocatable :: num_valid_eeps(:)
      type(isochrone), allocatable :: iso(:)
      end type isochrone_set

      contains

      include 'num_binary_search.inc'

      subroutine process_history_columns(history_columns_list,col,ierr)
      character(len=file_path), intent(in) :: history_columns_list
      character(len=col_width), pointer, intent(out) :: col(:)
      integer, intent(out) :: ierr
      integer :: i, io, ncol(2), nchar, column_length, pass
      character(len=file_path) :: line, column_name
      ierr=0
      io=alloc_iounit(ierr)
      open(io,file=trim(history_columns_list),action='read',status='old',iostat=ierr)
      if(ierr/=0) then
         write(*,*) 'failed to open history columns list: ', trim(history_columns_list)
         call free_iounit(io)
         return
      endif
      ncol=0
      do pass=1,2
         if(pass==2) allocate(col(ncol(1)))
         inner_loop: do while(.true.)
            read(io,'(a)',iostat=ierr) line
            if(ierr/=0) exit inner_loop
            nchar=len_trim(line)
            if(nchar==0) cycle inner_loop ! ignore blank line

            line=adjustl(line)
            i=index(line,'!')-1
            if(i<0) i=len_trim(line)
            
            if(i==0) then       ! comment line
               if(verbose) write(*,*) ' comment: ', trim(line)
               cycle inner_loop
            else if(index(line(1:i),'number')>0) then
               if(verbose) write(*,*) '****** ', trim(line)
               if(verbose) write(*,*) '****** index of number > 0, cycling'
               cycle inner_loop !ignore integers
            else if(index(line(1:i),'num_')==1)then
               if(verbose) write(*,*) '****** ', trim(line)
               if(verbose) write(*,*) '****** index of num_ == 1, cycling'
               cycle inner_loop !ignore integers
            endif
            column_name = line
            ncol(pass)=ncol(pass)+1
            if(i==0) then
               column_length=len_trim(column_name)
            else
               column_length=len_trim(column_name(1:i))
            endif
            do i=1,column_length
               if(column_name(i:i)==' ') column_name(i:i)='_'
            enddo
            if(verbose) write(*,'(2i5,a32,i5)') pass, ncol(pass),trim(column_name(1:column_length)), column_length
            if(pass==2) then
               col(ncol(pass)) = trim(column_name(1:column_length))
            endif
         end do inner_loop
         if(pass==1) rewind(io)
         if(pass==2) close(io)
      end do 
      if(ncol(1)==ncol(2)) ierr=0
      call free_iounit(io)
      end subroutine process_history_columns

      integer function locate_column(col_name,columns)
      character(len=col_width), intent(in) :: col_name, columns(:)
      integer :: i
      locate_column = -1
      do i=1,size(columns)
         if(trim(columns(i))==trim(col_name)) then 
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
      open(io,file=trim(x% filename))
      write(io,'(a20,5a8)') 'initial_mass', 'N_pts', 'N_EEP', 'N_col', 'version', 'phase'
      write(io,'(1p1e20.10,4i8,a8)') x% initial_mass, x% ntrack, x% neep, x% ncol, x% version_number, 'NO'
      write(io,'(a10,20i8)') '   EEPs:  ', x% eep
      write(io,'(299(27x,i5))') (j,j=1,x% ncol + 1)
      write(io,'(299a32)') adjustr(x% cols), 'distance'
      do j=x% eep(1),x% ntrack
         write(io,'(299(1pes32.16e3))') x% tr(:,j), x% dist(j)
      enddo
      close(io)
      call free_iounit(io)
      end subroutine write_track_orig

      subroutine write_track_phase(x)
      type(track), intent(in) :: x
      integer :: io,ierr,j
      io=alloc_iounit(ierr)
      write(*,*) '    ', trim(x% filename)
      open(io,file=trim(x% filename))
      write(io,'(a20,5a8)') 'initial_mass', 'N_pts', 'N_EEP', 'N_col', 'version', 'phase'
      write(io,'(1p1e20.10,4i8,a8)') x% initial_mass, x% ntrack, x% neep, x% ncol, x% version_number, 'YES'
      write(io,'(a10,20i8)') '   EEPs:  ', x% eep
      write(io,'(299(27x,i5))') (j,j=1,x% ncol+2)
      write(io,'(299a32)') adjustr(x% cols), 'distance', 'phase'
      do j=x% eep(1),x% ntrack
         write(io,'(299(1pes32.16e3))') x% tr(:,j), x% dist(j), x% phase(j)
      enddo
      close(io)
      call free_iounit(io)
      end subroutine write_track_phase

      subroutine read_eep(x)
      type(track), intent(inout) :: x
      integer :: ierr, io, j
      character(len=8) :: phase_info
      character(len=file_path) :: eepfile
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
      read(io,'(1p1e20.10,4i8,a8)') x% initial_mass, x% ntrack, x% neep, x% ncol, x% version_number, phase_info
      !if(x% ncol /= ncol) write(*,*) '  WARNING: NCOL != DEFAULT  '
      allocate(x% tr(x% ncol, x% ntrack), x% dist(x% ntrack), x% eep(x% neep), x% cols(x% ncol))
      if(index(phase_info,'YES')/=0) then
         x% has_phase = .true.
         allocate(x% phase(x% ntrack))
      else
         x% has_phase = .false.
      endif
      read(io,'(10x,299i8)') x% eep
      read(io,*) ! column numbers
      if(x% has_phase) then
         read(io,'(299a32)') x% cols ! column names
         do j=x% eep(1), x% ntrack
            read(io,'(299(1pes32.16e3))')  x% tr(:,j), x% dist(j), x% phase(j)
         enddo
      else
         read(io,'(299a32)') x% cols
         do j=x% eep(1), x% ntrack
            read(io,'(299(1pes32.16e3))') x% tr(:,j), x% dist(j)
         enddo
      endif
      close(io)
      call free_iounit(io)
      end subroutine read_eep

      subroutine read_history_file(t)
      type(track), intent(inout) :: t
      character(len=8192) :: line
      character(len=file_path) :: binfile
      character(len=3) :: type_string(2)
      integer :: i, ilo, ihi, io, j, imass, iversion
      integer :: ierr
      integer, pointer :: output(:) !ncol
      logical :: binfile_exists

      ierr = 0
      if(verbose)then
         write(*,*)  '    main = ', main
         write(*,*)  '    head = ', head
         write(*,*)  '    xtra = ', xtra
      endif

      type_string = (/ 'flt', 'int' /)
      
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

      ! if the binfile does not exist, then we read the .data files and write new
      ! .bins.  slow.
      io=alloc_iounit(ierr)
      open(unit=io,file=trim(trim(history_dir) // '/' // t% filename),status='old')
      !read first 3 lines of header
      !currently don't use all of this, but could...
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
      call dict(line,ncol,cols,output)

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

      t% cols = cols(:)

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

      deallocate(output)

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
      open(io,file=trim(binfile),form='unformatted')
      read(io) t% filename
      read(io) t% ncol, t% ntrack, t% neep, t% version_number
      read(io) t% initial_mass
      allocate(t% tr(t% ncol, t% ntrack),t% cols(t% ncol),t% dist(t% ntrack))
      read(io) t% cols
      read(io) t% tr
      read(io) t% dist
      close(io) 
      call free_iounit(io)
      end subroutine read_history_bin

      subroutine write_history_bin(t)
      type(track), intent(in) :: t
      integer :: io, ierr
      character(len=file_path) :: binfile
      io=alloc_iounit(ierr)
      binfile = trim(history_dir) // '/' // trim(t% filename) // '.bin'
      open(io,file=trim(binfile),form='unformatted')
      write(io) t% filename
      write(io) t% ncol, t% ntrack, t% neep, t% version_number
      write(io) t% initial_mass
      write(io) t% cols
      write(io) t% tr
      write(io) t% dist
      close(io) 
      call free_iounit(io)
      end subroutine write_history_bin

      subroutine distance_along_track(t)
      type(track), intent(inout) :: t

      !real(dp), parameter :: Teff_scale=1d1, logL_scale=1.25d0
      !real(dp), parameter :: age_scale=4d0, Rhoc_scale=0d0, Tc_scale=0d0

      real(dp), parameter :: Teff_scale=1d2, logL_scale=12.5d0
      real(dp), parameter :: age_scale=4d0, Rhoc_scale=0d0, Tc_scale=0d0

      integer :: j
      
      t% dist(1) = 0d0
      if(t% ntrack > 1)then
         do j=2, t% ntrack
            t% dist(j) = t% dist(j-1) + sqrt( &
                                        Teff_scale*(t% tr(i_logTe,j) - t% tr(i_logTe,j-1))**2  &
                                      + logL_scale*(t% tr(i_logL, j) - t% tr(i_logL, j-1))**2  &
                                      + Rhoc_scale*(t% tr(i_Rhoc, j) - t% tr(i_Rhoc, j-1))**2  &
                                      + Tc_scale*  (t% tr(i_Tc,   j) - t% tr(i_Tc,   j-1))**2  &
                                      + age_scale* (log10(t% tr(i_age,j)) - log10(t% tr(i_age,j-1)))**2  &
                                      )
         enddo
      endif
      end subroutine distance_along_track

      subroutine dict(input,ncol,cols,output)
      character(len=*), intent(in) :: input
      integer :: ncol
      character(len=col_width), intent(in) :: cols(ncol)
      !character(len=*), intent(in) :: cols(ncol)
      integer, pointer, intent(out) :: output(:)
      integer :: i,ihi,ilo,j
      logical :: have_col(ncol)
      have_col = .false.
      allocate(output(ncol))
      output = 0
      do j=1,ncol
         iloop: do i=1,max_col
            ilo =   1 + main*(i-1) + xtra*(i-1)
            ihi = ilo + main-1
            if(adjustl(adjustr(input(ilo:ihi)))==trim(cols(j)))then
               output(j)=i
               have_col(j) = .true.
               exit iloop
            endif
         enddo iloop
      enddo
      do j=1,ncol
         if(.not.have_col(j)) write(*,*) 'do not have ', trim(cols(j))
      enddo
      end subroutine dict

      subroutine set_eep_interval
      integer :: ierr, i, j, io, my_eep_interval, my_num_eep
      io=alloc_iounit(ierr)
      open(unit=io,file='input.eep',iostat=ierr,status='old')
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
      !to the .eep files; the identifies those columns that are required to
      !locate eeps in the evolutionary tracks
      character(len=file_path) :: history_columns_list
      integer, intent(out) :: ierr
      character(len=col_width) :: col_name
      call process_history_columns(history_columns_list,cols,ierr)
      if(ierr/=0) then
         write(*,*) 'failed in process_history_columns'
         return
      endif
      ncol = size(cols) 
      col_name = 'star_age'; i_age = locate_column(col_name,cols)
      col_name = 'star_mass'; i_mass= locate_column(col_name,cols)
      col_name='log_LH'; i_logLH=locate_column(col_name,cols)
      col_name='log_LHe'; i_logLHe=locate_column(col_name,cols)
      col_name='log_Teff'; i_logTe=locate_column(col_name,cols)
      col_name='log_L'; i_logL=locate_column(col_name,cols)
      col_name='log_g'; i_logg=locate_column(col_name,cols)
      col_name='log_center_T'; i_Tc=locate_column(col_name,cols)
      col_name='log_center_Rho'; i_Rhoc=locate_column(col_name,cols)
      col_name='center_h1'; i_Xc=locate_column(col_name,cols)
      col_name='center_he4'; i_Yc=locate_column(col_name,cols)
      if(old_core_mass_names)then
         col_name='h1_boundary_mass'; i_he_core = locate_column(col_name,cols)
         col_name='he4_boundary_mass'; i_co_core = locate_column(col_name,cols)
      else
         col_name='he_core_mass'; i_he_core = locate_column(col_name,cols)
         col_name='c_core_mass'; i_co_core = locate_column(col_name,cols)
      endif
      end subroutine setup_columns

      end module iso_eep_support
