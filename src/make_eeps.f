      program make_eeps

      !MESA modules
      use utils_lib, only: alloc_iounit, free_iounit

      !local modules
      use iso_eep_support
      use eep
      use phases

      implicit none

      integer :: i, ierr, io, num
      character(len=file_path) :: input_file, history_columns_list
      character(len=file_path), allocatable :: history_files(:)
      type(track), pointer :: t=>NULL(), s=>NULL()
      logical :: do_track_output = .false., do_eep_output = .true.
      logical :: do_phases = .false.

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
         call read_history_file(t)
         write(*,*) trim(t% filename), t% neep, t% version_number
         call primary_eep(t)
         write(*,'(99i8)') t% eep
         allocate(s)
         call secondary_eep(t,s)
         if(do_phases) call set_track_phase(s)
         if(do_track_output) call write_new_track(t)
         if(do_eep_output) call write_new_eep(s)
         deallocate(t,s)
      enddo

      deallocate(cols)

      contains

      subroutine read_input(ierr)
      integer, intent(out) :: ierr
      ierr=0
      call get_command_argument(1,input_file)
      io=alloc_iounit(ierr)
      open(unit=io,file=trim(input_file))
      read(io,*) !skip first line
      read(io,'(a)') data_dir
      read(io,*) !skip comment line
      read(io,'(a)') history_columns_list
      read(io,*) !skip comment line
      read(io,*) num
      allocate(history_files(num))
      do i=1,num
         read(io,'(a)',iostat=ierr) history_files(i)
         if(ierr/=0) exit
      enddo
      close(io)
      call free_iounit(io)
      !set number of secondary EEPs between each primary EEP
      call set_eep_interval 
      !read history file format specs
      io=alloc_iounit(io)
      open(unit=io,file='input.format')
      read(io,*) 
      read(io,*) head
      read(io,*) main
      read(io,*) xtra
      close(io)
      call free_iounit(io)
      !set up columns to be used
      call setup_columns(history_columns_list,ierr)
      end subroutine read_input

      subroutine write_new_track(t)
      type(track) :: t
      !overwrite filename for new output
      t% filename = trim(data_dir) // '/' // trim(t% filename) // '.trk'
      call write_track(t)
      end subroutine write_new_track

      subroutine write_new_EEP(s)
      type(track) :: s
      !overwrite filename for new output
      s% filename = trim(data_dir)// '/' // trim(s% filename) // '.eep'
      call write_track(s) 
      end subroutine write_new_EEP

      end program make_eeps
