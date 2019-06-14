program bc_table_convert

  use const_def, only: sp

  implicit none
  type BC
    real(sp), pointer :: data(:)
  end type BC

  type BC_table
     character(len=64) :: photometric_system_string = ''
     character(len=256) :: filename
     character(len=20), allocatable :: labels(:)
     integer :: num_Av, num_Rv, num_filter, num_T, num_g, iling, ilinT
     real(sp) :: Rv, FeH, alphaFe
     real(sp), allocatable :: logT(:), logg(:), Av(:)
     type(BC), allocatable :: bcs(:,:) !Av, filters
  end type BC_table
  type(BC_table) :: new

contains

    subroutine read_one_ascii(t,ierr)
      type(BC_table), intent(inout) :: t
      integer, intent(out) :: ierr
      integer :: i, j, k, pass, ng, nT, num_lines, io
      real(sp) :: Teff, logg, logT
      !integer :: ibcTmin, ibcTmax, ibcgmin, ibcgmax
      real(sp), allocatable :: grid(:,:), BCdata(:) !, bcTmin(:), bcTmax(:), bcgmin(:), bcgmax(:)

      write(*,*) ' opening ', trim(t% filename)
      open(newunit=io,file=trim(t% filename),status='old',action='read',iostat=ierr)
      if(ierr/=0) then
         write(*,*) ' problem opening ascii file, ierr = ', ierr
         return
      endif
      read(io,'(2x,a64)') t% photometric_system_string
      read(io,*) !skip the header
      read(io,'(2x,3i8)') t% num_filter, num_lines, t% num_Av
      read(io,*)

      allocate(t% Av(t% num_Av), grid(2+t% num_filter,num_lines))
      allocate(t% labels(t% num_filter), t% bcs(t% num_filter, t% num_Av))
      allocate(BCdata(t% num_filter))

      do i=1,t% num_Av
         read(io,*) !skip the column numbers
         read(io,'(37x,99a20)') t% labels(1:t% num_Filter)

         do k=1,t% num_filter
            allocate(t% bcs(k,i)% data(num_lines))
         enddo

         do j=1, num_lines
            read(io,'(f8.0,f5.1,4f6.2,99f20.6)') Teff, logg, t% FeH, t% alphaFe, t% Av(i), t% Rv, BCdata
            grid(1,j) = log10(Teff)
            grid(2,j) = logg
            do k=1,t% num_filter
               t% bcs(k,i)% data(j) = BCdata(k)
            enddo
         enddo
      enddo
      close(io)

      write(*,*) ' finished reading . . .'

      !these loops determine a set of unique logT and logg values
      nT=0; ng=0
      do pass=1,2

         if(pass==2)then
            if( ng*nT == num_lines) then !pass
               t% num_T = nT; t% num_g = ng
               allocate(t% logg(t% num_g), t% logT(t% num_T))
            else
               write(0,*) 'nT = ', t% num_T
               write(0,*) 'ng = ', t% num_g
               write(0,*) 'lines = ', num_lines
               stop 'fail!'
            endif
         endif
         logT = -99.; nT = 0
         logg = -99.; ng = 0

         do i = 1, num_lines

            if( logT < grid(1,i) )then
               nT = nT + 1
               logT = grid(1,i)
            endif

            if( logg < grid(2,i) )then
               ng = ng + 1
               logg = grid(2,i)
            endif

            if(pass==2)then
               t% logT(nT) = logT
               t% logg(ng) = logg
            endif

         enddo
      enddo

    end subroutine read_one_ascii

    subroutine write_one_ascii(t,ierr)
      type(BC_table), intent(inout) :: t
      integer, intent(out) :: ierr
      integer :: i, j, k, pass, ng, nT, num_lines, io
      real(sp) :: Teff, logg, logT
      !integer :: ibcTmin, ibcTmax, ibcgmin, ibcgmax
      real(sp), allocatable :: grid(:,:), BCdata(:) !, bcTmin(:), bcTmax(:), bcgmin(:), bcgmax(:)

      write(*,*) ' opening ', trim(t% filename)
      open(newunit=io,file=trim(t% filename),status='unknown',action='write',iostat=ierr)
      if(ierr/=0) then
         write(*,*) ' problem opening ascii file, ierr = ', ierr
         return
      endif
      write(io,'(2x,a64)') t% photometric_system_string
      write(io,'(a1,1x,4a8)') '#', 'filters', 'spectra', 'num Av', 'num Rv'
      write(io,'(a1,1x,4i8)') '#', t% num_filter, t% num_T*t% num_g, t% num_Av, t% num_Rv
      write(io,'(a1)') '#'

      do k=1,t% num_Av

        write(io,'(a1,i7, i5, 3i6, 99(17x,i3))') '#', i,i=1,t% num_filters+5
        write(io,'(a31, a5, 3a6 ,99a20)') '#  Teff logg [Fe/H]   Av    Rv', t% labels
        do i=1,t% num_lines
           write(io,'(f8.0,f5.1,3f6.2,99f20.6)') t% logT(i), t% logg(i), t% FeH, t% Av(k), &
                t% Rv, t% bcs(k)% data(i)
        enddo
        if(k<t% num_Av)then
           write(io,*)
           write(io,*)
        endif
      enddo
      close(io)

      write(*,*) ' finished writing . . .'

    end subroutine write_one_ascii

end program
