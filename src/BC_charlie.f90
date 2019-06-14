module BC_tables

  implicit none

  type BC
     real, pointer :: data(:)
  end type BC

  type BC_table
     character(len=64) :: photometric_system_string = ''
     logical :: is_loaded = .false.
     character(len=256) :: filename
     character(len=20), allocatable :: labels(:)
     integer :: num_Av, num_Rv, num_filter, num_T, num_g, iling, ilinT
     real :: Rv, FeH, alphaFe
     real, allocatable :: logT(:), logg(:), Av(:)
     type(BC), allocatable :: bcs(:,:) !Av, filters
  end type BC_table
  
contains

  subroutine BC_table_init(phot_string,table_list,t,ierr)
    character(len=*), intent(in) :: phot_string,table_list
    type(BC_table), allocatable, intent(inout) :: t(:)
    integer, intent(out) :: ierr
    integer :: i,n,io
    character(len=256) :: line
    n=0
    open(newunit=io,file=trim(table_list),iostat=ierr)
    if(ierr/=0) return
    read(io,*,iostat=ierr) n
    if(ierr/=0) return
    if(n<1)then
       ierr=-1
       write(*,*) 'BC_table_init: no tables to open, n < 1'
       return
    endif
    allocate(t(n))
    i=1
    do while(i<=n)
       read(io,'(a)',iostat=ierr) line
       if(ierr/=0) exit
       if(line=='' .or. line(1:1)=='#' .or. line(1:1)=='!') cycle
       read(line,'(a)') t(i)% filename
       t(i)% filename = trim(t(i)% filename) // '.' // trim(phot_string)
       call load_one_BC(t(i),ierr)
       if(ierr/=0) exit
       if(i>1)then !check to make sure the tables are ordered properly
          if(t(i)% FeH <= t(i-1)% FeH)then
             write(0,*) ' [Fe/H] must be in ascending order in bc_table.list '
             ierr=-1
             return
          endif
       endif
       i=i+1
    enddo
    close(io)
  end subroutine BC_table_init

  subroutine load_bcs(t,ierr)
    type(BC_table), intent(inout) :: t(:)
    integer, intent(out) :: ierr
    integer :: i
    do i=1,size(t)
       call load_one_bc(t(i),ierr)
       if(ierr/=0) exit
    enddo
  end subroutine load_bcs

  subroutine load_one_bc(t,ierr)
    type(BC_table), intent(inout) :: t
    integer, intent(out) :: ierr
    character(len=256) :: binfile
    logical :: binary_exists
    ierr=0
    !if .bin file exists, read it that way
    binfile = trim(t% filename) // '.bin'
    inquire(file=trim(binfile),exist=binary_exists)
    if(binary_exists)then
       call read_one_bin(t,ierr)
    else !read it the old-fashioned way and write a binary
       call read_one_ascii(t,ierr)
       if(ierr==0) call write_one_bin(t,ierr)
    endif
  end subroutine load_one_bc

  subroutine read_one_bin(t,ierr)
    type(BC_table), intent(inout) :: t
    integer, intent(out) :: ierr
    character(len=256) :: binfile
    integer :: i, j, num_lines, io
    ierr=0
    binfile = trim(t% filename) // '.bin'
    open(newunit=io,file=trim(binfile),status='old',form='unformatted',iostat=ierr)
    read(io) t% photometric_system_string
    read(io) t% num_g, t% num_T, t% num_Av, t% num_filter
    allocate(t% Av(t% num_Av))
    allocate(t% bcs(t% num_filter,t% num_Av))
    allocate(t% labels(t% num_filter))
    allocate(t% logg(t% num_g), t% logT(t% num_T))
    read(io) t% labels
    read(io) t% FeH, t% alphaFe, t% Rv
    read(io) t% Av
    read(io) t% logg
    read(io) t% logT
    num_lines = t% num_g * t% num_T
    do i=1,t% num_filter
       do j=1,t% num_Av
          allocate(t% bcs(i,j)% data(num_lines))
          read(io) t% bcs(i,j)% data
       enddo
    enddo
    close(io)
  end subroutine read_one_bin

  subroutine write_one_bin(t,ierr)
    type(BC_table), intent(in) :: t
    integer, intent(out) :: ierr
    character(len=256) :: binfile
    integer :: i, j, io
    ierr=0
    binfile = trim(t% filename) // '.bin'
    open(newunit=io,file=trim(binfile),action='write',form='unformatted',iostat=ierr)
    write(io) t% photometric_system_string
    write(io) t% num_g, t% num_T, t% num_Av, t% num_filter
    write(io) t% labels
    write(io) t% FeH, t% alphaFe, t% Rv
    write(io) t% Av
    write(io) t% logg
    write(io) t% logT
    do i=1,t% num_Filter
       do j=1,t% num_Av
          write(io) t% bcs(i,j)% data
       enddo
    enddo
    close(io)
  end subroutine write_one_bin

  subroutine read_one_ascii(t,ierr)
    type(BC_table), intent(inout) :: t
    integer, intent(out) :: ierr
    integer :: i, j, k, pass, ng, nT, num_lines, io
    real :: Teff, logg, logT
    !integer :: ibcTmin, ibcTmax, ibcgmin, ibcgmax
    real, allocatable :: grid(:,:), BCdata(:) !, bcTmin(:), bcTmax(:), bcgmin(:), bcgmax(:)

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

    t% alphaFe = 0.0
    
    do i=1,t% num_Av
       read(io,*) !skip the column numbers
       read(io,'(31x,99a20)') t% labels(1:t% num_Filter)

       do k=1,t% num_filter
          allocate(t% bcs(k,i)% data(num_lines))
       enddo

       do j=1, num_lines
          read(io,'(f8.0,f5.1,3f6.2,99f20.6)') Teff, logg, t% FeH, t% Av(i), t% Rv, BCdata
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

  function BC_interp_filter_fixed_Av(t,logg,logT,iflt,iAv) result(res)
    type(BC_table), intent(in) :: t
    real, intent(in) :: logT, logg
    integer, intent(in) :: iflt, iAv
    real :: res, res1, res2
    integer :: i, iT, ig, ng, nT
    real :: dT, dg, T1, T2, g1, g2, x(4), y(4), a(3), dx
    iT=0; ig=1

    ng=t% num_g
    nT=t% num_T

    if(logT > t% logT(nT)) then
       iT = nT
    elseif(logT < t% logT(1))then
       iT = 0
    else
       do i = 1, nT-1
          if(logT >= t% logT(i) .and. logT < t% logT(i+1)) iT=i
       enddo
    endif
    if(logg > t% logg(ng)) then
       ig = ng
    elseif(logg < t% logg(1))then
       ig = 0
    else
       do i = 1, ng-1
          if(logg >= t% logg(i) .and. logg < t% logg(i+1)) ig=i
       enddo
    endif

    g1=t% logg(ig)
    g2=t% logg(ig+1)
    dg=g2-g1

    if(ig==0.or.ig==ng) then

       if(ig==0) ig=1

       if(iT==0)then
          res = t% bcs(iflt,iAv)% data(1)
       elseif(iT==nT)then
          res = t% bcs(iflt,iAv)% data(nT*ng)
       elseif(iT==1.or.iT==nT-1)then
          T1=t% logT(iT)
          T2=t% logT(iT+1)
          dT=T2-T1
          res1=t% bcs(iflt,iAv)% data((iT-1)*ng+ig)
          res2=t% bcs(iflt,iAv)% data((iT  )*ng+ig)
          res= ((T2-logT)*res1 + (logT-T1)*res2)/dT
       else
          x=t% logT(iT-1:iT+2)
          dx = logT - x(2)
          y(1) = t% bcs(iflt,iAv)% data((iT-2)*ng + ig)
          y(2) = t% bcs(iflt,iAv)% data((iT-1)*ng + ig)
          y(3) = t% bcs(iflt,iAv)% data((iT  )*ng + ig)
          y(4) = t% bcs(iflt,iAv)% data((iT+1)*ng + ig)
          a = interp_4pt_pm(x,y)
          res = y(2) + dx*(a(1) + dx*(a(2) + dx*a(3)))
       endif

    else

       if(iT==0)then
          res1 = t% bcs(iflt,iAv)% data(ig)
          res2 = t% bcs(iflt,iAV)% data(ig+1)
       elseif(iT==nT)then
          res1 = t% bcs(iflT,iAv)% data((nT-1)*ng+ig)
          res2 = t% bcs(iflT,iAv)% data((nT-1)*ng+ig+1)
       elseif(iT==1.or.iT==nT-1)then
          T1=t% logT(iT)
          T2=t% logT(iT+1)
          dT=T2-T1
          res1 = ((T2-logT)*t% bcs(iflt,iAv)% data((iT-1)*ng+ig) &
               + (logT-T1)*t% bcs(iflt,iAv)% data((iT  )*ng+ig))/dT
          res2 = ((T2-logT)*t% bcs(iflt,iAv)% data((iT-1)*ng+ig+1) &
               + (logT-T1)*t% bcs(iflt,iAv)% data((iT  )*ng+ig+1))/dT
       else
          x = t% logT(iT-1:iT+2)
          dx = logT - x(2)
          !at g1
          y(1) = t% bcs(iflt,iAv)% data((iT-2)*ng + ig)
          y(2) = t% bcs(iflt,iAv)% data((iT-1)*ng + ig)
          y(3) = t% bcs(iflt,iAv)% data((iT  )*ng + ig)
          y(4) = t% bcs(iflt,iAv)% data((iT+1)*ng + ig)
          a = interp_4pt_pm(x,y)
          res1 = y(2) + dx*(a(1) + dx*(a(2) + dx*a(3)))
          !at g2
          y(1) = t% bcs(iflt,iAv)% data((iT-2)*ng + ig+1)
          y(2) = t% bcs(iflt,iAv)% data((iT-1)*ng + ig+1)
          y(3) = t% bcs(iflt,iAv)% data((iT  )*ng + ig+1)
          y(4) = t% bcs(iflt,iAv)% data((iT+1)*ng + ig+1)
          a = interp_4pt_pm(x,y)
          res2 = y(2) + dx*(a(1) + dx*(a(2) + dx*a(3)))
       endif

       res = ((g2-logg)*res1 + (logg-g1)*res2)/dg
    endif
  end function BC_interp_filter_fixed_Av

  function BC_interp_filter(t,logg,logT,Av,FeH,iflt) result(res)
    type(BC_table), intent(in) :: t(:)
    real, intent(in) :: logT, logg, FeH, Av
    integer, intent(in) :: iflt
    real :: res, alfa, beta, y(4), x(4), dx, a(3)
    integer :: i, n, iZ=0

    n=size(t)
    if(n==1 .or. FeH < t(1)% FeH ) then
       res=BC_interp_filter_fixed_metallicity(t(1),logg,logT,Av,iflt)
       return
    elseif(FeH > t(n)% FeH)then
       res=BC_interp_filter_fixed_metallicity(t(n),logg,logT,Av,iflt)
       return       
    endif

    do i=1,n-1
       if(FeH >= t(i)% FeH .and. FeH < t(i+1)% FeH)then
          iZ=i
          exit
       endif
    enddo

    if(iZ==1 .or. iZ >= n-2)then !linear
       alfa=(t(iZ+1)% FeH - FeH)/(t(iZ+1)% FeH-t(iZ)% FeH)
       beta=1.0-alfa
       res = alfa*BC_interp_filter_fixed_metallicity(t(iZ)  ,logg,logT,Av,iflt) &
           + beta*BC_interp_filter_fixed_metallicity(t(iZ+1),logg,logT,Av,iflt)
    else !cubic
       x(1) = t(iZ-1)% FeH
       x(2) = t(iZ  )% FeH
       x(3) = t(iZ+1)% FeH
       x(4) = t(iZ+2)% FeH
       dx = FeH - x(2)
       y(1) = BC_interp_filter_fixed_metallicity(t(iZ-1),logg,logT,Av,iflt)
       y(2) = BC_interp_filter_fixed_metallicity(t(iZ  ),logg,logT,Av,iflt)
       y(3) = BC_interp_filter_fixed_metallicity(t(iZ+1),logg,logT,Av,iflt)
       y(4) = BC_interp_filter_fixed_metallicity(t(iZ+2),logg,logT,Av,iflt)
       a = interp_4pt_pm(x,y)
       res = y(2) + dx*(a(1) + dx*(a(2) + dx*a(3)))
    endif
  end function BC_interp_filter
  
  function BC_interp_filter_fixed_metallicity(t,logg,logT,Av,iflt) result(res)
    type(BC_table), intent(in) :: t
    real, intent(in) :: logT, logg, Av
    integer, intent(in) :: iflt
    real :: res, alfa, beta, y(4), x(4), dx, a(3)
    integer :: i, iAv
    res = 0.
    if(t% num_Av==1 .or. Av < t% Av(1))then
       res=BC_interp_filter_fixed_Av(t,logg,logT,iflt,1)
       return
    elseif(Av > t% Av(t% num_Av))then
       res=BC_interp_filter_fixed_Av(t,logg,logT,iflt,t% num_Av)
       return
    endif
    !general case
    do i=1,t% num_Av-1
       if(Av>=t% Av(i).and. Av < t% Av(i+1)) then
          iAv=i
          exit
       endif
    enddo
    if(iAv==1 .or. iAv>=t% num_Av-2)then !linear
       alfa=(t% Av(iAv+1)-Av)/(t% Av(iAv+1)-t% Av(iAv))
       beta=1.0-alfa
       res  = alfa*BC_interp_filter_fixed_Av(t,logg,logT,iflt,iAv  ) &
            + beta*BC_interp_filter_fixed_Av(t,logg,logT,iflt,iAv+1)
    else !cubic
       x=t% Av(iAv-1:iAv+2)
       dx=Av-x(2)
       y(1) = BC_interp_filter_fixed_Av(t,logg,logT,iflt,iAv-1)
       y(2) = BC_interp_filter_fixed_Av(t,logg,logT,iflt,iAv  )
       y(3) = BC_interp_filter_fixed_Av(t,logg,logT,iflt,iAv+1)
       y(4) = BC_interp_filter_fixed_Av(t,logg,logT,iflt,iAv+2)
       a = interp_4pt_pm(x,y)
       res = y(2) + dx*(a(1) + dx*(a(2) + dx*a(3)))
    endif
  end function BC_interp_filter_fixed_metallicity
  
  subroutine unload_bcs(t)
    type(BC_table), intent(inout) :: t(:)
    integer :: i
    do i=1,size(t)
       call unload_one_bc(t(i))
    enddo
  end subroutine unload_bcs

  subroutine unload_one_bc(t)
    type(BC_table), intent(inout) :: t
    deallocate(t% logT, t% logg, t% Av, t% bcs)
    t% is_loaded = .false.
  end subroutine unload_one_bc
  
  function interp_4pt_pm(x, y) result(a)
    ! returns coefficients for monotonic cubic interpolation from x(2) to x(3)
    real, intent(in)    :: x(4)    ! junction points, strictly monotonic
    real, intent(in)    :: y(4)    ! data values at x's
    real :: a(3), h1, h2, h3, s1, s2, s3, p2, p3, as2, ss2, yp2, yp3
    ! for x(2) <= x <= x(3) and dx = x-x(2), 
    ! y(x) = y(2) + dx*(a(1) + dx*(a(2) + dx*a(3)))
    h1 = x(2)-x(1); h2 = x(3)-x(2); h3 = x(4)-x(3)
    s1 = (y(2)-y(1))/h1; s2 = (y(3)-y(2))/h2; s3 = (y(4)-y(3))/h3
    p2 = (s1*h2+s2*h1)/(h1+h2); p3 = (s2*h3+s3*h2)/(h2+h3)
    as2 = abs(s2); ss2 = sign(1.0, s2)
    yp2 = (sign(1.0, s1)+ss2)*min(abs(s1), as2, 0.5*abs(p2))
    yp3 = (ss2+sign(1.0, s3))*min(as2, abs(s3), 0.5*abs(p3))
    a(1) = yp2
    a(2) = (3.0*s2-2.0*yp2-yp3)/h2
    a(3) = (yp2+yp3-2.0*s2)/(h2*h2)
  end function interp_4pt_pm
  
end module BC_tables


program test
  use BC_tables
  implicit none
  character(len=256) :: phot_string, table_list
  integer :: ierr
  type(BC_table), allocatable :: bcs(:)
  integer, parameter :: nflt=57 !for WFC3
  integer, parameter :: i275=4
  integer, parameter :: i336=7
  integer, parameter :: i438=15
  integer, parameter :: i555=23
  integer, parameter :: i606=25
  integer, parameter :: i814=39
  real :: logT, logg, FeH, Av, Mbol, mag

  Mbol=4.74
  phot_string='HST_WFC3'
  table_list='/home/dotter/science/iso/bc_table.list'

  call BC_table_init(phot_string, table_list, bcs, ierr)

  logT=log10(5778.0)
  logg=4.45
  FeH = 0.0
  Av = 0.0

  !example of one call
  mag = BC_interp_filter(bcs,logg,logT,Av,FeH,i275)

  !various calls
  write(*,*) trim(bcs(1)% photometric_system_string)
  write(*,'(1x,a20,a13)') 'filter', 'magnitude'
  write(*,*) bcs(1)% labels(i275), Mbol - BC_interp_filter(bcs,logg,logT,Av,FeH,i275)
  write(*,*) bcs(1)% labels(i336), Mbol - BC_interp_filter(bcs,logg,logT,Av,FeH,i336)
  write(*,*) bcs(1)% labels(i438), Mbol - BC_interp_filter(bcs,logg,logT,Av,FeH,i438)
  write(*,*) bcs(1)% labels(i555), Mbol - BC_interp_filter(bcs,logg,logT,Av,FeH,i555)
  write(*,*) bcs(1)% labels(i606), Mbol - BC_interp_filter(bcs,logg,logT,Av,FeH,i606)
  write(*,*) bcs(1)% labels(i814), Mbol - BC_interp_filter(bcs,logg,logT,Av,FeH,i814)

end program test

