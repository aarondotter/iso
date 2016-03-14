module iso_color

  !MESA modules
  use interp_1d_lib

  !local modules
  use iso_eep_support
  use BC_tables

  implicit none

  private

  real(sp), parameter :: SolBol=4.74
  real(sp), parameter :: Z_sol = 0.0134
  real(sp), parameter :: log_Z_sol = log10(Z_sol)
  real(sp) :: BC_Fe_div_H
  type(BC_table), allocatable :: b(:), c(:)
  logical :: BC_do_Cstars = .false., do_fixed_Z = .false.
  public iso_color_init, write_cmds_to_file

contains

  subroutine iso_color_init(bc_table_list,do_Cstars,cstar_table_list,set_fixed_Fe_div_H,Fe_div_H,ierr)
    character(len=file_path), intent(in) :: bc_table_list, cstar_table_list
    logical, intent(in) :: do_Cstars
    logical, intent(in) :: set_fixed_Fe_div_H
    real(sp), intent(in) :: Fe_div_H
    integer, intent(out) :: ierr
    BC_do_Cstars = do_Cstars
    call BC_table_init(bc_table_list,b,ierr)
    if(set_fixed_Fe_div_H)then
       BC_Fe_div_H = Fe_div_H
       do_fixed_Z = .true.
    endif
    if(BC_do_Cstars) call BC_table_init(cstar_table_list,c,ierr)
    if(ierr/=0) write(0,*) 'BC tables initialized!'
  end subroutine iso_color_init

  subroutine get_mags(iso,ierr)
    type(isochrone), intent(inout) :: iso
    integer, intent(out) :: ierr
    integer :: i, nb, nc, j, n, jZ
    real(sp), allocatable :: res(:)
    real(sp) :: logT, logg, logL, X, Y, Z, FeH
    real(sp) :: c_min_logT, c_max_logT, c_min_logg, c_max_logg
    real(dp) :: C_div_O
    logical :: Cstar_ok
    integer :: iT, ig, iL, iH, iHe, iC, iO, isurfZ, iCdivO
    !a few initializations
    c_min_logT=0; c_max_logT=0; c_min_logg=0; c_max_logg=0
    iT=0; ig=0; iL=0; iH=0; iHe=0; iC=0; iO=0; jZ=0; isurfZ=0; iCdivO=0

    !locate columns
    do i = 1, iso% ncol
       if(trim(adjustl(iso% cols(i)% name)) == 'log_Teff') then
          iT=i
       else if(trim(adjustl(iso% cols(i)% name)) == 'log_g') then
          ig=i
       else if(trim(adjustl(iso% cols(i)% name))== 'log_L') then
          iL=i
       else if(trim(adjustl(iso% cols(i)% name))=='surface_h1')then
          iH=i
       else if(trim(adjustl(iso% cols(i)% name))=='surface_he4')then
          iHe=i   
       else if(trim(adjustl(iso% cols(i)% name))=='surface_c12')then
          iC=i
       else if(trim(adjustl(iso% cols(i)% name))=='surface_o16')then
          iO=i
       else if(trim(adjustl(iso% cols(i)% name))=='log_surf_z')then
          isurfZ=i
       else if(trim(adjustl(iso% cols(i)% name))=='surf_num_c12_div_num_o16')then
          iCdivO=i
       endif
    enddo

    if(BC_do_Cstars)then
       c_min_logT = minval(c(1)% logT); c_min_logg = minval(c(1)% logg)
       c_max_logT = maxval(c(1)% logT); c_max_logg = maxval(c(1)% logg)
    endif

    iso% nfil = b(1)% num_filter
    allocate(iso% mags(iso% nfil, iso% neep),res(iso% nfil))
    allocate(iso% labels(iso% nfil))
    iso% labels = b(1)% labels
    res = 0.
    iso% mags = 0.
    nb=size(b)
    if(BC_do_Cstars)then
       nc=size(c)
    else
       nc=0
    endif


    do i=1,iso% neep
       logT = real(iso% data(iT,i),kind=sp)
       logg = real(iso% data(ig,i),kind=sp)
       logL = real(iso% data(iL,i),kind=sp)
       if(do_fixed_Z)then
          FeH = BC_Fe_div_H
       else
          if(isurfZ>0)then
             FeH = real(iso% data(isurfZ,i),kind=sp) - log_Z_sol
          else
             X    = real(iso% data(iH,i),kind=sp)
             Y    = real(iso% data(iHe,i),kind=sp)
             Z    = 1.0 - X - Y 
             FeH = log10(Z) - log_Z_sol
          endif
       endif
       if(iCdivO>0)then
          C_div_O  =log10(iso% data(iCdivO,i))
       else
          C_div_O = log10((16d0/12d0) * iso% data(iC,i) / iso% data(iO,i))
       endif

       Cstar_OK = (BC_do_Cstars) .and. (C_div_O > 0d0) .and. (logT > c_min_logT) .and. &
            (logT < c_max_logT) .and. (logg > c_min_logg) .and. (logg < c_max_logg)

       if( Cstar_OK )then !use Cstar grid
          n=nc
          !interpolation is either linear, quadratic, or cubic
          if(FeH < c(1)% FeH .or. n==1) then
             call BC_interp_filters(c(1),logg,logT,iso% Av, res,ierr)
          else if(FeH >= c(n)% FeH) then
             call BC_interp_filters(c(n),logg,logT,iso% Av, res,ierr)

          else if(n==2) then
             call linear(c(1), c(2), logg, logT, iso% Av, FeH, res, ierr)

          else if(n==3)then
             call quadratic(c(1), c(2), c(3), logg, logT, iso% Av, FeH, res, ierr)

          else 
             do j=1,n-1
                if(FeH >= c(j)% FeH .and. FeH < c(j+1)% FeH) jZ=j
             enddo

             if(jZ==1) then
                call quadratic(c(1), c(2), c(3), logg, logT, iso% Av, FeH, res, ierr)
             else if(jZ==n-1) then
                call quadratic(c(n-2), c(n-1), c(n), logg, logT, iso% Av, FeH, res, ierr)
             else
                call cubic(c(jZ-1),c(jZ),c(jZ+1),c(jZ+2), logg, logT, iso% Av, FeH, res, ierr)
             endif

          endif

       else !use standard grid

          n=nb
          !interpolation is either linear, quadratic, or cubic
          if(FeH < b(1)% FeH .or. n==1) then
             call BC_interp_filters(b(1),logg,logT,iso% Av, res,ierr)
          else if(FeH >= b(n)% FeH) then
             call BC_interp_filters(b(n),logg,logT,iso% Av, res,ierr)

          else if(n==2) then
             call linear(b(1), b(2), logg, logT, iso% Av, FeH, res, ierr)

          else if(n==3)then
             call quadratic(b(1), b(2), b(3), logg, logT, iso% Av, FeH, res, ierr)

          else 
             do j=1,n-1
                if(FeH >= b(j)% FeH .and. FeH < b(j+1)% FeH) jZ=j
             enddo

             if(jZ==1) then
                call quadratic(b(1), b(2), b(3), logg, logT, iso% Av, FeH, res, ierr)
             else if(jZ==n-1) then
                call quadratic(b(n-2), b(n-1), b(n), logg, logT, iso% Av, FeH, res, ierr)
             else
                call cubic(b(jZ-1),b(jZ),b(jZ+1),b(jZ+2), logg, logT, iso% Av, FeH, res, ierr)
             endif

          endif

       endif
       
       if(has_bad_real(iso% nfil, res))then
          write(0,*) ' logT =', logT
          write(0,*) ' logg =', logg
          write(0,*) ' logL =', logL
          write(0,*) '  FeH =', FeH
          write(0,*) '    X =', X
          write(0,*) '    Y =', Y
          write(0,*) '    Z =', Z
          stop 'get_mags: bad value!'
       endif

       iso% mags(:,i) = SolBol - 2.5*logL - res
    enddo


  contains

    subroutine linear(t1,t2,logg,logT,Av,FeH,res,ierr)
      type(bc_table), intent(inout) :: t1, t2
      real(sp), intent(in) :: logg, logT, Av, FeH
      real(sp), intent(out) :: res(:) !iso% nfil
      integer, intent(out) :: ierr
      real(sp), allocatable :: res_tmp(:,:), slope(:), intercept(:)
      integer :: n
      n=size(res)
      allocate(res_tmp(n,2),slope(n),intercept(n))
      call BC_interp_filters(t1, logg, logT, Av, res_tmp(:,1), ierr)
      call BC_interp_filters(t2, logg, logT, Av, res_tmp(:,2), ierr)
      !linear interpolation
      slope = (res_tmp(:,2)-res_tmp(:,1))/(b(2)% FeH-b(1)% FeH)
      intercept = (res_tmp(:,1)*b(2)% FeH-res_tmp(:,2)*b(1)% FeH)/(b(2)% FeH-b(1)% FeH)
      res = slope*FeH + intercept
    end subroutine linear

    subroutine quadratic(t1,t2,t3,logg,logT,Av,FeH,res,ierr)
      type(bc_table), intent(inout) :: t1,t2,t3
      real(sp), intent(in) :: logg, logT, Av, FeH
      real(sp), intent(out) :: res(:) !iso% nfil
      integer, intent(out) :: ierr
      real(sp), allocatable :: res_tmp(:,:)
      real(sp) :: x12, x23, x0
      integer :: k, n
      character(len=file_path) :: str
      n=size(res)
      allocate(res_tmp(n,3))
      x12 = t2% FeH - t1% FeH
      x23 = t3% FeH - t2% FeH
      x0  =     FeH - t2% FeH
      call BC_interp_filters(t1, logg, logT, Av, res_tmp(:,1), ierr)
      call BC_interp_filters(t2, logg, logT, Av, res_tmp(:,2), ierr)
      call BC_interp_filters(t3, logg, logT, Av, res_tmp(:,3), ierr)
      do k = 1, iso% nfil
         call interp_3_to_1_sg(x12,x23,x0,res_tmp(k,1),res_tmp(k,2),res_tmp(k,3),res(k),str,ierr)
         if(ierr/=0) then
            write(*,*) trim(str)
            return
         endif
      enddo
    end subroutine quadratic

    subroutine cubic(t1,t2,t3,t4,logg,logT,Av,FeH,res,ierr)
      type(bc_table), intent(inout) :: t1,t2,t3,t4
      real(sp), intent(in) :: logg, logT, Av, FeH
      real(sp), intent(out) :: res(:) !iso% nfil
      integer, intent(out) :: ierr
      real(sp), allocatable :: res_tmp(:,:)
      real(sp) :: x12, x23, x34, x0
      integer :: k, n
      character(len=file_path) :: str
      n=size(res)
      allocate(res_tmp(n,4))
      x12 = t2% FeH - t1% FeH
      x23 = t3% FeH - t2% FeH
      x34 = t4% FeH - t3% FeH
      x0  =     FeH - t2% FeH
      call BC_interp_filters(t1, logg, logT, Av, res_tmp(:,1), ierr)
      call BC_interp_filters(t2, logg, logT, Av, res_tmp(:,2), ierr)
      call BC_interp_filters(t3, logg, logT, Av, res_tmp(:,3), ierr)
      call BC_interp_filters(t4, logg, logT, Av, res_tmp(:,4), ierr)
      do k = 1, iso% nfil
         call interp_4_to_1_sg(x12,x23,x34,x0,res_tmp(k,1),res_tmp(k,2),res_tmp(k,3),res_tmp(k,4),res(k),str,ierr)
         if(ierr/=0) then
            write(*,*) trim(str)
            return
         endif
      enddo
    end subroutine cubic

  end subroutine get_mags


  subroutine write_cmds_to_file(set)
    type(isochrone_set), intent(inout) :: set
    character(len=256) :: output
    integer :: i, io, n, ierr
    if(set% cmd_suffix/='') output = trim(set% filename) // '.' // trim(set% cmd_suffix)
    n=size(set% iso)
    write(0,*) ' cmd output file = ', trim(output)
    io = alloc_iounit(ierr)
    if(ierr/=0) return
    open(io,file=trim(output),action='write',status='unknown',iostat=ierr)
    if(ierr/=0) return
    write(io,'(a25,i5)')    '# number of isochrones = ', n
    write(io,'(a25,i5)')    '# MESA version number  = ', set% version_number
    write(io,'(a25,2f6.3)') '# CCM89 extinction: Av = ', set% Av
    do i=1,n
       set% iso(i)% Av = set% Av
       call write_cmd_to_file(io, set% iso(i))
       if(i<n) write(io,*)
       if(i<n) write(io,*)
    enddo
    close(io)
    call free_iounit(io)
  end subroutine write_cmds_to_file
  

  subroutine write_cmd_to_file(io,iso)
    integer, intent(in) :: io
    type(isochrone), intent(inout) :: iso
    integer :: i, iT, ig, iL, ierr, iM
    ierr=0; iT=0; ig=0; iL=0; iM=0

    do i=1, iso% ncol
       if(trim(adjustl(iso% cols(i)% name)) == 'log_Teff') then
          iT=i
       else if(trim(adjustl(iso% cols(i)% name)) == 'log_g') then
          ig=i
       else if(trim(adjustl(iso% cols(i)% name)) == 'log_L') then
          iL=i
       else if(trim(adjustl(iso% cols(i)% name)) == 'initial_mass')then
          iM=i
       endif
    enddo

    call get_mags(iso,ierr)
    if(ierr/=0) write(0,*) ' problem in get_mags '

    write(io,'(a25,2i5)') '# number of EEPs, cols = ', iso% neep, iso% nfil + 6
    write(io,'(a1,i4,5i32,299(17x,i3))') '#    ', (i,i=1,iso% nfil+6)

    if(iso% age_scale==age_scale_linear)then
       write(io,'(a5,5a32,299a20)') '# EEP', 'isochrone_age_yr', 'initial_mass', 'log_Teff', &
            'log_g', 'log_L', adjustr(iso% labels)
    elseif(iso% age_scale==age_scale_log10)then
       write(io,'(a5,5a32,299a20)') '# EEP', 'log10_isochrone_age_yr', 'initial_mass', 'log_Teff', &
            'log_g', 'log_L', adjustr(iso% labels)
    endif

    do i = 1,iso% neep
       write(io,'(i5,5(1pes32.16e3),299(0pf20.6))') iso% eep(i), iso% age, iso% data(iM,i), &
            iso% data(iT,i), iso% data(ig,i), iso% data(iL,i), iso% mags(:,i)
    enddo
  end subroutine write_cmd_to_file


end module iso_color
