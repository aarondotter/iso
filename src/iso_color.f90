module iso_color

  !MESA modules
  use const_def, only: sp
  use utils_lib, only: alloc_iounit, free_iounit
  use interp_1d_lib

  !local modules
  use iso_eep_support
  use BC_tables
  
  implicit none

  private

  real(sp), parameter :: SolBol=4.74
  real(sp), parameter :: FeH_sol = -1.8729 !log10(0.0134)
  type(BC_table), allocatable :: b(:)


  public iso_color_init, get_mags

contains

  subroutine iso_color_init(bc_table_list,ierr)
    character(len=file_path), intent(in) :: bc_table_list
    integer, intent(out) :: ierr
    call BC_table_init(bc_table_list,b,ierr)
    if(ierr/=0) write(0,*) 'BC tables initialized!'
  end subroutine iso_color_init

  subroutine get_mags(iso,iT,ig,iL,iH,iHe,ierr)
    type(isochrone), intent(inout) :: iso
    integer, intent(in) :: iT, ig, iL, iH, iHe
    integer, intent(out) :: ierr
    integer :: i, n, j, jZ=0
    real(sp), allocatable :: res(:)
    real(sp) :: logT, logg, logL, X, Y, Z, FeH
    iso% nfil = b(1)% num_filter
    allocate(iso% mags(iso% nfil, iso% neep),res(iso% nfil))
    allocate(iso% labels(iso% nfil))
    iso% labels = b(1)% labels
    res = 0.
    iso% mags = 0.0
    n=size(b)
    do i=1,iso% neep
       logT = real(iso% data(iT,i),kind=sp)
       logg = real(iso% data(ig,i),kind=sp)
       logL = real(iso% data(iL,i),kind=sp)
       X    = real(iso% data(iH,i),kind=sp)
       Y    = real(iso% data(iHe,i),kind=sp)
       Z    = 1.0 - X - Y
       FeH  = log10(Z/X) - FeH_sol

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

end module iso_color
