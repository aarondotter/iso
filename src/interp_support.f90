module interp_support

  !MESA modules
  use interp_1d_def
  use interp_1d_lib

  !local modules
  use iso_eep_support

  implicit none

  logical :: force_linear=.false.

  real(dp), parameter :: Z_proto=1.42d-2
  real(dp), parameter :: Y_proto=2.703d-1
  real(dp), parameter :: X_proto=7.155d-1
  real(dp), parameter :: Y_BBN=2.49d-1
  real(dp), parameter :: Z_div_X_proto=-1.7023212937127388 !log10(Z_proto/X_proto)
  real(dp), parameter :: dY_dZ = (Y_proto-Y_BBN)/Z_proto

contains

  function get_initial_Y_and_Z(Fe_div_H) result(Xi)
    real(dp), intent(in) :: Fe_div_H
    real(dp) :: Z_div_X, top, bottom, Xi(3)
    Z_div_X = pow10(Fe_div_H + Z_div_X_proto)
    top = 1d0 - Y_BBN 
    bottom = 1d0 + Z_div_X*(1d0+dY_dZ)
    Xi(1) = top/bottom           !X
    Xi(3) = Z_div_X*Xi(1)        !Z
    Xi(2) = 1d0 - Xi(1) - Xi(3)  !Y
  end function get_initial_Y_and_Z

  subroutine set_initial_Y_and_Z_for_eep(t)
    type(track), intent(inout) :: t
    real(dp) :: Xi(3)
    Xi = get_initial_Y_and_Z(t% Fe_div_H)
    t% initial_Y = Xi(2)
    t% initial_Z = Xi(3)
  end subroutine set_initial_Y_and_Z_for_eep

  subroutine set_initial_Y_and_Z_for_iso(t)
    type(isochrone_set), intent(inout) :: t
    real(dp) :: Xi(3)
    Xi = get_initial_Y_and_Z(t% Fe_div_H)
    t% initial_Y = Xi(2)
    t% initial_Z = Xi(3)
  end subroutine set_initial_Y_and_Z_for_iso

  subroutine consistency_check(s,ierr)
    type(isochrone_set), intent(in) :: s(:)
    integer, intent(out) :: ierr
    integer :: n, i

    ierr = 0

    n=size(s)

    if(any(s(:)% number_of_isochrones /= s(1)% number_of_isochrones))then
       write(0,*) ' iso_interp_met : failed consistency check - inconsistent number of isochrones '
       ierr=-1
       do i=1,n
          write(0,*) trim(s(i)% filename), s(i)% MESA_revision_number
       enddo
    endif

    if(any(s(:)% MESA_revision_number /= s(1)% MESA_revision_number))then
       write(0,*) ' iso_interp_met : failed consistency check - inconsistent version number '
       ierr=-1
       do i=1,n
          write(0,*) trim(s(i)% filename), s(i)% number_of_isochrones
       enddo
    endif

    !more checks here as needed

  end subroutine consistency_check

  function linear(ncol,Z,newZ,x,y) result(res)
    integer, intent(in) :: ncol
    real(dp), intent(in) :: Z(2), newZ, x(ncol), y(ncol)
    real(dp) :: res(ncol), alfa, beta
    alfa = (Z(2)-newZ)/(Z(2)-Z(1))
    beta = 1d0 - alfa
    res = alfa*x + beta*y
  end function linear

  function cubic_pm(ncol,Z,newZ,v,w,x,y) result(res)
    integer, intent(in) :: ncol
    real(dp), intent(in) :: Z(4), newZ, v(ncol), w(ncol), x(ncol), y(ncol)
    real(dp) :: res(ncol), Q(4), a(3), dZ
    integer :: i
    dZ = newZ - Z(2)
    do i=1,ncol
       Q=[v(i),w(i),x(i),y(i)]
       call interp_4pt_pm(Z, Q, a)
       res(i) = Q(2) + dZ*(a(1) + dZ*(a(2) + dZ*a(3)))
    enddo
  end function cubic_pm

  function quadratic(ncol,Z,newZ,w,x,y) result(res)
    integer, intent(in) :: ncol
    real(dp), intent(in) :: Z(3), newZ, w(ncol), x(ncol), y(ncol)
    real(dp) :: res(ncol), x12, x23, x00
    integer :: i, ierr
    character(len=file_path) :: str
    ierr= 0
    x00 = newZ - Z(1) 
    x12 = Z(2) - Z(1)
    x23 = Z(3) - Z(2)
    do i=1,ncol
       call interp_3_to_1(x12, x23, x00, w(i), x(i), y(i), res(i), str, ierr)
    enddo
  end function quadratic

  function cubic(ncol,Z,newZ,v,w,x,y) result(res)
    integer, intent(in) :: ncol
    real(dp), intent(in) :: Z(4), newZ, v(ncol), w(ncol), x(ncol), y(ncol)
    real(dp) :: res(ncol), x12, x23, x34, x00
    integer :: i, ierr
    character(len=file_path) :: str
    ierr= 0
    x00 = newZ - Z(1)
    x12 = Z(2) - Z(1)
    x23 = Z(3) - Z(2)
    x34 = Z(4) - Z(3)
    do i=1,ncol
       call interp_4_to_1(x12, x23, x34, x00, v(i), w(i), x(i), y(i), res(i), str, ierr)
    enddo
  end function cubic

  function cubic_generic(ncol,Z,newZ,v,w,x,y) result(res)
    integer, intent(in) :: ncol
    real(dp), intent(in) :: Z(4), newZ, v(ncol), w(ncol), x(ncol), y(ncol)
    real(dp) :: res(ncol), Q(4), x_new(1), y_new(1)
    integer, parameter :: nwork = max(pm_work_size, mp_work_size)        
    real(dp), target :: work_ary(4*nwork)
    real(dp), pointer :: work(:)
    integer :: i, ierr
    character(len=file_path) :: str
    !can use: interp_m3a, interp_m3b, interp_m3q
    work => work_ary
    x_new(1) = newZ
    ierr= 0
    do i=1,ncol
       Q=[v(i),w(i),x(i),y(i)]
       call interpolate_vector(4, Z, 1, x_new, Q, y_new, interp_m3q, nwork, work, str, ierr)
       res(i) = y_new(1)
    enddo
  end function cubic_generic


  subroutine interpolate_mass_Z(a1,a2,b,ierr)
    type(track), intent(in) :: a1(:), a2(:)
    type(track), intent(inout) :: b
    integer, intent(out) :: ierr
    type(track) :: t1, t2
    integer :: i
    real(dp) :: alfa, beta
    ierr=0

    t1% initial_mass = b% initial_mass
    t2% initial_mass = b% initial_mass

    call interpolate_mass(a1,t1,ierr)
    if(ierr/=0) return
    call interpolate_mass(a2,t2,ierr)
    if(ierr/=0) return

    alfa = (t2% Fe_div_H - b% Fe_div_H)/(t2% Fe_div_H - t1% Fe_div_H)
    beta = 1d0 - alfa

    b% star_type = t1% star_type
    b% alpha_div_Fe = t1% alpha_div_Fe
    b% v_div_vcrit = t1% v_div_vcrit
    b% version_string = t1% version_string
    b% MESA_revision_number = t1% MESA_revision_number
    b% initial_mass = t1% initial_mass
    b% has_phase = t1% has_phase
    b% ncol = t1% ncol
    allocate(b% cols(b% ncol))
    b% cols = t1% cols
    call set_initial_Y_and_Z_for_eep(b)

    b% neep = min(t1% neep, t2% neep)
    allocate(b% eep(b% neep))
    b% eep = t1% eep(1:b% neep)
    b% ntrack = min(t1% ntrack, t2% ntrack)
    allocate(b% tr(b% ncol, b% ntrack))

    do i=1,b% ntrack
       b% tr(:,i) = alfa*t1% tr(:,i) + beta*t2% tr(:,i)      
    enddo
    if(b% has_phase)then
       allocate(b% phase(b% ntrack))
       b% phase = t1% phase(1:b% ntrack)
    endif
  end subroutine interpolate_mass_Z

  !takes a set of EEP-defined tracks and interpolates a new
  !track for the desired initial mass
  subroutine interpolate_mass(a,b,ierr)
    type(track), intent(in) :: a(:)
    type(track), intent(inout) :: b
    integer, intent(out) :: ierr
    real(dp) :: f(3), dx, x(4), y(4), alfa, beta
    real(dp), pointer :: initial_mass(:)=>NULL() !(n)
    integer :: i, j, k, m, mlo, mhi, n, iage, imass

    ierr=0
    dx=0d0; alfa=0d0; beta=0d0; x=0d0; y=0d0; iage=0; imass=0
    n = size(a)

    allocate(initial_mass(n))
    initial_mass = a(:)% initial_mass
    m = binary_search(n, initial_mass, 1, b% initial_mass)
    if(m<1 .or. m>=n)then
       ierr=-1
       return
    endif

    if(force_linear)then
       mlo=max(m,1)
       mhi=min(mlo+1,n)
       if(mhi==mlo) mlo=mhi-1
       if(a(mlo)% neep < a(mhi)% neep)then
          k=mlo
       else
          k=mhi
       endif
    else
       mlo = min(max(1,m-1),n-3)
       mhi = max(min(m+2,n),4)
       k = minloc(a(mlo:mhi)% neep,dim=1) + mlo - 1
    endif

    b% neep = a(k)% neep
    b% ntrack = a(k)% ntrack
    b% star_type = a(m)% star_type
    b% version_string = a(m)% version_string
    b% initial_Z = a(m)% initial_Z
    b% initial_Y = a(m)% initial_Y
    b% Fe_div_H = a(m)% Fe_div_H
    b% alpha_div_Fe = a(m)% alpha_div_Fe
    b% v_div_vcrit = a(m)% v_div_vcrit
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

    if(force_linear)then
       alfa = (b% initial_mass - a(mlo)% initial_mass)/(a(mhi)% initial_mass - a(mlo)% initial_mass)
       beta = 1d0 - alfa
    else
       x = initial_mass(mlo:mhi)
       dx = b% initial_mass - x(2)
    endif

    do i=1,b% ntrack
       do j=1,b% ncol
          if(force_linear)then
             b% tr(j,i) = alfa*a(mhi)% tr(j,i) + beta*a(mlo)% tr(j,i)
          else
             do k=1,4
                y(k) = a(mlo-1+k)% tr(j,i)
             enddo
             call interp_4pt_pm(x, y, f)
             b% tr(j,i) = y(2) + dx*(f(1) + dx*(f(2) + dx*f(3)))
          endif
       enddo
    enddo

    do i=1,b% ncol
       if(adjustl(adjustr(b% cols(i)% name)) == 'star_age') iage=i
       if(adjustl(adjustr(b% cols(i)% name)) == 'star_mass') imass=i
    enddo

    !fix any problems
    call PAV(b% tr(iage,:))
    do i=2,b% ntrack
       b% tr(imass,i) = min(b% tr(imass,i), b% tr(imass,i-1))
    enddo

    deallocate(initial_mass)
  end subroutine interpolate_mass

end module interp_support
