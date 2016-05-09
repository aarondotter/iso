module iso_interp_support

  !MESA modules
  use interp_1d_def
  use interp_1d_lib

  !local modules
  use iso_eep_support

  implicit none

  real(dp), parameter :: Z_proto=1.42d-2
  real(dp), parameter :: Y_proto=2.703d-1
  real(dp), parameter :: X_proto=7.155d-1
  real(dp), parameter :: Y_BBN=2.49d-1
  real(dp), parameter :: Z_div_X_proto=log10(Z_proto/X_proto)
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
    integer :: i, ierr
    ierr= 0
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

end module iso_interp_support
