program make_track
  !MESA modules
  use const_def, only: dp
  use utils_lib
  use interp_1d_def
  use interp_1d_lib

  !local modules
  use iso_eep_support

  implicit none

contains

  !takes a set of EEP-defined tracks and interpolates a new
  !track for the desired initial mass
  subroutine interpolate_track(a,b,ierr)
    type(track), intent(in) :: a(:)
    type(track), intent(inout) :: b
    integer, intent(out) :: ierr
    real(dp) :: f(3), dm, x(4), y(4)
    real(dp), pointer :: initial_mass(:)=>NULL() !(n)
    integer :: i, j, k, m, mlo, mhi, n

    n = size(a)

    allocate(initial_mass(n))
    initial_mass = a(:)% initial_mass
    m = binary_search(n, initial_mass, 1, b% initial_mass)
    write(*,*) b% initial_mass
    write(*,*) initial_mass(m:m+1)

    mlo = min(max(1,m-1),n-3)
    mhi = max(min(m+2,n),4)
    m = mlo+1

    write(*,*) '   mlo, m, mhi = ', mlo, m, mhi
    write(*,*) initial_mass(mlo:mhi)
    write(*,*) a(mlo:mhi)% neep

    b% neep = minval(a(mlo:mhi)% neep)

    write(*,*) b% neep

    b% ntrack = a(m)% eep(b% neep)

    write(*,*) b% ntrack

    b% ncol = a(1)% ncol
    b% has_phase = a(1)% has_phase
    allocate(b% cols(b% ncol))
    b% cols = a(1)% cols
    b% version_number = a(1)% version_number
    allocate(b% eep(b% neep))
    allocate(b% tr(b% ncol, b% ntrack))
    allocate(b% dist(b% ntrack))
    if(b% has_phase) allocate(b% phase(b% ntrack))
    b% eep = a(m)% eep(1:b% neep)
    if(a(m)% has_phase) b% phase = a(m)% phase(1:b% ntrack)
    b% tr = 0d0
    b% dist = 0d0

    dm = b% initial_mass - initial_mass(m)
    x = initial_mass(mlo:mhi)
    do i=1,b% ntrack
       do j=1,b% ncol
          forall(k=1:4) y(k) = a(mlo-1+k)% tr(j,i)
          call interp_4pt_pm(x, y, f)
          b% tr(j,i) = y(2) + dm*(f(1) + dm*(f(2) + dm*f(3)))
       enddo
    enddo

    write(*,*) ' ierr = ', ierr

    deallocate(initial_mass)
  end subroutine interpolate_track


end program make_track
