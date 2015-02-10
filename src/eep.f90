module eep

  !MESA modules
  use const_def, only: dp
  use interp_1d_lib, only: interp_4pt_pm

  !local modules
  use iso_eep_support

  implicit none

  logical, parameter :: eep_verbose = .false.

contains

  subroutine secondary_eep(t,s)
    ! takes a regular track (read from history) and creates an
    ! EEP track with a fixed number of regularly spaced points
    type(track), intent(in) :: t
    type(track), intent(out) :: s
    integer :: i, j, k, num_p, num_s

    s% ncol = t% ncol
    !determine total number of EEPs in this track
    num_p = 0
    num_s = 0
    do j=1,primary
       if(t% eep(j) > 0) then
          num_p = num_p + 1
          if(j > 1) then
             num_s = num_s + eep_interval(j-1)
          endif
       endif
    enddo
    if(eep_verbose) write(*,*) ' track, number of p, s eeps = ', i, num_p, num_s
    ! identify the attributes of the new EEP track
    s% filename = t% filename
    s% initial_mass = t% initial_mass
    s% version_number = t% version_number
    s% ntrack = num_p + num_s
    s% neep = num_p
    s% ncol = t% ncol
    allocate(s% cols(s% ncol))
    s% cols = t% cols

    ! allocate and fill new track

    allocate(s% tr(s% ncol,s% ntrack))
    allocate(s% dist(s% ntrack))
    allocate(s% eep(s% neep))
    s% tr = 0d0
    s% dist = 0d0

    k = 1
    ! fill with primary and secondary EEPs
    do j=1,num_p
       if(eep_verbose .and. t% eep(j) > 0) write(*,*) t% tr(i_age,t% eep(j))
       s% tr(:,k) = t% tr(:,t% eep(j))
       s% dist(k) = t% dist(t% eep(j))
       s% eep(j) = k
       !interpolation on distance to fill secondary EEPs
       if(j < num_p) call eep_interpolate(t,j,k,s)
       if(j<num_p) k = k + eep_interval(j) + 1
    enddo
  end subroutine secondary_eep

  subroutine eep_interpolate(t,j,k,s)
    type(track), intent(in) :: t
    integer, intent(in) :: j, k
    type(track), intent(inout) :: s
    real(dp) :: alpha, beta, dist, delta, new_dist, dx, a(3), x(4), y(4)
    real(dp), pointer :: vec(:)
    integer :: i, jj, j0, j1, n, num_eep
    logical, parameter :: linear = .false.

    if(t% eep(j) == t% eep(j+1)) return

    num_eep = eep_interval(j)

    ! determine distance between primary eeps
    dist = t% dist(t% eep(j+1)) - t% dist(t% eep(j))

    ! interval = (total distance) / (number of secondary eeps)
    delta = dist / real(num_eep+1,kind=dp)

    j0 = t% eep(j)
    n = t% ntrack
    do i=1,num_eep
       new_dist = t% dist(t% eep(j)) + delta*real(i,kind=dp)
       allocate(vec(n))
       vec = t% dist
       j0=binary_search(n, vec, j0, new_dist)
       deallocate(vec); nullify(vec)
       j1=j0+1
       if(linear)then
          alpha = (t% dist(t% eep(j+1)) - new_dist)/dist
          beta = 1d0-alpha
          s% tr(:,k+i) = alpha*t% tr(:,j0) + beta*t% tr(:,j1)
       else
          dx = new_dist - t% dist(j0)

          !for interpolation purposes, reset j0 if near the boundaries
          if(j0 > t% ntrack-2) j0 = t% ntrack-2
          if(j0 < 2) j0=2

          do jj=1,t% ncol
             x = t% dist(j0-1:j0+2)
             y = t% tr(jj,j0-1:j0+2)
             call interp_4pt_pm( x, y, a)
             s% tr(jj,k+i) = t% tr(jj,j0) + dx*(a(1) + dx*(a(2) + dx*a(3)))
          enddo
       endif
       s% dist(k+i) = new_dist
    enddo
  end subroutine eep_interpolate


  subroutine primary_eep(t)
    ! sets the locations of the primary EEPs in a track read from a .log
    type(track), intent(inout) :: t
    integer :: ieep
    t% EEP = 0 !initialize
    ieep=1
    t% EEP(ieep) = PreMS_Tc(t,5.0d0,1); if(check(t,ieep)) return; ieep=ieep+1
    t% EEP(ieep) = ZAMS(t,10); if(check(t,ieep)) return; ieep=ieep+1
    t% EEP(ieep) = TAMS(t,1d-1,t% EEP(ieep-1)+2); if(check(t,ieep)) return; ieep=ieep+1
    t% EEP(ieep) = TAMS(t,1d-4,t% EEP(ieep-1)+2); if(check(t,ieep)) return; ieep=ieep+1
    t% EEP(ieep) = RGBTip(t,t% EEP(ieep-1)+2); if(check(t,ieep)) return; ieep=ieep+1
    t% EEP(ieep) = ZAHB(t,t% EEP(ieep-1)+2); if(check(t,ieep)) return; ieep=ieep+1
    t% EEP(ieep) = TAHB(t,1d-1,t% EEP(ieep-1)+2); if(check(t,ieep)) return; ieep=ieep+1
    t% EEP(ieep) = TAHB(t,1d-4,t% EEP(ieep-1)+2); if(check(t,ieep)) return; ieep=ieep+1
    t% EEP(ieep) = TPAGB(t,t% EEP(ieep-1)+2); if(check(t,ieep)) return; ieep=ieep+1
    t% EEP(ieep) = PostAGB(t,t% EEP(ieep-1)+2); if(check(t,ieep)) return; ieep=ieep+1
    t% EEP(ieep) = WDCS(t, t% EEP(ieep-1)+2)
  end subroutine primary_eep

  logical function check(t,i)
    type(track), intent(in) :: t
    integer, intent(in) :: i
    check = (t% EEP(i)==0 .or. t% EEP(i)==t% ntrack)
  end function check

  function PreMS_Tc(t,logTc,guess) result(PreMS)
    type(track), intent(in) :: t
    real(dp), intent(in) :: logTc
    integer, intent(in) :: guess
    integer :: i, my_guess, PreMS
    real(dp) :: my_logTc
    PreMS = 0
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    else
       my_guess = guess
    endif

    ! check for high mass stars, which may
    ! have hotter Tc than the input logTc
    ! in that case try the first value + 0.05
    if(t% tr(i_Tc,1) > logTc) then
       my_logTc = t% tr(i_Tc,1) + 0.05d0
    else
       my_logTc = logTc
    endif

    do i=my_guess,t% ntrack
       if(t% tr(i_Tc,i) > my_logTc) then
          PreMS = i
          return
       endif
    enddo

  end function PreMS_Tc

  !an alternative to constant central temperature is constant age
  function PreMS_age(t,logAge,guess) result(PreMS)
    type(track), intent(in) :: t
    real(dp), intent(in) :: logAge
    integer, intent(in) :: guess
    integer :: i, my_guess, PreMS
    PreMS = 0

    if(t% initial_mass > 3d1)then
       PreMS=1
       return
    endif

    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    else
       my_guess = guess
    endif

    do i=my_guess,t% ntrack
       if(log10(t% tr(i_age,i)) > logAge) then
          PreMS = i
          return
       endif
    enddo
  end function PreMS_age


  integer function ZAMS(t,guess)
    type(track), intent(in) :: t
    real(dp) :: Xmax, Xmin, Xc, LH, Lmin !g_max, 
    real(dp), parameter :: Lfac = 9.99d-1 !Xfac = 9.8d-1
    integer, intent(in) :: guess
    integer :: i, my_guess
    ZAMS = 0
    Xmax = t% tr(i_Xc,1)
    Xmin = Xmax - 1.5d-3
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    else
       my_guess = guess
    endif

    i = my_guess
    Xc = t% tr(i_Xc,i)
    do while(Xc > Xmin .and. i < t% ntrack )
       i = i+1
       Xc = t% tr(i_Xc,i)
    enddo

    ZAMS = i

    LH = pow10(t% tr(i_logLH,i))
    Lmin = Lfac * pow10(t% tr(i_logL,i))
    do while(LH > Lmin)
       i = i-1
       LH = pow10(t% tr(i_logLH,i))
       Lmin = Lfac * pow10(t% tr(i_logL,i))
    enddo

    !set upper limit in case M < 0.1 and no H-burning occurs
    if(t% initial_mass < 0.1d0) ZAMS = min(i+1,t% ntrack)
  end function ZAMS

  integer function TAMS(t,Xmin,guess)
    type(track), intent(in) :: t
    real(dp), intent(in) :: Xmin
    integer, intent(in) :: guess
    integer :: i, my_guess
    real(dp), parameter :: max_age = 2d10
    TAMS = 0
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    elseif(guess == t% ntrack)then
       return
    else
       my_guess = guess
    endif
    do i=my_guess,t% ntrack
       if(t% tr(i_Xc,i) < Xmin) then
          TAMS = i
          return
       endif
    enddo
    ! Xc test fails so consider instead the age for low mass tracks
    ! if the age of the last point > Max_Age, accept it
    if(t% initial_mass <= 0.5d0 .and. t% tr(i_age,t% ntrack) >= max_age) TAMS = t% ntrack
    !if(t% initial_mass <= 0.5d0) TAMS = t% ntrack
  end function TAMS

  integer function RGBTip(t,guess)
    type(track), intent(in) :: t
    integer, intent(in) :: guess
    integer :: i, my_guess
    real(dp) :: Ymin, Lmax
    RGBTip = 0
    if(guess < 1 .or. guess >= t% ntrack) then 
       return
    else
       my_guess = guess
    endif
    if(eep_verbose)then
       write(*,*) '    guess = ', guess
       write(*,*) ' my_guess = ', my_guess
       write(*,*) '   ntrack = ', t% ntrack
    endif
    Ymin = 9d-1
    Lmax = -99d0
    do i=my_guess, t% ntrack
       if(t% tr(i_Yc,i) > Ymin .and. t% tr(i_logL,i) > Lmax)then
          Lmax = t% tr(i_logL,i)
          RGBTip = i
       endif
    enddo
  end function RGBTip

  integer function ZAHB(t,guess)
    type(track), intent(in) :: t
    integer, intent(in) :: guess
    integer :: i, my_guess, my_guess_2
    real(dp) :: Ymin, Tmin, LHemax
    ZAHB = 0
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    elseif(guess == t% ntrack)then
       return
    else
       my_guess = guess
    endif
    Ymin = 9d-1
    LHemax = -99d0
    my_guess_2 = my_guess

    do i=my_guess, t% ntrack
       if(t% tr(i_Yc,i) < 0.98d0)then
          ZAHB=i
          return
       endif
    enddo

    do i=my_guess, t% ntrack
       if(t% tr(i_Yc,i) > Ymin .and. t% tr(i_logLHe,i) > LHemax)then
          LHemax = t% tr(i_logLHe,i)
          my_guess_2 = i
       endif
    enddo
    Tmin = 99d0
    do i=my_guess_2,t% ntrack
       if(t% tr(i_Yc,i) > Ymin .and. t% tr(i_Tc,i) < Tmin) then
          Tmin = t% tr(i_Tc,i)
          ZAHB = i
       endif
    enddo
    if(ZAHB==guess) ZAHB=0
  end function ZAHB

  integer function TAHB(t,Ymin,guess)
    type(track), intent(in) :: t
    real(dp), intent(in) :: Ymin
    integer, intent(in) :: guess
    integer :: i, my_guess
    TAHB = 0
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    elseif(guess == t% ntrack)then
       return
    else
       my_guess = guess
    endif
    do i=my_guess,t% ntrack
       if(t% tr(i_Yc,i) < Ymin)then
          TAHB = i
          return
       endif
    enddo
  end function TAHB

  integer function TPAGB(t,guess)
    type(track), intent(in) :: t
    integer, intent(in) :: guess
    real(dp), parameter :: Ymin = 1d-6
    real(dp), parameter :: HeShellMin = 1d-1
    real(dp) :: Yc, HeShell
    integer :: i, my_guess
    TPAGB = 0
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    elseif(guess == t% ntrack)then
       return
    else
       my_guess = guess
    endif
    do i=my_guess, t% ntrack
       Yc = t% tr(i_Yc,i)
       HeShell = t% tr(i_He_Core,i) - t% tr(i_CO_Core,i)
       if(Yc < Ymin .and. HeShell < HeShellMin)then
          TPAGB = i
          return
       endif
    enddo
  end function TPAGB

  integer function PostAGB(t,guess)
    type(track), intent(in) :: t
    integer, intent(in) :: guess
    real(dp), parameter :: env_mass_min = 8d-1
    real(dp) :: Tc_now, TC_end, env_mass_frac
    integer :: i, my_guess
    PostAGB=0
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    elseif(guess == t% ntrack)then
       return
    else
       my_guess = guess
    endif

    Tc_now = t% tr(i_Tc,my_guess)
    Tc_end = t% tr(i_Tc,t% ntrack)

    if(Tc_now > Tc_end)then
       !has TP-AGB
       do i=my_guess, t% ntrack
          env_mass_frac = t% tr(i_co_core,i) / t% tr(i_mass,i) 
          if(env_mass_frac > env_mass_min)then
             PostAGB = i
             exit
          endif
       enddo
    else
       !has no TP-AGB
       return
    endif
  end function PostAGB

  integer function WDCS(t,guess)
    type(track), intent(in) :: t
    integer, intent(in) :: guess
    real(dp) :: Tc_now, TC_end
    integer :: my_guess
    WDCS=0
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    elseif(guess == t% ntrack)then
       return
    else
       my_guess = guess
    endif

    Tc_now = t% tr(i_Tc,my_guess)
    Tc_end = t% tr(i_Tc,t% ntrack)

    if(Tc_now > Tc_end) WDCS = t% ntrack
  end function WDCS

end module eep
