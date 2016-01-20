module phases

  use const_lib, only: dp
  use iso_eep_support

  implicit none


!    t% EEP( 1) = PreMS_Tc(t,5.0d0,1); if(check(t,ieep)) return; ieep=ieep+1
!    t% EEP( 2) = ZAMS(t,10); if(check(t,ieep)) return; ieep=ieep+1
!    t% EEP( 3) = TAMS(t,3.5d-1,t% EEP(ieep-1)+inc); if(check(t,ieep)) return; ieep=ieep+1
!    t% EEP( 4) = TAMS(t,1d-12,t% EEP(ieep-1)+inc); if(check(t,ieep)) return; ieep=ieep+1
!    t% EEP( 5) = RGBTip(t,t% EEP(ieep-1)+inc); if(check(t,ieep)) return; ieep=ieep+1
!    t% EEP( 6) = ZAHB(t,t% EEP(ieep-1)+inc); if(check(t,ieep)) return; ieep=ieep+1
!    t% EEP( 7) = TAHB(t,1d-4,t% EEP(ieep-1)+inc); if(check(t,ieep)) return; ieep=ieep+1
! if(t% star_type == star_low_mass)then
!    t% EEP( 8) = TPAGB(t,t% EEP(ieep-1)+inc); if(check(t,ieep)) return; ieep=ieep+1
!    t% EEP( 9) = PostAGB(t,t% EEP(ieep-1)+inc); if(check(t,ieep)) return; ieep=ieep+1
!    t% EEP(10) = WDCS(t, t% EEP(ieep-1)+inc)
! elseif(t% star_type == star_high_mass)then
!    t% EEP( 8) = CarbonBurn(t,t% EEP(ieep-1)+inc)
! endif

!     FSPS phase information
!     0: main sequence
!     1: sub giant branch
!     2: red giant branch
!     3: red horizontal branch / red clump (core He burning)
!     4: early asymptotic giant branch
!     5: thermally-pulsating AGB (composition sets C-rich vs. O-rich)
!     6: Post-AGB
!     7: blue stragglers
!     8: blue horizontal branch
!     9: Wolf-Rayet (composition sets WC vs. WN)

  real(dp), parameter :: pre_MS = -1d0
  real(dp), parameter :: MS = 0d0
  real(dp), parameter :: SGB = 1d0
  real(dp), parameter :: RGB = 2d0
  real(dp), parameter :: red_HB_clump = 3d0
  real(dp), parameter :: early_AGB = 4d0
  real(dp), parameter :: TP_AGB = 5d0
  real(dp), parameter :: post_AGB = 6d0
  real(dp), parameter :: blue_straggler = 7d0
  real(dp), parameter :: blue_HB = 8d0
  real(dp), parameter :: Wolf_Rayet = 9d0
  real(dp), parameter :: undefined = -9d0

contains

  !translates between isochrone primary EEPs and FSPS phases
  subroutine set_track_phase(t)
    type(track), intent(inout) :: t
    integer :: i, i_p, j, n, next_eep 
    real(dp) :: phase

    n = t% ntrack

    allocate(t% phase(n))

    t% phase = undefined

    do i=1,n
       next_eep = 0
       do j=1,primary-1
          next_eep = next_eep + 1 + eep_interval(j)
          if (i <= next_eep) exit
       enddo

       i_p = j !set i_p to nearest primary EEP <= i

       select case(i_p)
       case(1)
          phase = pre_MS
       case(2:3)
          phase = MS
       case(4)
          phase = RGB
       case(5:6)
          phase = red_HB_clump
       case(7)
          phase = early_AGB
       case(8)
          phase = TP_AGB
       case(9:10)
          phase = post_AGB
       case default
          phase = undefined   
       end select

       if(t% tr(i_logTe,i) >= 4d0 .and. t% tr(i_surfH,i) <= 3d-1 .and. &
            t% star_type == star_high_mass) phase = Wolf_Rayet

       t% phase(i) = phase
    enddo

    t% has_phase = .true.

  end subroutine set_track_phase

end module phases
