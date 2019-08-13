module spotify

  use iso_eep_support
  
  implicit none
  
  real(dp) :: spot_scaling_factor

contains

  subroutine add_spots(s)
    !Chabrier et al. (2007) with minor extension
    !beta = fraction of surface covered by spots
    !gamma = temperature ratio, Tspot/Teff (Chabrier+2007 set gamma=0)
    !Teff and radius are changed while holding luminosity constant
    type(isochrone_set), intent(inout) :: s
    real(dp) :: alfa, beta, gamma, R, Teff, g, gnew, Rnew, Tnew, logR, logTeff, logg
    integer :: i, j, kT=0, kR=0, kg=0
    logical :: spot_function

    spot_function = .true.
    
    !first locate log_Teff and log_L
    do i=1,s% iso(1)% ncol
       if(trim(adjustl(s% iso(1)% cols(i)% name)) == 'log_Teff') then
          kT=i
       else if(trim(adjustl(s% iso(1)% cols(i)% name))== 'log_R') then
          kR=i
       else if(trim(adjustl(s% iso(1)% cols(i)% name))=='log_g') then
          kg=i
       endif
    enddo

    do i=1, s% number_of_isochrones
       do j=1, s% iso(i)% neep
          logTeff=s% iso(i)% data(kT,j)
          logg=s% iso(i)% data(kg,j)
          logR=s% iso(i)% data(kR,j)
          if(logg > 4.0d0)then
             Teff = pow10(logTeff)
             g    = pow10(logg)
             R    = pow10(logR)

             beta =  spot_scaling_factor*berdyugina_beta(Teff) !fraction of spot coverage
             alfa = 1.0d0 - beta !this is the fraction of the unspotted star

             if(spot_function)then
                gamma = berdyugina_gamma(Teff) !this is the ratio Tspot/Teff
             else
                gamma = 0.0d0
             endif
             
             Tnew = pow(alfa*pow(Teff,4.0d0) + beta*pow(Teff*gamma,4.0d0),2.5d-1)
             Rnew = R*pow(Teff/Tnew,2.0d0)
             gnew = g*pow(R/Rnew,2.0d0)
             s% iso(i)% data(kR,j) = log10(Rnew)
             s% iso(i)% data(kT,j) = log10(Tnew)
             s% iso(i)% data(kg,j) = log10(gnew)
          endif
       enddo
    enddo

  contains
    
    function berdyugina_beta(T) result(beta) !spot fraction vs. Teff
      real(dp), intent(in) :: T
      real(dp) :: beta
      real(dp), parameter :: A=-1.77514793d-7, B=1.63313609d-3, C=-3.35621302d0
      if( T > 3098.89d0 .and. T < 6101.11d0)then
         beta = A*T*T + B*T + C
      else
         beta = 0.0d0
      endif
    end function berdyugina_beta

    function berdyugina_gamma(T) result(gamma) !gamma = Tspot/Teff
      real(dp), intent(in) :: T
      real(dp) :: gamma
!!$      real(dp) :: gamma, deltaT
!!$      real(dp), parameter :: alpha=1.0d0
!!$      if (T > 2500.0d0 .and. T < 6000.0d0)then
!!$         if ( T <= 3500.0d0)then
!!$            deltaT = 0.5d0*(T-2500.0d0)
!!$         elseif (T > 3500.0d0 .and. T <= 5500.0d0)then
!!$            deltaT = 0.5d0*alpha*(T-3500.0d0) + 500.0d0
!!$         else
!!$            deltaT = -(1.0d0 + 2.0d0*alpha)*(T-6000.0d0)
!!$         endif
!!$      else
!!$         deltaT = 0.0d0
!!$      endif
!!$      gamma = 1.0d0 - deltaT/T
      gamma = 1.0d0 - 0.5d0*berdyugina_beta(T)
    end function berdyugina_gamma
    
  end subroutine add_spots

end module spotify
