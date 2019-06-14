module spotify

  use iso_eep_support
  
  implicit none
  
  real(dp) :: spot_beta, spot_gamma, Teff_spot_full_on, Teff_spot_full_off

contains

  subroutine add_spots(s,beta0,gamma)
    !Chabrier et al. (2007) with minor extension
    !beta = fraction of surface covered by spots
    !gamma = temperature ratio, Tspot/Teff (Chabrier+2007 set gamma=0)
    !Teff and radius are changed while holding luminosity constant
    type(isochrone_set), intent(inout) :: s
    real(dp), intent(in) :: beta0, gamma
    real(dp) :: a, Tfac, alfa, beta, R, Teff, g, log_Teff
    real(dp) :: log_Teff_spot_full_on, log_Teff_spot_full_off, gnew, Rnew, Tnew
    integer :: i, j, kT=0, kR=0, kg=0

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

    log_Teff_spot_full_on = log10(Teff_spot_full_on)
    log_Teff_spot_full_off = log10(Teff_spot_full_off)

    do i=1, s% number_of_isochrones
       do j=1, s% iso(i)% neep
          log_Teff = s% iso(i)% data(kT,j)
          if(log_Teff > log_Teff_spot_full_off) then
             a = 0.0d0
          elseif(log_Teff < log_Teff_spot_full_on) then
             a = 1.0d0
          else
             a = (log_Teff_spot_full_off - log_Teff)/(log_Teff_spot_full_off-log_Teff_spot_full_on)
          endif

          Tfac = 0.5d0*(1d0 - cos(pi*a))

          !combine all multiplicative factors into one
          beta = Tfac*beta0 !this is the fraction of the spots
          alfa = 1d0 - beta !this is the fraction of the unspotted star

          g = pow10(s% iso(i)% data(kg,j))
          R = pow10(s% iso(i)% data(kR,j))
          Teff = pow10(log_Teff)
          Tnew = pow(alfa*pow(Teff,4.0d0) + beta*pow(Teff*gamma,4.0d0),2.5d-1)
          Rnew = R*pow(Teff/Tnew,2.0d0)
          gnew = g*pow(R/Rnew,2.0d0)
          s% iso(i)% data(kR,j) = log10(Rnew)
          s% iso(i)% data(kT,j) = log10(Tnew)
          s% iso(i)% data(kg,j) = log10(gnew)
       enddo
    enddo
  end subroutine add_spots

end module spotify
