program eep_cmd

  use iso_eep_support
  
  implicit none

  integer :: ierr = 0
  type(track) :: t
  real(dp) :: alpha, beta
  
  call get_command_argument(1,t% filename)
  
  call read_eep(t,full_path=.true.,append_eep=.false.)

  call spotify(t)

  t% filename = trim(t% filename) // '_spot'
 
  call write_track_cmd_to_file(t)

contains

  subroutine spotify(t)
    type(track), intent(inout) :: t
    real(dp) :: Teff, R, g, spot, alfa, beta, Tnew, Rnew, gnew
    integer :: i, ilogTeff, ilogR, ilogg
    
    allocate(Teff(t% neep), R(t% neep))

    do i = 1, t% ncol
       if(trim(t% label(i)) == 'log_Teff') ilogTeff = i
       if(trim(t% label(i)) == 'log_R') ilogR = i
       if(trim(t% label(i)) == 'log_g') ilogg = i
    enddo

    do i = 1, t% neep
       Teff = pow10(t% tr(iTeff,i))
       R = pow10(t% tr(iR,i))
       g = pow10(t% tr(ilogg,i)
       !calculate spot parameters
       beta = spot_frac(Teff)
       alfa = 1.0d0 - beta
       Tspot = 0.0d0
       !caculate adjusted Teff, R, log(g)
       Tnew = pow(alfa*powi(Teff,4) + beta*powi(Tspot,4), 0.25d0)
       Rnew = R * powi(Teff/Tnew,2)
       gnew = g * powi(Tnew/Teff,4)
       !put the results back into the track
       t% tr(ilogR,i) = log10(Rnew)
       t% tr(ilogTeff,i) = log10(Tnew)
       t% tr(ilogg,i) = log10(g)
    enddo
    
  end subroutine spotify

  real(dp) function spot_frac(Teff) result(y)
    real(dp), intent(in) :: Teff
    real(dp) :: y, frac
    real(dp), parameter :: A = -1.77514793d-7
    real(dp), parameter :: B = 1.63313609d-3
    real(dp), parameter :: C = -3.35621302d0
    real(dp), parameter :: Tmin = 3.0988893333198866d3
    real(dp), parameter :: Tmax = 6.1011106351334465d3
    if(Teff > Tmin .and. Teff < Tmax) then
       frac = A*Teff*Teff + B*Teff + C
    else
       frac = 0.0d0
    endif
    y = 1.0d0 - frac
  end function spot_frac
  
end program eep_cmd
