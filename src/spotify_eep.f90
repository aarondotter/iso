program eep_cmd

  use iso_eep_support
  
  implicit none

  type(track) :: t
  character(len=32) :: arg
  real(dp) :: scaling_factor
  
  call get_command_argument(1,t% filename)
  
  if(command_argument_count() > 1) then
     call get_command_argument(2,arg)
     read(arg,*) scaling_factor
  else
     scaling_factor = 0.0d0
  endif
  
  call read_eep(t,full_path=.true.,append_eep=.false.)

  call spotify(t)

  t% filename = trim(t% filename) // '_spot'

  write(*,*) ' writing ', trim(t% filename)
  call write_track(t)

contains

  subroutine spotify(t)
    type(track), intent(inout) :: t
    real(dp) :: Teff, R, logg, g, alfa, beta, Tnew, Rnew, gnew, Tspot
    integer :: i, ilogTeff, ilogR, ilogg

    ilogTeff = 0; ilogR = 0; ilogg = 0
    
    do i = 1, t% ncol
       if(trim(adjustl(t% cols(i)% name))== 'log_Teff') ilogTeff = i
       if(trim(adjustl(t% cols(i)% name))== 'log_R') ilogR = i
       if(trim(adjustl(t% cols(i)% name))== 'log_g') ilogg = i
    enddo

    do i = 1, t% ntrack
       Teff = pow10(t% tr(ilogTeff,i))
       R = pow10(t% tr(ilogR,i))
       logg = t% tr(ilogg,i)
       g = pow10(logg)
       !calculate spot parameters
       beta = fd(logg)*spot_frac(Teff)
       alfa = 1.0d0 - beta
       Tspot = scaling_factor*Teff
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

  function spot_frac(Teff) result(y)
    real(dp), intent(in) :: Teff
    real(dp) :: y, x
    real(dp), parameter :: A = -1.77514793d-7
    real(dp), parameter :: B = 1.63313609d-3
    real(dp), parameter :: C = -3.35621302d0
    real(dp), parameter :: xmin = 3.0988893333198866d3
    real(dp), parameter :: xmax = 6.1011106351334465d3
    x = Teff + 400.0_dp
    
    if(x > xmin .and. x < xmax) then
       y = A*x*x + B*x + C
    else
       y = 0.0d0
    endif
  end function spot_frac

  function fd(x) result(y)
    real(dp), intent(in) :: x
    real(dp), parameter :: x0 = 4.0d0
    real(dp), parameter :: sigma = 0.25d0 !0.2d0
    real(dp) :: y
    y = 1.0d0 - 1.0d0/(1.0d0 + exp( (x-x0)/sigma))
  end function fd
end program eep_cmd
