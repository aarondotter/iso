program iso_phase_fix_v12

  use iso_eep_support

  implicit none
  
  type(isochrone_set) :: s
  integer :: i, j, ierr, iH, iMinit, ilogT
  character(len=file_path) :: output

  !if(command_argument_count()<3)then
  !   write(0,*) "*** Charlie's Angel ***"
  !   write(0,*) '   ./charlies_angel [eep file] [iso file] [output]'
  !   write(0,*) '   appends points to [iso file] from [eep file] at each age'
  !   write(0,*) '   [output] is the new isochrone file name'
  !   stop '   badness: not enough command arguments provided'
  !endif
  
  !call get_command_argument(1,m% filename)
  !call read_eep(m,.true.)

  call get_command_argument(1,s% filename)
  call read_isochrone_file(s,ierr)
  if(ierr/=0) stop 'failed reading isochrone file'

  iMinit = find_col('initial_mass')
  iH     = find_col('surface_h1')
  ilogT  = find_col('log_Teff')
  
  do i=1, s% number_of_isochrones
     do j=1,s% iso(i)% neep
        s% iso(i)% phase(j) = eep_to_phase(s% iso(i)% eep(j))
        !except override for W-R stars
        if(s% iso(i)% data(ilogT,j) > 4d0 .and. s% iso(i)% data(iH,j) <= 3d-1 .and. &
             s% iso(i) % data(iMinit,j) >= 1d1) s% iso(i)% phase(j) = 9d0
     enddo
  enddo

  call write_isochrones_to_file(s)


  !cols hold columns names
  !want initial mass, logTeff and surface H

contains

  function find_col(col_name) result(icol)
    character(len=*) :: col_name
    character(len=20) :: iso_col_name
    integer :: icol, k
    icol = 0
    do k=1,s% iso(1)% ncol
       if( adjustl(adjustr( s% iso(1)% cols(k)% name )) == trim(col_name)) icol=k
    enddo 
  end function find_col
  
  function eep_to_phase(eep) result(ph)
    integer, intent(in) :: eep
    real(dp) :: ph

    if (eep < 202) then
       ph = -1d0
    elseif (eep >= 202 .and. eep < 454) then
       ph = 0d0
    elseif (eep >= 454 .and. eep < 605) then
       ph = 2d0
    elseif (eep >= 605 .and. eep < 707) then
       ph = 3d0
    elseif (eep >= 707 .and. eep < 808) then
       ph = 4d0
    elseif (eep >= 808 .and. eep < 1409) then
       ph = 5d0
    elseif (eep >= 1409) then
       ph = 6d0
    endif
    
  end function eep_to_phase
    
  
end program iso_phase_fix_v12
