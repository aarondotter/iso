program blend_eeps

  !local modules
  use iso_eep_support

  implicit none

  integer :: i, ierr, io, num
  character(len=file_path) :: input_file, output
  character(len=file_path), allocatable :: eep_files(:)
  real(dp), allocatable :: weights(:)
  type(track), allocatable :: s_old(:)
  type(track) :: s_new

  ierr=0
  if(command_argument_count()<1) then
     write(*,*) '   blend_eeps                 '
     write(*,*) '   usage: ./blend_eeps <input>'
     write(*,*) '   ERROR:                     '
     write(*,*) '   no input argument supplied '
     stop
  endif

  call read_input

  allocate(s_old(num))

  !process existing files
  do i=1,num
     s_old(i)% filename = trim(eep_files(i))
     call read_eep(s_old(i))
  enddo

  s_new = eep_blender(s_old,weights)
  s_new% filename = trim(eep_dir) // '/' // trim(output)

  call write_track(s_new)

  deallocate(s_old)

contains

  subroutine read_input
    call get_command_argument(1,input_file)
    open(newunit=io,file=trim(input_file),status='old')
    read(io,*) !skip first line
    read(io,'(a)') eep_dir
    read(io,*) !skip comment line
    read(io,*) num
    read(io,*) !skip comment line
    allocate(eep_files(num))
    do i=1,num
       read(io,'(a)',iostat=ierr) eep_files(i)
       if(ierr/=0) exit
    enddo
    read(io,*) !skip comment line
    allocate(weights(num))
    do i=1,num
       read(io,*) weights(i)
    enddo

    !make sure weights sum to 1.
    weights = weights/sum(weights)

    read(io,*) !skip comment line
    read(io,'(a)') output
    close(io)
    write(*,*) ' blending ', num, ' tracks:'
    write(*,'(99a32)') eep_files
    write(*,*) weights, ' => ', sum(weights)
  end subroutine read_input

  function eep_blender(s_old,weights) result(s)
    type(track), intent(in) :: s_old(:)
    real(dp), intent(in) :: weights(:)
    integer :: j,k,n
    type(track) :: s
    !make sure everything is OK
    if(size(weights) /= size(s_old)) stop 'size(weights) != size(s_old)'
    n = size(weights)
    do k = 2,n
       if(abs(s_old(1)% initial_mass - s_old(k)% initial_mass) > mass_eps)then
          write(*,*) ' initial masses not equal: '
          write(*,*) s_old(:)% initial_mass
       endif
       !other checks go here:
    enddo

    s% initial_mass = s_old(1)% initial_mass
    s% Fe_div_H = s_old(1)% Fe_div_H
    s% alpha_div_Fe = s_old(1)% alpha_div_Fe
    s% version_string = s_old(1)% version_string
    s% initial_Z = s_old(1)% initial_Z
    s% initial_Y = s_old(1)% initial_Y
    s% v_div_vcrit = s_old(1)% v_div_vcrit
    s% has_phase = s_old(1)% has_phase
    s% ncol = s_old(1)% ncol
    allocate(s% cols(s% ncol))
    s% cols(:)% type = s_old(1)% cols(:)% type
    s% cols(:)% loc  = s_old(1)% cols(:)% loc
    s% cols(:)% name = s_old(1)% cols(:)% name
    s% neep = s_old(1)% neep
    s% ntrack = s_old(1)% ntrack
    s% MESA_revision_number = s_old(1)% MESA_revision_number
    s% star_type = s_old(1)% star_type
    allocate(s% eep(s% neep))
    s% eep = s_old(1)% eep
    allocate(s% tr(s% ncol, s% ntrack), s% phase(s% ntrack))
    if(s_old(1)% has_phase)then
       s% phase = s_old(1)% phase
    endif

    s% tr = 0d0

    do j = 1, s% ntrack
       do k=1,n
          s% tr(:,j) = s% tr(:,j) + weights(k)*s_old(k)% tr(:,j)
       enddo
    enddo
  end function eep_blender

end program blend_eeps
