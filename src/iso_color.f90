module iso_color

  !MESA modules
  use const_def, only: sp
  use utils_lib, only: alloc_iounit, free_iounit

  !other modules
  use color_def
  use color_lib

  !local modules
  use iso_eep_support

  implicit none

  type(bc_table) :: bc
  character(len=file_path) :: color_input = 'input.cmd', bc_table_list, color_suffix

contains

  subroutine read_color_input(ierr)
    integer, intent(out) :: ierr
    integer :: io
    type(bc_table), allocatable :: bc0(:)
    real(sp) :: Av, Rv
    io=alloc_iounit(ierr)
    if(ierr/=0) return
    open(io,file=trim(color_input),action='read',status='old',iostat=ierr)
    if(ierr/=0) return
    read(io,'(a)') bc_table_list
    read(io,'(a)') color_suffix
    read(io,'(3x,f6.3)') Av
    read(io,'(3x,f6.3)') Rv
    close(io)
    call free_iounit(io)
    call color_init(bc_table_list,bc0,ierr)
    call color_create_fixed_Av_Rv(bc0(1),bc,Av,Rv)
    deallocate(bc0)
  end subroutine read_color_input

  subroutine get_mags(iso,iT,ig,iL)
    type(isochrone), intent(inout) :: iso
    integer :: i, ierr
    real(sp), allocatable :: res(:)
    integer :: iT, ig, iL
    real(sp) :: logT, logg, logL
    iso% nfil = bc% num_filter
    allocate(iso% mags(iso% nfil, iso% neep),res(iso% nfil))
    allocate(iso% labels(iso% nfil))
    iso% labels = bc% labels
    iso% Av = bc% Av(1)
    iso% Rv = bc% Rv(1)
    res = 0.
    iso% mags = 0.0         
    do i=1,iso% neep
       logT = real(iso% data(iT,i),kind=sp)
       logg = real(iso% data(ig ,i),kind=sp)
       logL = real(iso% data(iL ,i),kind=sp)
       call color_get(bc, logT, logg, 1, 1, res, ierr)
       iso% mags(:,i) = SolBol - 2.5*logL - res
    enddo
    deallocate(res)
  end subroutine get_mags

end module iso_color
