module iso_color

  !MESA modules
  use const_def, only: sp
  use utils_lib, only: alloc_iounit, free_iounit

  !local modules
  use iso_eep_support
  use BC_tables
  
  implicit none

  private

  real(sp), parameter :: SolBol=4.74
  type(BC_table), allocatable :: b(:)

  public read_color_input, get_mags

contains

  subroutine read_color_input(s,ierr)
    type(isochrone_set), intent(inout) :: s
    integer, intent(out) :: ierr
    integer :: io
    character(len=file_path) :: bc_table_list
    io=alloc_iounit(ierr)
    if(ierr/=0) return
    open(io,file=trim('input.cmd'),action='read',status='old',iostat=ierr)
    if(ierr/=0) return
    read(io,'(a)') bc_table_list
    read(io,'(a)') s% cmd_suffix
    read(io,'(3x,f6.3)') s% Av
    close(io)
    call free_iounit(io)
    call BC_table_init(bc_table_list,b,ierr)
  end subroutine read_color_input

  subroutine get_mags(iso,iT,ig,iL)
    type(isochrone), intent(inout) :: iso
    integer :: i, ierr
    real(sp), allocatable :: res(:)
    integer :: iT, ig, iL
    real(sp) :: logT, logg, logL
    iso% nfil = b(1)% num_filter
    allocate(iso% mags(iso% nfil, iso% neep),res(iso% nfil))
    allocate(iso% labels(iso% nfil))
    iso% labels = b(1)% labels
    res = 0.
    iso% mags = 0.0         
    do i=1,iso% neep
       logT = real(iso% data(iT,i),kind=sp)
       logg = real(iso% data(ig,i),kind=sp)
       logL = real(iso% data(iL,i),kind=sp)
       call BC_interp_filters(b(1),logg,logT,iso% Av, res,ierr)
       iso% mags(:,i) = SolBol - 2.5*logL - res
    enddo
    deallocate(res)
  end subroutine get_mags

end module iso_color
