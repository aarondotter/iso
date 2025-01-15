module iso_eep_color

  !MESA modules
  use interp_1d_lib_sg

  !local modules
  use iso_eep_support
  use BC_tables

  implicit none

  private

  real(sp), parameter :: SolBol=4.74
  real(sp), parameter :: Z_sol = 0.0134, Xsol = 0.7381
  real(sp), parameter :: Z_div_X_sol = 0.0181
  real(sp), parameter :: log_Z_sol = log10(Z_sol)
  real(sp) :: BC_Fe_div_H
  type(BC_table), allocatable :: b(:)
  logical ::do_fixed_Z = .false.
  public color_init, write_cmds_to_file, write_track_cmd_to_file

contains

  subroutine color_init(phot_string,bc_table_list,set_fixed_Fe_div_H,Fe_div_H,ierr)
    character(len=32) :: phot_string
    character(len=file_path), intent(in) :: bc_table_list
    logical, intent(in) :: set_fixed_Fe_div_H
    real(sp), intent(in) :: Fe_div_H
    integer, intent(out) :: ierr
    call BC_table_init(phot_string,bc_table_list,b,ierr)
    if(set_fixed_Fe_div_H)then
       BC_Fe_div_H = Fe_div_H
       do_fixed_Z = .true.
    endif
    if(ierr/=0) write(0,*) 'color_init: failed to initialize BC tables'
  end subroutine color_init

  subroutine get_iso_mags(iso,log_Z_div_Zsol,ierr)
    type(isochrone), intent(inout) :: iso
    real(sp), allocatable, intent(out) :: log_Z_div_Zsol(:)
    integer, intent(out) :: ierr
    integer :: i, nb, j, n, jZ
    real(sp), allocatable :: res(:)
    real(sp) :: logT, logg, logL, X, Y, Z, FeH
    integer :: iT, ig, iL, iH, iHe3, iHe4, isurfZ, ilogR, ire, irp, iw
    !a few initializations
    ilogR=0; ire=0; irp=0
    iT=0; ig=0; iL=0; iH=0; iHe3=0; iHe4=0; jZ=0; isurfZ=0; iw=0
    
    !locate columns
    do i = 1, iso% ncol
       if(trim(adjustl(iso% cols(i)% name)) == 'log_Teff') then
          iT=i
       else if(trim(adjustl(iso% cols(i)% name)) == 'log_g') then
          ig=i
       else if(trim(adjustl(iso% cols(i)% name))== 'log_L') then
          iL=i
       else if(trim(adjustl(iso% cols(i)% name))=='surface_h1')then
          iH=i
       else if(trim(adjustl(iso% cols(i)% name))=='surface_he3')then
          iHe3=i   
       else if(trim(adjustl(iso% cols(i)% name))=='surface_he4')then
          iHe4=i   
       else if(trim(adjustl(iso% cols(i)% name))=='log_surf_z')then
          isurfZ=i
       else if(trim(adjustl(iso% cols(i)% name))=='r_equatorial_div_r')then
          ire=i
       else if(trim(adjustl(iso% cols(i)% name))=='r_polar_div_r')then
          irp=i
       else if(trim(adjustl(iso% cols(i)% name))=='log_R')then
          ilogR=i
       else if(trim(adjustl(iso% cols(i)% name))=='surf_avg_omega_div_omega_crit')then
          iw=i
       endif
    enddo
    
    iso% nfil = b(1)% num_filter
    allocate(iso% mags(iso% nfil, iso% neep),res(iso% nfil))
    allocate(iso% labels(iso% nfil))
    iso% labels = b(1)% labels
    res = 0.
    iso% mags = 0.
    nb=size(b)

    allocate(log_Z_div_Zsol(iso% neep))

    do i=1,iso% neep
       logT = real(iso% data(iT,i),kind=sp)
       logg = real(iso% data(ig,i),kind=sp)
       logL = real(iso% data(iL,i),kind=sp)

       if(do_fixed_Z)then
          FeH = BC_Fe_div_H
       else
          if(isurfZ>0)then
             FeH = real(iso% data(isurfZ,i),kind=sp) - log_Z_sol
          else
             X    = real(iso% data(iH,i),kind=sp)
             Y    = real(iso% data(iHe3,i)+iso% data(iHe4,i),kind=sp)
             Z    = 1.0 - X - Y 
             FeH = log10(Z) - log_Z_sol
          endif
       endif

       !this sets a hard upper limit to combat 3DUP
       FeH = min(FeH,real(iso% Fe_div_H,kind=sp) + 0.1)

       log_Z_div_Zsol(i) = FeH

       n=nb
       !interpolation is either linear, quadratic, or cubic
       if(FeH < b(1)% FeH .or. n==1) then
          call BC_interp_filters(b(1),logg,logT,iso% Av, iso% Rv, res,ierr)
       else if(FeH >= b(n)% FeH) then
          call BC_interp_filters(b(n),logg,logT,iso% Av, iso% Rv, res,ierr)

       else if(n==2) then
          call linear(b(1), b(2), logg, logT, iso% Av, iso% Rv, FeH, res, ierr)

       else if(n==3)then
          call quadratic(b(1), b(2), b(3), logg, logT, iso% Av, iso% Rv, FeH, res, ierr)

       else 
          do j=1,n-1
             if(FeH >= b(j)% FeH .and. FeH < b(j+1)% FeH) jZ=j
          enddo

          if(jZ==1) then
             call quadratic(b(1), b(2), b(3), logg, logT, iso% Av, iso% Rv, FeH, res, ierr)
          else if(jZ==n-1) then
             call quadratic(b(n-2), b(n-1), b(n), logg, logT, iso% Av, iso% Rv, FeH, res, ierr)
          else
             call cubic(b(jZ-1),b(jZ),b(jZ+1),b(jZ+2), logg, logT, iso% Av, iso% Rv, FeH, res, ierr)
          endif
       endif

       iso% mags(:,i) = SolBol - 2.5*logL - res
    enddo
    
  end subroutine get_iso_mags


  subroutine get_eep_mags(t,log_Z_div_Zsol,ierr)
    type(track), intent(inout) :: t
    real(sp), allocatable, intent(out) :: log_Z_div_Zsol(:)
    integer, intent(out) :: ierr
    integer :: i, nb, j, n, jZ
    real(sp), allocatable :: res(:)
    real(sp) :: logT, logg, logL, X, Y, Z, FeH
    integer :: iT, ig, iL, iH, iHe3, iHe4, isurfZ
    !a few initializations
    iT=0; ig=0; iL=0; iH=0; iHe3=0; iHe4=0; jZ=0; isurfZ=0

    !locate columns
    do i = 1, t% ncol
       if(trim(adjustl(t% cols(i)% name)) == 'log_Teff') then
          iT=i
       else if(trim(adjustl(t% cols(i)% name)) == 'log_g') then
          ig=i
       else if(trim(adjustl(t% cols(i)% name))== 'log_L') then
          iL=i
       else if(trim(adjustl(t% cols(i)% name))=='surface_h1')then
          iH=i
       else if(trim(adjustl(t% cols(i)% name))=='surface_he3')then
          iHe3=i   
       else if(trim(adjustl(t% cols(i)% name))=='surface_he4')then
          iHe4=i   
       else if(trim(adjustl(t% cols(i)% name))=='log_surf_z')then
          isurfZ=i
       endif
    enddo

    t% nfil = b(1)% num_filter
    allocate(t% mags(t% nfil, t% ntrack),res(t% nfil))
    allocate(t% labels(t% nfil))
    t% labels = b(1)% labels
    res = 0.
    t% mags = 0.
    nb=size(b)

    allocate(log_Z_div_Zsol(t% ntrack))

    do i=1,t% ntrack
       logT = real(t% tr(iT,i),kind=sp)
       logg = real(t% tr(ig,i),kind=sp)
       logL = real(t% tr(iL,i),kind=sp)
       if(do_fixed_Z)then
          FeH = BC_Fe_div_H
       else
          if(isurfZ>0)then
             FeH = real(t% tr(isurfZ,i),kind=sp) - log_Z_sol
          else
             X    = real(t% tr(iH,i),kind=sp)
             Y    = real(t% tr(iHe3,i)+t% tr(iHe4,i),kind=sp)
             Z    = 1.0 - X - Y 
             FeH = log10(Z) - log_Z_sol
          endif
       endif

       !this sets a hard upper limit to combat 3DUP
       FeH = min(FeH,real(t% Fe_div_H,kind=sp) + 0.1)

       log_Z_div_Zsol(i) = FeH
       
       n=nb
       !interpolation is either linear, quadratic, or cubic
       if(FeH < b(1)% FeH .or. n==1) then
          call BC_interp_filters(b(1),logg,logT,t% Av, t% Rv, res,ierr)
       else if(FeH >= b(n)% FeH) then
          call BC_interp_filters(b(n),logg,logT,t% Av, t% Rv, res,ierr)

       else if(n==2) then
          call linear(b(1), b(2), logg, logT, t% Av, t% Rv, FeH, res, ierr)

       else if(n==3)then
          call quadratic(b(1), b(2), b(3), logg, logT, t% Av, t% Rv, FeH, res, ierr)

       else 
          do j=1,n-1
             if(FeH >= b(j)% FeH .and. FeH < b(j+1)% FeH) jZ=j
          enddo

          if(jZ==1) then
             call quadratic(b(1), b(2), b(3), logg, logT, t% Av, t% Rv,  FeH, res, ierr)
          else if(jZ==n-1) then
             call quadratic(b(n-2), b(n-1), b(n), logg, logT, t% Av, t% Rv, FeH, res, ierr)
          else
             call cubic(b(jZ-1),b(jZ),b(jZ+1),b(jZ+2), logg, logT, t% Av, t% Rv, FeH, res, ierr)
          endif
       endif

       t% mags(:,i) = SolBol - 2.5*logL - res
    enddo

  end subroutine get_eep_mags


  subroutine linear(t1,t2,logg,logT,Av,Rv,FeH,res,ierr)
    type(bc_table), intent(inout) :: t1, t2
    real(sp), intent(in) :: logg, logT, Av, Rv, FeH
    real(sp), intent(out) :: res(:) !iso% nfil
    integer, intent(out) :: ierr
    real(sp), allocatable :: res_tmp(:,:), slope(:), intercept(:)
    integer :: n
    n=size(res)
    allocate(res_tmp(n,2),slope(n),intercept(n))
    call BC_interp_filters(t1, logg, logT, Av, Rv, res_tmp(:,1), ierr)
    call BC_interp_filters(t2, logg, logT, Av, Rv, res_tmp(:,2), ierr)
    !linear interpolation
    slope = (res_tmp(:,2)-res_tmp(:,1))/(b(2)% FeH-b(1)% FeH)
    intercept = (res_tmp(:,1)*b(2)% FeH-res_tmp(:,2)*b(1)% FeH)/(b(2)% FeH-b(1)% FeH)
    res = slope*FeH + intercept
  end subroutine linear

  subroutine quadratic(t1,t2,t3,logg,logT,Av,Rv,FeH,res,ierr)
    type(bc_table), intent(inout) :: t1,t2,t3
    real(sp), intent(in) :: logg, logT, Av, Rv, FeH
    real(sp), intent(out) :: res(:) !iso% nfil
    integer, intent(out) :: ierr
    real(sp), allocatable :: res_tmp(:,:)
    real(sp) :: x12, x23, x0
    integer :: k, n
    character(len=file_path) :: str
    n=size(res)
    allocate(res_tmp(n,3))
    x12 = t2% FeH - t1% FeH
    x23 = t3% FeH - t2% FeH
    x0  =     FeH - t1% FeH
    call BC_interp_filters(t1, logg, logT, Av, Rv, res_tmp(:,1), ierr)
    call BC_interp_filters(t2, logg, logT, Av, Rv, res_tmp(:,2), ierr)
    call BC_interp_filters(t3, logg, logT, Av, Rv, res_tmp(:,3), ierr)
    do k = 1, n
       call interp_3_to_1_sg(x12,x23,x0,res_tmp(k,1),res_tmp(k,2),res_tmp(k,3),res(k),str,ierr)
       if(ierr/=0) then
          write(*,*) trim(str)
          return
       endif
    enddo
  end subroutine quadratic

  subroutine cubic(t1,t2,t3,t4,logg,logT,Av,Rv,FeH,res,ierr)
    type(bc_table), intent(inout) :: t1,t2,t3,t4
    real(sp), intent(in) :: logg, logT, Av, Rv, FeH
    real(sp), intent(out) :: res(:) !iso% nfil
    integer, intent(out) :: ierr
    real(sp), allocatable :: res_tmp(:,:)
    real(sp) :: x12, x23, x34, x0
    integer :: k, n
    character(len=file_path) :: str
    n=size(res)
    allocate(res_tmp(n,4))
    x12 = t2% FeH - t1% FeH
    x23 = t3% FeH - t2% FeH
    x34 = t4% FeH - t3% FeH
    x0  =     FeH - t1% FeH
    call BC_interp_filters(t1, logg, logT, Av, Rv, res_tmp(:,1), ierr)
    call BC_interp_filters(t2, logg, logT, Av, Rv, res_tmp(:,2), ierr)
    call BC_interp_filters(t3, logg, logT, Av, Rv, res_tmp(:,3), ierr)
    call BC_interp_filters(t4, logg, logT, Av, Rv, res_tmp(:,4), ierr)
    do k = 1, n
       call interp_4_to_1_sg(x12,x23,x34,x0,res_tmp(k,1),res_tmp(k,2),res_tmp(k,3),res_tmp(k,4),res(k),str,ierr)
       if(ierr/=0) then
          write(*,*) trim(str)
          return
       endif
    enddo
  end subroutine cubic

  subroutine write_track_cmd_to_file(t)
    type(track), intent(inout) :: t
    character(len=file_path) :: filename
    integer :: io, ierr, iT, ig, iL, iA, i, ncol
    real(sp), allocatable :: log_Z_div_Zsol(:), Zsurf(:)
    character(len=8) :: have_phase
    ierr=0; iT=0; ig=0; iL=0; iA=0
    if(t% cmd_suffix/='') filename = trim(t% filename) // '.' // trim(t% cmd_suffix)
    write(*,*) ' writing ', trim(filename)
    open(newunit=io,file=trim(filename),action='write',status='unknown',iostat=ierr)
    if(ierr/=0) return
    if(t% has_phase)then
       have_phase = 'YES'
    else
       have_phase = 'NO'
    endif
    have_phase = adjustr(have_phase)
    
    do i=1, t% ncol
       if(trim(adjustl(t% cols(i)% name)) == 'log_Teff') then
          iT=i
       else if(trim(adjustl(t% cols(i)% name)) == 'log_g') then
          ig=i
       else if(trim(adjustl(t% cols(i)% name)) == 'log_L') then
          iL=i
       else if(trim(adjustl(t% cols(i)% name)) == 'star_age')then
          iA=i
       endif
    enddo

    call get_eep_mags(t,log_Z_div_Zsol,ierr)
    if(ierr/=0) write(0,*) ' problem in get_eep_mags '

    if(t% has_phase)then
       ncol = t% nfil + 6
    else
       ncol = t% nfil + 5
    endif

    allocate(Zsurf(size(log_Z_div_Zsol)))
    Zsurf = pow10_sg(log_Z_div_Zsol + log_Z_sol)

    write(io,'(a25,a8)') '# MIST version         : ', t% version_string
    write(io,'(a25,i8)') '# MESA revision number : ', t% MESA_revision_number
    write(io,'(a25,a)') '# photometric system    : ', b(1)% photometric_system_string
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a88)') '#  Yinit        Zinit   [Fe/H]   [a/Fe]  v/vcrit                                        '
    write(io,'(a2,f6.4,1p1e13.5,0p3f9.2)') '# ', t% initial_Y, t% initial_Z, t% Fe_div_H, t% alpha_div_Fe, t% v_div_vcrit
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a1,1x,a16,4a8,2x,a10)') '#','initial_mass', 'N_pts', 'N_EEP', 'N_col', 'phase', 'type'
    write(io,'(a1,1x,1p1e16.10,3i8,a8,2x,a10)') '#', t% initial_mass, t% ntrack, t% neep, ncol, have_phase, &
         star_label(t% star_type)
    write(io,'(a8,20i8)') '# EEPs: ', t% eep
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a25,f6.3)') '# extinction:   Av, Rv = ', t% Av, t% Rv
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'


    write(io,'(a1,i31,4i32,299(17x,i3))') '#    ', (i,i=1,ncol)

    if(t% has_phase)then
       write(io,'(a1,a31,4a32,299a20)') '#', 'star_age', 'log_Teff', &
            'log_g', 'log_L', 'Z_surf', adjustr(t% labels), 'phase'
    else
       write(io,'(a1,a31,4a32,299a20)') '#', 'star_age', 'log_Teff', &
            'log_g', 'log_L', 'Z_surf', adjustr(t% labels)
    endif
    do i=1,t% ntrack
       if(t% has_phase) then
          write(io,'(5(1pes32.16e3),299(0pf20.6))') t% tr(iA,i), t% tr(iT,i), &
               t% tr(ig,i), t% tr(iL,i), Zsurf(i), t% mags(:,i), real(t% phase(i),kind=sp)
       else
          write(io,'(5(1pes32.16e3),299(0pf20.6))') t% tr(iA,i), t% tr(iT,i), &
               t% tr(ig,i), t% tr(iL,i), Zsurf(i), t% mags(:,i)
       endif
    enddo
    close(io)
  end subroutine write_track_cmd_to_file


  subroutine write_cmds_to_file(set)
    type(isochrone_set), intent(inout) :: set
    character(len=file_path) :: output
    integer :: i, io, n, ierr
    io=98
    if(set% cmd_suffix/='') output = trim(set% filename) // '.' // trim(set% cmd_suffix)
    n=set% number_of_isochrones
    write(0,*) ' cmd output file = ', trim(output)
    open(unit=io,file=trim(output),action='write',status='unknown',iostat=ierr)
    if(ierr/=0) return

    write(io,'(a25,a8)')  '# MIST version number  = ', set% version_string
    write(io,'(a25,i8)')  '# MESA revision number = ', set% MESA_revision_number
    write(io,'(a25, a)')  '# photometric system   = ', b(1)% photometric_system_string
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a88)') '#  Yinit        Zinit   [Fe/H]   [a/Fe]  v/vcrit                                        '
    write(io,'(a2,f6.4,1p1e13.5,0p3f9.2)') '# ', set% initial_Y, set% initial_Z, set% Fe_div_H, set% alpha_div_Fe, &
         set% v_div_vcrit
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a25,i5)')   '# number of isochrones = ', n
    write(io,'(a25,f6.3)') '# extinction:   Av, Rv = ', set% Av, set% Rv
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'

    do i=1,n
       set% iso(i)% Av = set% Av
       call write_cmd_to_file(io, set% iso(i))
       if(i<n) write(io,*)
       if(i<n) write(io,*)
    enddo
    close(io)
  end subroutine write_cmds_to_file


  subroutine write_cmd_to_file(io,iso)
    integer, intent(in) :: io
    type(isochrone), intent(inout) :: iso
    character(len=32) :: age_column_header
    real(dp) :: isochrone_age
    integer :: i, iT, ig, iR, iL, ierr, iH, iMi, iMc, num_cols
    real(sp), allocatable :: log_Z_div_Zsol(:), Zsurf(:), Fe_div_H(:)
    ierr=0; iT=0; ig=0; iL=0; iMi=0; iMc=0; iH=0; iR=0

    do i=1, iso% ncol
       if(trim(adjustl(iso% cols(i)% name)) == 'log_Teff') then
          iT=i
       else if(trim(adjustl(iso% cols(i)% name)) == 'log_g') then
          ig=i
       else if(trim(adjustl(iso% cols(i)% name)) == 'log_L') then
          iL=i
       else if(trim(adjustl(iso% cols(i)% name)) == 'log_R') then
          iR=i
       else if(trim(adjustl(iso% cols(i)% name)) == 'initial_mass')then
          iMi=i
       else if(trim(adjustl(iso% cols(i)% name)) == 'star_mass')then
          iMc=i
       else if(trim(adjustl(iso% cols(i)% name)) == 'surface_h1')then
          iH=i
       endif
    enddo

    call get_iso_mags(iso,log_Z_div_Zsol,ierr)
    if(ierr/=0) write(0,*) ' problem in get_iso_mags '

    allocate(Zsurf(size(log_Z_div_Zsol)))
    Zsurf = pow10_sg(log_Z_div_Zsol + log_Z_sol)

    allocate(Fe_div_H(size(log_Z_div_Zsol)))
    if(iH==0)then
       Fe_div_H = BC_Fe_div_H
    else
       Fe_div_H = log_Z_div_Zsol - log10( real(iso% data(iH,:),kind=sp) / Xsol)
    endif

    num_cols = iso% nfil + 10
    
    if(iso% has_phase) then !add extra column for phase
       num_cols = num_cols + 1
    endif
    
    write(io,'(a25,2i5)') '# number of EEPs, cols = ', iso% neep, num_cols
    write(io,'(a1,i4,8i32,299(17x,i3))') '#    ', (i,i=1,num_cols)
 
    if(iso% age_scale==age_scale_linear)then
       age_column_header='isochrone_age_yr'
       isochrone_age = pow10(iso% age)
    else
       age_column_header='log10_isochrone_age_yr'
       isochrone_age = iso% age
    endif

    if(iso% has_phase)then
       write(io,'(a5,8a32,299a20)') '# EEP', adjustr(age_column_header), 'initial_mass', 'star_mass', 'log_Teff', &
            'log_R', 'log_g', 'log_L', '[Fe/H]_init','[Fe/H]', adjustr(iso% labels), 'phase'
    else
       write(io,'(a5,8a32,299a20)') '# EEP', adjustr(age_column_header), 'initial_mass', 'star_mass', 'log_Teff', &
            'log_R', 'log_g', 'log_L', '[Fe/H]_init', '[Fe/H]', adjustr(iso% labels)
    endif
    do i = 1,iso% neep
       if(iso% has_phase)then
          write(io,'(i5,8(1pes32.16e3),299(0pf20.6))') iso% eep(i), isochrone_age, iso% data(iMi,i), iso% data(iMc,i), &
               iso% data(iT,i), iso% data(iR,i), iso% data(ig,i), iso% data(iL,i), iso% Fe_div_H, Fe_div_H(i), iso% mags(:,i), &
               real(iso% phase(i), kind=sp)
       else
          write(io,'(i5,8(1pes32.16e3),299(0pf20.6))') iso% eep(i), isochrone_age, iso% data(iMi,i), iso% data(iMc,i), &
               iso% data(iT,i), iso% data(ig,i), iso% data(iL,i), iso% Fe_div_H, Fe_div_H(i), iso% mags(:,i)
       endif
    enddo
  end subroutine write_cmd_to_file

 
end module iso_eep_color
