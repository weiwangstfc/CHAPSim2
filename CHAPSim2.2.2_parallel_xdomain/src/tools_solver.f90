module solver_tools_mod

  ! procedure
  private

  public  :: Calculate_massflux_from_velocity
  public  :: Check_cfl_convection
  public  :: Check_cfl_diffusion

  public  :: Update_RePrGr
  public  :: Calculate_xz_mean
  public  :: Calculate_xzmean_perturbation
  public  :: Check_maximum_velocity
  public  :: Get_volumetric_average_3d

contains
!===============================================================================
!> \brief The main code for initializing flow variables
!>
!> This subroutine is called once in \ref Initialize_chapsim.
!>
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[inout]  
!===============================================================================
  subroutine Update_RePrGr(dm, iter, fl, tm)
    use parameters_constant_mod
    use input_thermo_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in   ) :: dm
    integer(4),     intent(in   ) :: iter  
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm
    

    real(WP) :: u0, rtmp
  !-------------------------------------------------------------------------------
  !  1/Re                                   
  !-------------------------------------------------------------------------------
    if(iter < fl%nIterIniRen) then
      fl%rre = ONE / fl%renIni
    else
      fl%rre = ONE / fl%ren
    end if

    if(dm%ithermo == 1) then
  !-------------------------------------------------------------------------------
  !  1/(Re*Pr)                                   
  !-------------------------------------------------------------------------------
      tm%rPrRen = fl%rre * tpRef0%k / tpRef0%m / tpRef0%cp
  !-------------------------------------------------------------------------------
  !  gravity force                          
  !-------------------------------------------------------------------------------  
      u0 = ONE / fl%rre * tpRef0%m / tpRef0%d / tm%lenRef
      rtmp = tm%lenRef / u0 / u0 * GRAVITY
      fl%fgravity(:) = ZERO
      if (tm%igravity == 1 ) then ! flow/gravity same dirction - x
        fl%fgravity(1) =  rtmp
      else if (tm%igravity == 2 ) then ! flow/gravity same dirction - y
        fl%fgravity(2) =  rtmp
      else if (tm%igravity == 3 ) then ! flow/gravity same dirction - z
        fl%fgravity(3) =  rtmp
      else if (tm%igravity == -1 ) then ! flow/gravity opposite dirction - x
        fl%fgravity(1) =  - rtmp
      else if (tm%igravity == -2 ) then ! flow/gravity opposite dirction - y
        fl%fgravity(2) =  - rtmp
      else if (tm%igravity == -3 ) then ! flow/gravity opposite dirction - z
        fl%fgravity(3) =  - rtmp
      else ! no gravity
        fl%fgravity(:) = ZERO
      end if
      
    end if

    return
  end subroutine Update_RePrGr
!===============================================================================
!> \brief The main code for initializing flow variables
!>
!> not changing storage position, exclude b.c. values, for example, developing
!> flow.
!> MPI : x-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[inout]  none          NA
!===============================================================================
  subroutine Calculate_xz_mean(var, dtmp, varxz_work)
    use mpi_mod
    use udf_type_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP),          intent(in) :: var(:, :, :)
    real(WP),       intent(inout) :: varxz_work(:)

    real(wp) :: varxz( size(varxz_work) )
    integer(4) :: jj, kk, ii, ny, i, j, k
    integer(4) :: ist, ien, jst, jen, kst, ken
    integer(4) :: xst(3), xen(3), xsz(3)
!-------------------------------------------------------------------------------
!   Default X-pencil
!-------------------------------------------------------------------------------
    xst(:) = dtmp%xst(:)
    xen(:) = dtmp%xen(:)
    xsz(:) = dtmp%xsz(:)
    ist = 1
    ien = xsz(1)
    jst = 1
    jen = xsz(2)
    kst = 1
    ken = xsz(3)

    varxz(:) = ZERO
    varxz_work(:) = ZERO
    nk = 0
    ni = 0
    do j = jst, jen
      jj = j - 1 + xst(2)
      do k = kst, ken
        nk = nk + 1
        do i = ist, ien
          ni = ni + 1
          varxz(jj) = varxz(jj) + var(i, j, k) !
        end do
      end do
    end do
    varxz(:) = varxz(:) / real(nk * ni, wp)

    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(varxz, varxz_work, size(varxz_work), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    varxz_work(:) = varxz_work(:) / real(ncol, wp)

    return
  end subroutine
!===============================================================================
!> \brief : 
!> MPI : x-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[inout]            
!===============================================================================
  subroutine Calculate_xzmean_perturbation(var, dtmp, varxz, varxz_shift)
    use mpi_mod
    use udf_type_mod
    implicit none
    type(DECOMP_INFO),  intent(in) :: dtmp
    real(WP),           intent(inout) :: var(:, :, :)
    real(WP),           intent(in)    :: varxz(:)
    real(WP), optional, intent(in)    :: varxz_shift(:)

    integer(4) :: jj, i, j, k
    integer(4) :: ist, ien, jst, jen, kst, ken
    integer(4) :: xst(3), xen(3), xsz(3)

!-------------------------------------------------------------------------------
!   Default X-pencil
!-------------------------------------------------------------------------------
    xst(:) = dtmp%xst(:)
    xen(:) = dtmp%xen(:)
    xsz(:) = dtmp%xsz(:)
    ist = 1
    ien = xsz(1)
    jst = 1
    jen = xsz(2)
    kst = 1
    ken = xsz(3)

    do j = jst, jen
      jj = j - 1 + xst(2)
      do k = kst, ken
        do i = ist, ien
          if( present(varxz_shift) ) then
            var(:, j, :) = var(:, j, :) - varxz(jj) + varxz_shift(jj)
          else
            var(:, j, :) = var(:, j, :) - varxz(jj)
          end if
        end do
      end do
    end do

    return
  end subroutine
!===============================================================================
!===============================================================================
!> \brief Calculate the conservative variables from primary variable.     
!>
!> This subroutine is called to update $\rho u_i$ from $u_i$.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!> \param[in]     f             flow
!_______________________________________________________________________________
  subroutine Calculate_massflux_from_velocity(f, d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod
    use operations
    implicit none
    type(t_domain), intent(in )   :: d
    type(t_flow  ), intent(inout) :: f
    real(WP), dimension( dm%nc(1) ) :: fix
    real(WP), dimension( dm%np(1) ) :: fox
    real(WP), dimension( dm%nc(2) ) :: fiy
    real(WP), dimension( dm%np(2) ) :: foy
    real(WP), dimension( dm%nc(3) ) :: fiz
    real(WP), dimension( dm%np(3) ) :: foz
    real(WP), dimension( dm%uy_ysz(1), dm%uy_ysz(2), dm%uy_ysz(3)) ::  uy_ypencil
    real(WP), dimension( dm%uy_ysz(1), dm%uy_ysz(2), dm%uy_ysz(3)) :: duy_ypencil
    real(WP), dimension( dm%uz_ysz(1), dm%uz_ysz(2), dm%uz_ysz(3)) ::  uz_ypencil
    real(WP), dimension( dm%uz_ysz(1), dm%uz_ysz(2), dm%uz_ysz(3)) :: duz_ypencil
    real(WP), dimension( dm%ps_ysz(1), dm%ps_ysz(2), dm%ps_ysz(3)) ::   d_ypencil
    
    real(WP), dimension( dm%uz_zsz(1), dm%uz_zsz(2), dm%uz_zsz(3)) ::  uz_zpencil
    real(WP), dimension( dm%uz_zsz(1), dm%uz_zsz(2), dm%uz_zsz(3)) :: duz_zpencil
    real(WP), dimension( dm%ps_zsz(1), dm%ps_zsz(2), dm%ps_zsz(3)) ::   d_zpencil
    
    integer(4) :: i, j, k
!-------------------------------------------------------------------------------
! Default x-pencil
!-------------------------------------------------------------------------------
    if(dm%ux_xsz(1) /= dm%np(1)) call Print_error_msg("Error, not X-pencil")
!-------------------------------------------------------------------------------
! x-pencil : u1 -> g1
!-------------------------------------------------------------------------------
    do k = 1, dm%ux_xsz(3)
      do j = 1, dm%ux_xsz(2)
        fix(:) = f%dDens(:, j, k)
        call Get_midp_interpolation_1D( 'x', 'C2P', d, fix(:), fox(:) )
        f%gx(:, j, k) = fox(:) * f%qx(:, j, k)
      end do
    end do
!-------------------------------------------------------------------------------
! x-pencil --> y-pencil
!-------------------------------------------------------------------------------
    call transpose_x_to_y(f%qy,    uy_ypencil, dm%dcpc)
    call transpose_x_to_y(f%dDens,  d_ypencil, dm%dccc)
    call transpose_x_to_y(f%qz,    uz_ypencil, dm%dccp)
!-------------------------------------------------------------------------------
! y-pencil : u2 -> g2
!-------------------------------------------------------------------------------
    do k = 1, dm%uy_ysz(3)
      do i = 1, dm%uy_ysz(1)
        fiy(:) = d_ypencil(i, :, k)
        call Get_midp_interpolation_1D( 'y', 'C2P', d, fiy(:), foy(:) )
        duy_ypencil(i, :, k) = foy(:) * uy_ypencil
      end do
    end do
!-------------------------------------------------------------------------------
! y-pencil --> z-pencil
!-------------------------------------------------------------------------------
    call transpose_y_to_z( d_ypencil,  d_zpencil, dm%dccc)
    call transpose_y_to_z(uz_ypencil, uz_zpencil, dm%dccp)
!-------------------------------------------------------------------------------
! Z-pencil : u3 -> g3
!-------------------------------------------------------------------------------
    do j = 1, dm%uz_zsz(2)
      do i = 1, dm%uz_zsz(1)
        fiz(:) = d_zpencil(i, j, :)
        call Get_midp_interpolation_1D( 'z', 'C2P', d, fiz(:), foz(:) )
        duz_zpencil(i, j, :) = foz(:) * uz_zpencil(i, j, :)
      end do
    end do
!-------------------------------------------------------------------------------
! z-pencil --> y-pencil
!-------------------------------------------------------------------------------
    call transpose_z_to_y(duz_zpencil, duz_ypencil, dm%dccp)
!-------------------------------------------------------------------------------
! y-pencil --> x-pencil
!-------------------------------------------------------------------------------
    call transpose_y_to_x(duz_ypencil, f%gz, dm%dccp)
    call transpose_y_to_x(duy_ypencil, f%gy, dm%dcpc)

    return
  end subroutine Calculate_massflux_from_velocity
!===============================================================================
!> \brief : 
!> MPI : x-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[inout]         
!===============================================================================
  subroutine Check_cfl_diffusion(x2r, rre)
    use input_general_mod, only : dt
    use parameters_constant_mod, only : TWO, ONE
    use precision_mod
    implicit none
    real(WP), intent(in) :: x2r(3)
    real(WP), intent(in) :: rre
    real(WP) :: cfl_diff

    ! check, ours is two times of the one in xcompact3d.
    cfl_diff = sum(x2r) * TWO * dt * rre

    if(cfl_diff > ONE) call Print_warning_msg("Warning: Diffusion number is larger than 1.")
    write (OUTPUT_UNIT,*) "  Diffusion number :"
    write (OUTPUT_UNIT,"(12X, F13.8)") cfl_diff
    
    return
  end subroutine
!===============================================================================
!> \brief : 
!> MPI : x-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[inout]         
!===============================================================================
  subroutine Check_cfl_convection(u, v, w, d)
    use parameters_constant_mod, only : ZERO, ONE
    use precision_mod
    use input_general_mod, only : dt
    use udf_type_mod, only : t_domain
    use operations, only : Get_midp_interpolation_1D
    use domain_decomposition_mod
    implicit none

    type(t_domain),               intent(in) :: d
    real(WP), dimension(:, :, :), intent(in) :: u, v, w

    real(WP), allocatable :: fi(:), fo(:)
    real(WP) :: udx_xpencil (dm%ps_xsz(1), dm%ps_xsz(2), dm%ps_xsz(3))
    real(WP) :: udx_ypencil (dm%ps_ysz(1), dm%ps_ysz(2), dm%ps_ysz(3))
    real(WP) :: udx_ypencil (dm%ps_zsz(1), dm%ps_zsz(2), dm%ps_zsz(3))
    real(WP) ::   v_ypencil (dm%uy_ysz(1), dm%uy_ysz(2), dm%uy_ysz(3))
    real(WP) ::   w_ypencil (dm%uz_ysz(1), dm%uz_ysz(2), dm%uz_ysz(3))
    real(WP) ::   w_zpencil (dm%uz_zsz(1), dm%uz_zsz(2), dm%uz_zsz(3))
    real(WP)   :: cfl_convection, cfl_convection_work
    integer(4) :: i, j, k


!-------------------------------------------------------------------------------
!   Ensure it is in x-pencil
!-------------------------------------------------------------------------------
    if(dm%ux_xsz(1) /= dm%np(1)) call Print_error_msg("Error, not X-pencil")
!-------------------------------------------------------------------------------
! X-pencil
!-------------------------------------------------------------------------------
    allocate ( udx_xpencil( dm%ps_xsz(1), dm%ps_xsz(2), dm%ps_xsz(3) ) )
    udx_xpencil = ZERO
!-------------------------------------------------------------------------------
! X-pencil : \overline{u}^x/dx at cell centre
!-------------------------------------------------------------------------------
    allocate ( fi( dm%np(1) ) ); fi = ZERO
    allocate ( fo( dm%nc(1) ) ); fo = ZERO
    udx_pencil(:, :, :) = ZERO
    do k = 1, dm%ux_xsz(3)
      do j = 1, dm%ux_xsz(2)
        fi(:) = u(:, j, k)
        call Get_midp_interpolation_1D('x', 'P2C', d, fi(:), fo(:))
        udx_xpencil(:, j, k) = fo(:) * dm%h1r(3) * dt
      end do
    end do
    deallocate (fi)
    deallocate (fo)
!-------------------------------------------------------------------------------
! Convert X-pencil to Y-Pencil
!-------------------------------------------------------------------------------
    call transpose_x_to_y(udx_xpencil, udx_ypencil, dm%dccc)
    call transpose_x_to_y(v,             v_ypencil, dm%dcpc)
    call transpose_x_to_y(w,             w_ypencil, dm%dccp)
!-------------------------------------------------------------------------------
! Y-pencil : \overline{v}^y/dy at cell centre
!-------------------------------------------------------------------------------
    allocate ( fi( dm%np(2) ) ); fi = ZERO
    allocate ( fo( dm%nc(2) ) ); fo = ZERO
    do k = 1, dm%uy_ysz(3)
      do i = 1, dm%uy_ysz(1)
        fi(:) = v(i, :, k)
        call Get_midp_interpolation_1D('y', 'P2C', d, fi(:), fo(:))
        udx_ypencil(i, :, k) = udx_ypencil(i, :, k) + fo(:) * dm%h1r(2) * dt
      end do
    end do
    deallocate (fi)
    deallocate (fo)
!-------------------------------------------------------------------------------
! Convert Y-pencil to Z-Pencil
!-------------------------------------------------------------------------------
    call transpose_y_to_z(udx_ypencil, udx_zpencil, dm%dccc)
    call transpose_y_to_z(  w_ypencil,   w_zpencil, dm%dccp)
!-------------------------------------------------------------------------------
! Z-pencil : \overline{w}^z/dz at cell centre
!-------------------------------------------------------------------------------
    allocate ( fi( dm%np(3) ) ); fi = ZERO
    allocate ( fo( dm%nc(3) ) ); fo = ZERO
    do j = 1, dm%uz_zsz(2)
      do i = 1, dm%uz_zsz(1)
        fi(:) = w_zpencil(i, j, :)
        call Get_midp_interpolation_1D('z', 'P2C', d, fi(:), fo(:))
        udx_zpencil(i, j, :) = udx_zpencil(i, j, :) + fo(:) * dm%h1r(3) * dt
      end do
    end do
    deallocate (fi)
    deallocate (fo)
!-------------------------------------------------------------------------------
! Z-pencil : Find the maximum 
!-------------------------------------------------------------------------------
    cfl_convection = MAXVAL(udx_zpencil(:, :, :))
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(cfl_convection, cfl_convection_work, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)

    if(nrank == 0) then
      if(cfl_convection_work > ONE) call Print_warning_msg("Warning: CFL is larger than 1.")
      write (OUTPUT_UNIT,*) "  CFL (convection) :"
      write (OUTPUT_UNIT,"(12X, F13.8)") cfl_convection_work
    end if
    
    return
  end subroutine

!===============================================================================
!>\brief : to calculate:
!>         fo = \int_1^nx \int_
!> This is based only y-direction stretching.
!> Here is 2nd order Trapezoid Method. Need to improve! Check!
!> MPI : 
!>     default x-pencil
!>     working in : y-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[inout]         
!===============================================================================
  subroutine Get_volumetric_average_3d(fi3d, str, d, fo_work)
    ! how to get a high order bulk value?
    use parameters_constant_mod, only : ZERO, HALF
    use udf_type_mod,            only : t_domain
    implicit none
  
    type(t_domain), intent(in) :: d
    real(WP),       intent(in) :: fi3d(:, :, :)
    real(WP),      intent(out) :: fo_work
    character(2),   intent(in) :: str
 
    type(DECOMP_INFO) :: decomp
    real(WP), allocatable   :: fo3dy_ypencil(:, :, :)
    real(WP), allocatable   :: fi3d_ypencil(:, :, :)
    real(WP)   :: vol, fo
    integer(4) :: i, j, k
    integer(4) :: nix, niy, niz
    integer(4) :: ncy
!-------------------------------------------------------------------------------
!   transpose to y pencil. Default is x-pencil.
!-------------------------------------------------------------------------------
    if(str=='ux') then
      if(dm%ux_xsz /= dm%np(1)) call Print_error_msg("Error, not X-pencil")
      ysz(1:3) = dm%ux_ysz(1:3)
      decomp = dm%dpcc
    else if(str=='uy') then
      if(dm%uy_xsz /= dm%nc(2)) call Print_error_msg("Error, not X-pencil")
      ysz(1:3) = dm%uy_ysz(1:3)
      decomp = dm%dcpc
    else if(str=='uz') then
      if(dm%uz_xsz /= dm%nc(3)) call Print_error_msg("Error, not X-pencil")
      ysz(1:3) = dm%uz_ysz(1:3)
      decomp = dm%dccp
    else if(str=='ps') then
      if(dm%ps_xsz /= dm%nc(3)) call Print_error_msg("Error, not X-pencil")
      ysz(1:3) = dm%ps_ysz(1:3)
      decomp = dm%dccc
    else
      call Print_error_msg("No such variables defined.")
    end if
    allocate ( fi3d_ypencil(ysz(1), ysz(2), ysz(3)) )
    fi3d_ypencil = ZERO

    call transpose_x_to_y(fi3d, fi3d_ypencil, decomp)
!-------------------------------------------------------------------------------
!   In Y-pencil now
!-------------------------------------------------------------------------------
    if(str=='uy')then
!-------------------------------------------------------------------------------
!   if variable is stored in y-nodes, extend them to y-cell centres
!   for example, uy
!-------------------------------------------------------------------------------
      if( dm%is_periodic(2) ) then
        noy = ysz(2)
      else
        noy = ysz(2) - 1
      end if
      allocate( fo3dy_ypencil(ysz(1), noy, ysz(3)) )
      fo3dy = ZERO
      call Get_y_midp_P2C_3dArray ( fi3d_ypencil, d, fo3dy_ypencil)
      fo = ZERO
      vol = ZERO
      do k = 1, ysz(3)
        do i = 1, ysz(1)
          do j = 1, noy
            ! fo = fo + &
            !     ( dm%yp(j + 1) - dm%yp(j) ) / SIX * &
            !     ( fi3d_ypencil(i, j, k) + &
            !       FOUR * fo3dy(i, j, k) + &
            !       fi3d_ypencil(i, dm%jNeighb(3, j), k)) ! Simpson 2nd order 
            fo = fo + &      
                ( fi3d_ypencil(i, dm%jNeighb(3, j), k) + fo3dy_ypencil(i, j, k) ) * &
                ( dm%yp(j + 1) - dm%yc(j) ) * HALF + &
                ( fi3d_ypencil(i, j,               k) + fo3dy_ypencil(i, j, k) ) * &
                ( dm%yc(j    ) - dm%yp(j) ) * HALF
            vol = vol + ( dm%yp(j + 1) - dm%yp(j) )
          end do
        end do
      end do
      deallocate(fo3dy)
    else
!-------------------------------------------------------------------------------
!   if variable is not stored in y-nodes, extends them to y-nodes.
!   for example, ux, density, etc.
!-------------------------------------------------------------------------------
      if( dm%is_periodic(2) ) then
        noy = ysz(2)
      else
        noy = ysz(2) + 1
      end if
      allocate( fo3dy_ypencil(ysz(1), noy, ysz(3)) )
      fo3dy = ZERO
      call Get_y_midp_C2P_3dArray ( fi3d_ypencil, d, fo3dy_ypencil)
      fo = ZERO
      vol = ZERO
      do k = 1, ysz(3)
        do i = 1, ysz(1)
          do j = 1, ysz(2)
            fo = fo + &
                ( fo3dy_ypencil(i, dm%jNeighb(3, j), k) + fi3d_ypencil(i, j, k) ) * &
                ( dm%yp(j + 1) - dm%yc(j) ) * HALF + &
                ( fo3dy_ypencil(i, j,               k) + fi3d_ypencil(i, j, k) ) * &
                ( dm%yc(j    ) - dm%yp(j) ) * HALF
            vol = vol + ( dm%yp(j + 1) - dm%yp(j) )
          end do
        end do
      end do
      deallocate(fo3dy_ypencil)
    end if
    deallocate(fi3d_ypencil)
    
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce( fo,  fo_work, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    call mpi_allreduce(vol, vol_work, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    fo_work = fo_work / vol_work

    if(nrank == 0) then
      Call Print_debug_mid_msg("  The bulk value is:")
      write (OUTPUT_UNIT, '(5X, A, 1ES13.5)') 'Variable bulk : ', fo_work
    end if

    return 
  end subroutine Get_volumetric_average_3d

end module
