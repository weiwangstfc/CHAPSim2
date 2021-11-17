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
    type(DECOMP_INFO),  intent(in)    :: dtmp
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
!> \brief Calculate the conservative variables from primary variable.     
!> This subroutine is called to update $\rho u_i$ from $u_i$.
!------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    needed         specified   pubic
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     dm             domain
!> \param[in]     fm             flow
!===============================================================================
  subroutine Calculate_massflux_from_velocity(dm, fl)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod
    use operations
    implicit none
    type(t_domain), intent(in )   :: dm
    type(t_flow  ), intent(inout) :: fl
    real(WP), dimension( dm%nc(1) ) :: fix
    real(WP), dimension( dm%np(1) ) :: fox
    real(WP), dimension( dm%nc(2) ) :: fiy
    real(WP), dimension( dm%np(2) ) :: foy
    real(WP), dimension( dm%nc(3) ) :: fiz
    real(WP), dimension( dm%np(3) ) :: foz
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) ::  uy_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: duy_ypencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) ::  uz_ypencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) :: duz_ypencil
    real(WP), dimension( dm%dccc%ysz(1), &
                         dm%dccc%ysz(2), &
                         dm%dccc%ysz(3)) ::   d_ypencil
    
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) ::  uz_zpencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: duz_zpencil
    real(WP), dimension( dm%dccc%zsz(1), &
                         dm%dccc%zsz(2), &
                         dm%dccc%zsz(3)) ::   d_zpencil
    
    integer(4) :: i, j, k
    type(DECOMP_INFO) :: dtmp

!-------------------------------------------------------------------------------
! x-pencil : u1 -> g1
!-------------------------------------------------------------------------------
    dtmp = dm%dpcc
    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
        fix(:) = fl%dDens(:, j, k)
        call Get_midp_interpolation_1D( 'x', 'C2P', dm, fix(:), fox(:) )
        fl%gx(:, j, k) = fox(:) * fl%qx(:, j, k)
      end do
    end do
!-------------------------------------------------------------------------------
! x-pencil --> y-pencil
!-------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qy,    uy_ypencil, dm%dcpc)
    call transpose_x_to_y(fl%dDens,  d_ypencil, dm%dccc)
    call transpose_x_to_y(fl%qz,    uz_ypencil, dm%dccp)
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
    call transpose_y_to_x(duz_ypencil, fl%gz, dm%dccp)
    call transpose_y_to_x(duy_ypencil, fl%gy, dm%dcpc)

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
!>
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
!------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    needed         specified   pubic
!-------------------------------------------------------------------------------
!> MPI : 
!>     default x-pencil
!>     working in : y-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!> Y: index arrangment
!>      j'-1   j'-1  j'    j'+1  j'+2
!>      _|__.__|__.__|__.__|__.__|__.__
!>         j-2   j-1   j     j+1    j+2
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[inout]         
!===============================================================================
  subroutine Get_volumetric_average_3d(var, dtmp, dm, fo_work)
    use mpi_mod
    use udf_type_mod
    implicit none
    type(t_domain),    intent(in) :: dm
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP),          intent(in) :: var(:, :, :)
    real(WP),         intent(out) :: fo_work
 
    real(WP), allocatable   :: varcp_ypencil(:, :, :)
    real(WP), allocatable   :: var_ypencil(:, :, :)
    real(WP)   :: vol, fo, vol_work
    integer(4) :: i, j, k, noy

!-------------------------------------------------------------------------------
!   transpose to y pencil. Default is x-pencil.
!-------------------------------------------------------------------------------
    allocate ( var_ypencil(dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3)) )
    var_ypencil = ZERO

    call transpose_x_to_y(var, var_ypencil, dtmp)
!-------------------------------------------------------------------------------
!   In Y-pencil now
!-------------------------------------------------------------------------------
    if( dtmp%ysz(2) == dm%np(2) )then
!-------------------------------------------------------------------------------
!   if variable is stored in y-nodes, extend them to y-cell centres
!   for example, uy
!-------------------------------------------------------------------------------
      if( dm%is_periodic(2) ) then
        noy = dtmp%ysz(2)
      else
        noy = dtmp%ysz(2) - 1
      end if

      allocate( varcp_ypencil(dtmp%ysz(1), noy, dtmp%ysz(3)) )
      varcp_ypencil = ZERO

      call Get_y_midp_P2C_3dArray ( var_ypencil, dm, varcp_ypencil)

      fo = ZERO
      vol = ZERO
      do k = 1, dtmp%ysz(3)
        do i = 1, dtmp%ysz(1)
          do j = 1, noy
!>       j'    j'+1
!>      _|__.__|_
!>         j     
            ! fo = fo + &
            !     ( dm%yp(j + 1) - dm%yp(j) ) / SIX * &
            !     ( fi3d_ypencil(i, j, k) + &
            !       FOUR * fo3dy(i, j, k) + &
            !       fi3d_ypencil(i, dm%jNeighb(3, j), k)) ! Simpson 2nd order 
            fo = fo + &      
                ( var_ypencil(i, dm%jNeighb(3, j), k) + varcp_ypencil(i, j, k) ) * &
                ( dm%yp(j + 1) - dm%yc(j) ) * HALF + &
                ( var_ypencil(i, j,                k) + varcp_ypencil(i, j, k) ) * &
                ( dm%yc(j    ) - dm%yp(j) ) * HALF
            vol = vol + ( dm%yp(j + 1) - dm%yp(j) )
          end do
        end do
      end do
      deallocate(fo3d)
    else
!-------------------------------------------------------------------------------
!   if variable is not stored in y-nodes, extends them to y-nodes.
!   for example, ux, density, etc.
!-------------------------------------------------------------------------------
      if( dm%is_periodic(2) ) then
        noy = dtmp%ysz(2)
      else
        noy = dtmp%ysz(2) + 1
      end if
      allocate( varcp_ypencil(dtmp%ysz(1), noy, dtmp%ysz(3)) )
      varcp_ypencil = ZERO
      call Get_y_midp_C2P_3dArray ( var_ypencil, dm, varcp_ypencil)

      fo = ZERO
      vol = ZERO
      do k = 1, dtmp%ysz(3)
        do i = 1, dtmp%ysz(1)
          do j = 1, dtmp%ysz(2)
!>       j'    j'+1
!>      _|__.__|_
!>         j  
            fo = fo + &
                ( varcp_ypencil(i, dm%jNeighb(3, j), k) + var_ypencil(i, j, k) ) * &
                ( dm%yp(j + 1) - dm%yc(j) ) * HALF + &
                ( varcp_ypencil(i, j,                k) + var_ypencil(i, j, k) ) * &
                ( dm%yc(j    ) - dm%yp(j) ) * HALF
            vol = vol + ( dm%yp(j + 1) - dm%yp(j) )
          end do
        end do
      end do
      deallocate(varcp_ypencil)
    end if
    deallocate(var_ypencil)
    
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

!===============================================================================
!>\brief : to find the maximum values of velocity
!------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    needed         N/A         pubic
!-------------------------------------------------------------------------------
!> MPI : 
!>     default x-pencil
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[inout]         
!===============================================================================
  subroutine Check_maximum_velocity(ux, uy, uz)
    use precision_mod
    use math_mod
    use mpi_mod
    implicit none

    real(WP), intent(in) :: ux(:, :, :), uy(:, :, :), uz(:, :, :)

    real(WP)   :: u(3), u_work(3)

    u(1) = MAXVAL( abs_wp( ux(:, :, :) ) )
    u(2) = MAXVAL( abs_wp( uy(:, :, :) ) )
    u(3) = MAXVAL( abs_wp( uz(:, :, :) ) )

    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(u, u_work, 3, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)

    if(nrank == 0) then
      Call Print_debug_mid_msg("  The maximum velocities are:")
      write (OUTPUT_UNIT, '(5X, A, 1ES13.5)') 'Umax : ', u_work(1)
      write (OUTPUT_UNIT, '(5X, A, 1ES13.5)') 'Umax : ', u_work(2)
      write (OUTPUT_UNIT, '(5X, A, 1ES13.5)') 'Umax : ', u_work(3)
    end if

    return
  end subroutine

end module
