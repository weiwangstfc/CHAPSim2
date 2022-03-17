module continuity_eq_mod
  use operations
  use decomp_2d
  private :: Calculate_drhodt
  private :: Get_divergence_vel
  public  :: Calculate_continuity_constrains
  public  :: Check_mass_conservation
contains
!===============================================================================
!===============================================================================
!> \brief To calculate d(\rho)/dt in the continuity eq.
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  dDens            density at the current time step
!> \param[in]  dDensm1          density at the t-1 time step
!> \param[in]  dDensm2          density at the t-2 time step
!> \param[out] drhodt           d(rho)/dt
!> \param[in]  itime            the sub-step in RK3
!_______________________________________________________________________________
  subroutine Calculate_drhodt(dm, dDens, dDensm1, dDensm2, drhodt)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in) :: dm
    real(WP), dimension(:, :, :), intent ( in  ) :: dDens, dDensm1, dDensm2
    real(WP), dimension(:, :, :), intent ( out ) :: drhodt

    integer :: i

    if(dm%iTimeScheme == ITIME_AB2) then

      drhodt(:, :, :) = HALF * dDens  (:, :, :) - &
                        TWO  * dDensm1(:, :, :) + &
                        HALF * dDensm2(:, :, :)
      drhodt(:, :, :) = drhodt(:, :, :) * dm%dt

    else if (dm%iTimeScheme == ITIME_RK3 .or. dm%iTimeScheme == ITIME_RK3_CN) then

      ! to check this part, is iteration necessary?
      drhodt(:, :, :) = dDens  (:, :, :)
      do i = 1, dm%nsubitr
        drhodt(:, :, :) = drhodt(:, :, :) + dm%tAlpha(i) * &
                          (dDensm1(:, :, :) - dDensm2(:, :, :))  * dm%dt
      end do

    else  

      ! default, Euler 1st order 
      drhodt(:, :, :) = dDens(:, :, :) - dDensm1(:, :, :)
      drhodt(:, :, :) = drhodt(:, :, :) * dm%dt

    end if


    return
  end subroutine Calculate_drhodt

!===============================================================================
!===============================================================================
!> \brief To calculate divergence of (rho * u) or divergence of (u)
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ux           ux or gx
!> \param[in]     uy           uy or gy
!> \param[in]     uz           uz or gz
!> \param[out]    div          div(u) or div(g)
!> \param[in]     d            domain
!_______________________________________________________________________________
  subroutine Get_divergence_vel(ux, uy, uz, div, dm)
    use parameters_constant_mod
    use udf_type_mod
    implicit none

    type(t_domain),               intent (in   ) :: dm
    real(WP), dimension(:, :, :), intent (in   ) :: ux
    real(WP), dimension(:, :, :), intent (in   ) :: uy
    real(WP), dimension(:, :, :), intent (in   ) :: uz
    real(WP), dimension(:, :, :), intent (inout) :: div

    real(WP), allocatable :: div0(:, :, :)
    integer :: nx, ny, nz

    nx = dm%nc(1)
    ny = dm%nc(2)
    nz = dm%nc(3)

    allocate( div0(nx, ny, nz) )

    div(:, :, :) = ZERO
!-------------------------------------------------------------------------------
! operation in x pencil, du/dx
!_______________________________________________________________________________
    !call Print_3d_array(ux, nx, ny, nz, 'ux:') ! test
    call Get_x_1st_derivative_P2C_3D(ux, div0, dm, dm%ibcx(:, 1))
    div(:, :, :) = div(:, :, :) + div0(:, :, :)
    !call Print_3d_array(div0, nx, ny, nz, 'du/dx:') ! test
!-------------------------------------------------------------------------------
! operation in y pencil, dv/dy
!_______________________________________________________________________________
    call Get_y_1st_derivative_P2C_3D(uy, div0, dm, dm%ibcy(:, 2))
    div(:, :, :) = div(:, :, :) + div0(:, :, :)
    !call Print_3d_array(div0, nx, ny, nz, 'dv/dy:')

!-------------------------------------------------------------------------------
! operation in z pencil, dv/dz
!_______________________________________________________________________________
    call Get_z_1st_derivative_P2C_3D(uz, div0, dm, dm%ibcz(:, 3))
    div(:, :, :) = div(:, :, :) + div0(:, :, :)
    !call Print_3d_array(div0, nx, ny, nz, 'dw/dz:')


    deallocate(div0)
    return
  end subroutine

!===============================================================================
!===============================================================================
!> \brief To calculate divergence of (rho * u) or divergence of (u)
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ux           ux or gx
!> \param[in]     uy           uy or gy
!> \param[in]     uz           uz or gz
!> \param[out]    div          div(u) or div(g)
!> \param[in]     d            domain
!_______________________________________________________________________________
  subroutine Check_mass_conservation(fl, dm)
    use precision_mod
    use udf_type_mod
    use input_general_mod    
    use parameters_constant_mod
    use math_mod                
    use mpi_mod
    implicit none

    type(t_domain), intent( in    ) :: dm
    type(t_flow),   intent( inout ) :: fl                  

    real(WP), allocatable  :: div (:, :, :)
    integer :: nx, ny, nz
    integer :: loc3d(3)
    real(WP)   :: divmax

    nx = size(fl%pcor, 1)
    ny = size(fl%pcor, 2)
    nz = size(fl%pcor, 3)
    allocate( div(nx, ny, nz) )

    fl%pcor = ZERO
    div  = ZERO
!-------------------------------------------------------------------------------
! $d\rho / dt$ at cell centre
!_______________________________________________________________________________
    if (dm%ithermo == 1) then
      call Calculate_drhodt(dm, fl%dDens, fl%dDensm1, fl%dDensm2, fl%pcor)
    end if
!-------------------------------------------------------------------------------
! $d(\rho u_i)) / dx_i $ at cell centre
!_______________________________________________________________________________
    if (dm%ithermo == 1) then
      call Get_divergence_vel(fl%gx, fl%gy, fl%gz, div, dm)
    else
      call Get_divergence_vel(fl%qx, fl%qy, fl%qz, div, dm)
    end if
  
    divmax = MAXVAL( abs_wp( div ) )
    loc3d  = MAXLOC( abs_wp( div ) )

    if(nrank == 0) then
      call Print_debug_mid_msg("  The maximum value of mass conservation and loc are:")
      write (OUTPUT_UNIT, "(12X, 3I8.1, 1ES13.5)") loc3d, divmax
    end if

    deallocate (div)

    return
  end subroutine Check_mass_conservation


!===============================================================================
!===============================================================================
  subroutine Calculate_continuity_constrains(fl, dm, isub)
    use parameters_constant_mod
    use udf_type_mod
    implicit none

    type(t_domain), intent( in    ) :: dm
    type(t_flow),   intent( inout ) :: fl                  
    integer,     intent( in    ) :: isub

    real(WP), allocatable  :: div (:, :, :)
    integer :: nx, ny, nz

    nx = dm%nc(1)
    ny = dm%nc(2)
    nz = dm%nc(3)
    allocate( div(nx, ny, nz) )

    fl%pcor = ZERO
    div  = ZERO
!-------------------------------------------------------------------------------
! $d\rho / dt$ at cell centre
!_______________________________________________________________________________
    if (dm%ithermo == 1) then
      call Calculate_drhodt(dm, fl%dDens, fl%dDensm1, fl%dDensm2, fl%pcor)
    end if
!-------------------------------------------------------------------------------
! $d(\rho u_i)) / dx_i $ at cell centre
!_______________________________________________________________________________
    if (dm%ithermo == 1) then
      call Get_divergence_vel(fl%gx, fl%gy, fl%gz, div, dm)
    else
      call Get_divergence_vel(fl%qx, fl%qy, fl%qz, div, dm)
    end if

    fl%pcor = fl%pcor + div
    fl%pcor = fl%pcor / (dm%tAlpha(isub) * dm%sigma2p * dm%dt)
    
    deallocate (div)

    return
  end subroutine Calculate_continuity_constrains

end module continuity_eq_mod
