module continuity_eq_mod

  private :: Calculate_drhodt
  public  :: Get_divergence_vel
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
  subroutine Calculate_drhodt(dDens, dDensm1, dDensm2, drhodt)
    use precision_mod
    use udf_type_mod, only : t_flow, t_domain
    use input_general_mod, only : dt, iTimeScheme, ITIME_RK3, ITIME_RK3_CN, ITIME_AB2, &
         nsubitr, tAlpha
    use parameters_constant_mod, only : ONEPFIVE, TWO, HALF
    implicit none

    real(WP), dimension(:, :, :), intent ( in  ) :: dDens, dDensm1, dDensm2
    real(WP), dimension(:, :, :), intent ( out ) :: drhodt

    integer(4) :: i

    if(iTimeScheme == ITIME_AB2) then

      drhodt(:, :, :) = HALF * dDens  (:, :, :) - &
                        TWO  * dDensm1(:, :, :) + &
                        HALF * dDensm2(:, :, :)
      drhodt(:, :, :) = drhodt(:, :, :) * dt

    else if (iTimeScheme == ITIME_RK3 .or. iTimeScheme == ITIME_RK3_CN) then

      ! to check this part, is iteration necessary?
      drhodt(:, :, :) = dDens  (:, :, :)
      do i = 1, nsubitr
        drhodt(:, :, :) = drhodt(:, :, :) + tAlpha(i) * &
                          (dDensm1(:, :, :) - dDensm2(:, :, :))  * dt
      end do

    else  

      ! default, Euler 1st order 
      drhodt(:, :, :) = dDens(:, :, :) - dDensm1(:, :, :)
      drhodt(:, :, :) = drhodt(:, :, :) * dt

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
  subroutine Get_divergence_vel(ux, uy, uz, div, d)
    use precision_mod
    use udf_type_mod, only : t_domain, t_flow
    use parameters_constant_mod, only : ZERO
    use operations
    implicit none

    type(t_domain),               intent (in   ) :: d
    real(WP), dimension(:, :, :), intent (in   ) :: ux
    real(WP), dimension(:, :, :), intent (in   ) :: uy
    real(WP), dimension(:, :, :), intent (in   ) :: uz
    real(WP), dimension(:, :, :), intent (inout) :: div

    real(WP), allocatable :: div0(:, :, :)
    integer(4) :: nx, ny, nz

    nx = size(div, 1)
    ny = size(div, 2)
    nz = size(div, 3)

    allocate( div0(nx, ny, nz) )

    div(:, :, :) = ZERO
!-------------------------------------------------------------------------------
! operation in x pencil
!_______________________________________________________________________________
    call Get_x_1st_derivative_P2C_3dArray( d, ux, div0)
    div(:, :, :) = div(:, :, :) + div0(:, :, :)

!-------------------------------------------------------------------------------
! operation in y pencil
!_______________________________________________________________________________
    call Get_y_1st_derivative_P2C_3dArray( d, uy, div0)
    div(:, :, :) = div(:, :, :) + div0(:, :, :)

!-------------------------------------------------------------------------------
! operation in z pencil
!_______________________________________________________________________________
    call Get_z_1st_derivative_P2C_3dArray( d, uz, div0)
    div(:, :, :) = div(:, :, :) + div0(:, :, :)


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
  subroutine Check_mass_conservation(f, d)
    use precision_mod,           only : WP
    use udf_type_mod,            only : t_domain, t_flow
    use input_general_mod,       only : ithermo
    use parameters_constant_mod, only : ZERO
    use math_mod,                only : abs_wp
    implicit none

    type(t_domain), intent( in    ) :: d
    type(t_flow),   intent( inout ) :: f                  

    real(WP), allocatable  :: div (:, :, :)
    integer(4) :: nx, ny, nz
    integer(4) :: loc3d(3)
    real(WP)   :: divmax

    nx = size(f%pcor, 1)
    ny = size(f%pcor, 2)
    nz = size(f%pcor, 3)
    allocate( div(nx, ny, nz) )

    f%pcor = ZERO
    div  = ZERO
!-------------------------------------------------------------------------------
! $d\rho / dt$ at cell centre
!_______________________________________________________________________________
    if (ithermo == 1) then
      call Calculate_drhodt(f%dDens, f%dDensm1, f%dDensm2, f%pcor)
    end if
!-------------------------------------------------------------------------------
! $d(\rho u_i)) / dx_i $ at cell centre
!_______________________________________________________________________________
    if (ithermo == 1) then
      call Get_divergence_vel(f%gx, f%gy, f%gz, div, d)
    else
      call Get_divergence_vel(f%qx, f%qy, f%qz, div, d)
    end if
  
    divmax = MAXVAL( abs_wp( div ) )
    loc3d  = MAXLOC( abs_wp( div ) )

    call Print_debug_mid_msg("  The maximum value of mass conservation and loc are:")
    write(*, "(12X, 3I8.1, 1ES13.5)") loc3d, divmax
    deallocate (div)

    return
  end subroutine Check_mass_conservation


!===============================================================================
!===============================================================================
  subroutine Calculate_continuity_constrains(f, d, isub)
    use precision_mod,           only : WP
    use udf_type_mod,            only : t_domain, t_flow
    use input_general_mod,       only : ithermo, dt, sigma2p, tAlpha
    use parameters_constant_mod, only : ZERO
    use math_mod,                only : abs_wp
    implicit none

    type(t_domain), intent( in    ) :: d
    type(t_flow),   intent( inout ) :: f                  
    integer(4),     intent( in    ) :: isub

    real(WP), allocatable  :: div (:, :, :)
    integer(4) :: nx, ny, nz

    nx = size(f%pcor, 1)
    ny = size(f%pcor, 2)
    nz = size(f%pcor, 3)
    allocate( div(nx, ny, nz) )

    f%pcor = ZERO
    div  = ZERO
!-------------------------------------------------------------------------------
! $d\rho / dt$ at cell centre
!_______________________________________________________________________________
    if (ithermo == 1) then
      call Calculate_drhodt(f%dDens, f%dDensm1, f%dDensm2, f%pcor)
    end if
!-------------------------------------------------------------------------------
! $d(\rho u_i)) / dx_i $ at cell centre
!_______________________________________________________________________________
    if (ithermo == 1) then
      call Get_divergence_vel(f%gx, f%gy, f%gz, div, d)
    else
      call Get_divergence_vel(f%qx, f%qy, f%qz, div, d)
    end if
    f%pcor = f%pcor + div
    f%pcor = f%pcor / (tAlpha(isub) * sigma2p * dt)
    deallocate (div)

    return
  end subroutine Calculate_continuity_constrains

end module continuity_eq_mod
