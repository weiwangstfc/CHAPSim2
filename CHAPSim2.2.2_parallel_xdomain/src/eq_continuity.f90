module continuity_eq_mod
  use operations
  use decomp_2d
  private :: Calculate_drhodt
  private :: Get_divergence_vel
  public  :: Calculate_continuity_constrains
  public  :: Check_mass_conservation
contains
!=============================================================================================================================================
!=============================================================================================================================================
!> \brief To calculate d(\rho)/dt in the continuity eq.
!---------------------------------------------------------------------------------------------------------------------------------------------
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
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent ( in  ) :: dDens, dDensm1, dDensm2
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent ( out ) :: drhodt

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

!=============================================================================================================================================
!=============================================================================================================================================
!> \brief To calculate divergence of (rho * u) or divergence of (u)
!---------------------------------------------------------------------------------------------------------------------------------------------
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
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent (in   ) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent (in   ) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent (in   ) :: uz
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent (inout) :: div

    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div0
    real(WP), dimension(dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3)) :: uy_ypencil
    real(WP), dimension(dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3)) :: uz_ypencil
    real(WP), dimension(dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3)) :: uz_zpencil
    real(WP), dimension(dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3)) :: div0_ypencil
    real(WP), dimension(dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3)) :: div0_zpencil

    div(:, :, :) = ZERO
!---------------------------------------------------------------------------------------------------------------------------------------------
! operation in x pencil, du/dx
!---------------------------------------------------------------------------------------------------------------------------------------------
    div0(:, :, :) = ZERO
    call Get_x_1st_derivative_P2C_3D(ux, div0, dm, dm%ibcx(:, 1))
    div(:, :, :) = div(:, :, :) + div0(:, :, :)
!---------------------------------------------------------------------------------------------------------------------------------------------
! operation in y pencil, dv/dy
!---------------------------------------------------------------------------------------------------------------------------------------------
    div0(:, :, :) = ZERO
    call transpose_x_to_y(uy, uy_ypencil, dm%dcpc)
    call Get_y_1st_derivative_P2C_3D(uy_ypencil, div0_ypencil, dm, dm%ibcy(:, 2))
    call transpose_y_to_x(div0_ypencil, div0, dm%dccc)
    div(:, :, :) = div(:, :, :) + div0(:, :, :)
!---------------------------------------------------------------------------------------------------------------------------------------------
! operation in z pencil, dv/dz
!---------------------------------------------------------------------------------------------------------------------------------------------
    div0(:, :, :) = ZERO
    call transpose_x_to_y(uz,         uz_ypencil, dm%dccp)
    call transpose_y_to_z(uz_ypencil, uz_zpencil, dm%dccp)
    call Get_z_1st_derivative_P2C_3D(uz_zpencil, div0_zpencil, dm, dm%ibcz(:, 3))
    call transpose_z_to_y(div0_zpencil, div0_ypencil, dm%dccp)
    call transpose_y_to_x(div0_ypencil, div0,         dm%dccp)
    div(:, :, :) = div(:, :, :) + div0(:, :, :)

    return
  end subroutine

!=============================================================================================================================================
!=============================================================================================================================================
!> \brief To calculate divergence of (rho * u) or divergence of (u)
!---------------------------------------------------------------------------------------------------------------------------------------------
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
    use solver_tools_mod
    use wtformat_mod
    implicit none

    type(t_domain), intent( in    ) :: dm
    type(t_flow),   intent( inout ) :: fl                  

    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div
    real(WP)   :: divmax

    fl%pcor = ZERO
    div(:, :, :)  = ZERO
!---------------------------------------------------------------------------------------------------------------------------------------------
! $d\rho / dt$ at cell centre
!_______________________________________________________________________________
    if (dm%ithermo == 1) then
      call Calculate_drhodt(dm, fl%dDens, fl%dDensm1, fl%dDensm2, fl%pcor)
    end if
!---------------------------------------------------------------------------------------------------------------------------------------------
! $d(\rho u_i)) / dx_i $ at cell centre
!_______________________________________________________________________________
    if (dm%ithermo == 1) then
      call Get_divergence_vel(fl%gx, fl%gy, fl%gz, div, dm)
    else
      call Get_divergence_vel(fl%qx, fl%qy, fl%qz, div, dm)
    end if
  
    call Find_maximum_absvar3d(div, divmax)

    if(nrank == 0) then
      write (OUTPUT_UNIT, wrtfmt1e) "  Check Mass Conservation:", divmax
    end if

    return
  end subroutine Check_mass_conservation


!=============================================================================================================================================
!=============================================================================================================================================
  subroutine Calculate_continuity_constrains(fl, dm, isub)
    use parameters_constant_mod
    use udf_type_mod
    implicit none

    type(t_domain), intent( in    ) :: dm
    type(t_flow),   intent( inout ) :: fl                  
    integer,     intent( in    ) :: isub

    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div

    fl%pcor = ZERO
    div  = ZERO
!---------------------------------------------------------------------------------------------------------------------------------------------
! $d\rho / dt$ at cell centre
!_______________________________________________________________________________
    if (dm%ithermo == 1) then
      call Calculate_drhodt(dm, fl%dDens, fl%dDensm1, fl%dDensm2, fl%pcor)
    end if
!---------------------------------------------------------------------------------------------------------------------------------------------
! $d(\rho u_i)) / dx_i $ at cell centre
!_______________________________________________________________________________
    if (dm%ithermo == 1) then
      call Get_divergence_vel(fl%gx, fl%gy, fl%gz, div, dm)
    else
      call Get_divergence_vel(fl%qx, fl%qy, fl%qz, div, dm)
    end if

    fl%pcor = fl%pcor + div
    fl%pcor = fl%pcor / (dm%tAlpha(isub) * dm%sigma2p * dm%dt)

    return
  end subroutine Calculate_continuity_constrains

end module continuity_eq_mod
