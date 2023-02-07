module eq_energy_mod
  use operations
  use decomp_2d
  implicit none

  private :: Compute_energy_rhs
  private :: Calculate_energy_fractional_step
  public  :: Solve_energy_eq
  public  :: Calculate_massflux_from_velocity
  public  :: Calculate_velocity_from_massflux
contains
!==========================================================================================================
  subroutine Calculate_energy_fractional_step(rhs0, rhs1, dm, isub)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in) :: dm
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent(inout) :: rhs0, rhs1
    integer,  intent(in) :: isub
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: rhs_dummy


  ! add explicit terms
    rhs_dummy(:, :, :) = rhs1(:, :, :)
    rhs1(:, :, :) = dm%tGamma(isub) * rhs1(:, :, :) + &
                    dm%tZeta (isub) * rhs0(:, :, :)
    rhs0(:, :, :) = rhs_dummy(:, :, :)

  ! times the time step 
    rhs1(:, :, :) = dm%dt * rhs1(:, :, :)

    return
  end subroutine
!==========================================================================================================
  subroutine Compute_energy_rhs(fl, tm, dm, isub)
    use operations
    use udf_type_mod
    use thermo_info_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(in) :: fl
    type(t_thermo), intent(inout) :: tm
    integer,        intent(in) :: isub    

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: apcc
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil
    
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: gy_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: gz_zpencil 

    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dppc%xsz(3) ) :: hEnth_pcc
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dppc%xsz(3) ) :: kCond_pcc
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: hEnth_cpc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: hEnth_ccp_zpencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: Ttemp_ypencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: ene_rhs_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: kCond_cpc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: kCond_ccp_zpencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: Ttemp_zpencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: kCond_zpencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: ene_rhs_zpencil
    
    real(WP) :: fbcx(2), fbcy(2), fbcz(2)
    integer  :: ibcx(2), ibcy(2), ibcz(2)
    integer  :: i
!==========================================================================================================
!   preparation
!==========================================================================================================
    call transpose_x_to_y(fl%gy,        gy_ypencil,   dm%dcpc)   ! for d(g_y h)/dy
    call transpose_x_to_y(fl%gz,        accp_ypencil, dm%dccp)   ! intermediate, accp_ypencil = gz_ypencil
    call transpose_y_to_z(accp_ypencil, gz_zpencil,   dm%dccp)   ! for d(g_z h)/dz
!----------------------------------------------------------------------------------------------------------
!    h --> h_pcc
!      --> h_ypencil --> h_cpc_ypencil
!                    --> h_zpencil --> h_ccp_zpencil
!----------------------------------------------------------------------------------------------------------
    call Get_x_midp_C2P_3D(tm%hEnth,     hEnth_pcc,         dm, dm%ibcx(:, 5), dm%ftpbcx(:, :, :)%h ) ! for d(g_x h_pcc))/dy
    call transpose_x_to_y (tm%hEnth,     accc_ypencil, dm%dccc)                     !intermediate, accc_ypencil = hEnth_ypencil
    call Get_y_midp_C2P_3D(accc_ypencil, hEnth_cpc_ypencil, dm, dm%ibcy(:, 5), dm%ftpbcy(:, :, :)%h)! for d(g_y h_cpc)/dy
    call transpose_y_to_z (accc_ypencil, accc_zpencil, dm%dccc) !intermediate, accc_zpencil = hEnth_zpencil
    call Get_z_midp_C2P_3D(accc_zpencil, hEnth_ccp_zpencil, dm, dm%ibcz(:, 5), dm%ftpbcz(:, :, :)%h) ! for d(g_z h_ccp)/dz
!----------------------------------------------------------------------------------------------------------
!    k --> k_pcc
!      --> k_ypencil --> k_cpc_ypencil
!                    --> k_zpencil --> k_ccp_zpencil              
!----------------------------------------------------------------------------------------------------------
    do i = 1, 2
      fbcx(i) = tm%ftpbcx(i)%k
      fbcy(i) = tm%ftpbcy(i)%k
      fbcz(i) = tm%ftpbcz(i)%k
    end do
    call Get_x_midp_C2P_3D(tm%kCond,      kCond_pcc,         dm, dm%ibcx(:, 5), dm%ftpbcx(:, :, :)%k) ! for d(k_pcc * (dT/dx) )/dx
    call transpose_x_to_y (tm%kCond,      accc_ypencil, dm%dccc)  ! for k d2(T)/dy^2
    call Get_y_midp_C2P_3D(accc_ypencil,  kCond_cpc_ypencil, dm, dm%ibcy(:, 5), dm%ftpbcy(:, :, :)%k)
    call transpose_y_to_z (accc_ypencil,  kCond_zpencil, dm%dccc) 
    call Get_z_midp_C2P_3D(kCond_zpencil, kCond_ccp_zpencil, dm, dm%ibcz(:, 5), dm%ftpbcz(:, :, :)%k)
!----------------------------------------------------------------------------------------------------------
!    T --> T_ypencil --> T_zpencil
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y (tm%Ttemp,      Ttemp_ypencil, dm%dccc)   ! for k d2(T)/dy^2
    call transpose_y_to_z (Ttemp_ypencil, Ttemp_zpencil, dm%dccc)   ! for k d2(T)/dz^2
!==========================================================================================================
! the RHS of energy equation
! x-pencil : the RHS terms of energy (derivative) operating in the x direction
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
! x-pencil : d (gx * h_pcc) / dx 
!----------------------------------------------------------------------------------------------------------
    tm%ene_rhs = ZERO
    call Get_x_1st_derivative_P2C_3D( - fl%gx * hEnth_pcc, accc, dm, dm%ibcx(i, 5), dm%ftpbcx(:, :, :)%dh * dm%fbcx(:, :, :, 1 + NBC) ) ! accc = -d(gx * h)/dx
    tm%ene_rhs = tm%ene_rhs + accc
!----------------------------------------------------------------------------------------------------------
! x-pencil : d (T) / dx 
!----------------------------------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3D(tm%tTemp, apcc, dm, dm%ibcx(i, 5), dm%ftpbcx(:, :, :)%t )
!----------------------------------------------------------------------------------------------------------
! x-pencil : k_pcc * d (T) / dx 
!----------------------------------------------------------------------------------------------------------
    apcc = apcc * kCond_pcc
    do n = 1, 2
      if (dm%ibcx(1, 5) == IBC_NEUMANN) then
        fbcx(:, :, :) = dm%fbcx_var(:, :, :, 5)
      else if (dm%ibcx(1, 5) == IBC_DIRICHLET) then
        fbcx(:, :, :) = dm%ftpbcx(:, :, :)%t * dm%ftpbcx(:, :, :)%k
      end if
    end if
!----------------------------------------------------------------------------------------------------------
! x-pencil : d ( k_pcc * d (T) / dx ) dx
!----------------------------------------------------------------------------------------------------------
    call Get_x_1st_derivative_P2C_3D(apcc, accc, dm, dm%ibcx(:, 5), fbcx )

    tm%ene_rhs = tm%ene_rhs + accc
!==========================================================================================================
! the RHS of energy equation
! y-pencil : the RHS terms of energy (derivative) operating in the y direction
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
! y-pencil : d (gy * h_cpc) / dy 
!----------------------------------------------------------------------------------------------------------
    ene_rhs_ypencil = ZERO
    call Get_y_1st_derivative_P2C_3D( - gy_ypencil * hEnth_cpc_ypencil, accc_ypencil, dm, dm%ibcy(:, 2), dm%ftpbcy(:, :, :)%dh * dm%fbcy(:, :, :, 2 + NBC) )
    ene_rhs_ypencil = ene_rhs_ypencil + accc_ypencil
!----------------------------------------------------------------------------------------------------------
! y-pencil : d (T) / dy
!----------------------------------------------------------------------------------------------------------
    call Get_y_1st_derivative_C2P_3D(tTemp_ypencil, acpc_ypencil, dm, dm%ibcy(i, 5), dm%ftpbcy(:, :, :)%t )
!----------------------------------------------------------------------------------------------------------
! y-pencil : k_cpc * d (T) / dy 
!----------------------------------------------------------------------------------------------------------
    acpc_ypencil = acpc_ypencil * kCond_cpc_ypencil
    do n = 1, 2
      if (dm%ibcy(1, 5) == IBC_NEUMANN) then
        fbcy(:, :, :) = dm%fbcy_var(:, :, :, 5)
      else if (dm%ibcy(1, 5) == IBC_DIRICHLET) then
        fbcy(:, :, :) = dm%ftpbcy(:, :, :)%t * dm%ftpbcy(:, :, :)%k
      end if
    end if
!----------------------------------------------------------------------------------------------------------
! y-pencil : d ( k_cpc * d (T) / dy ) dy
!----------------------------------------------------------------------------------------------------------
    call Get_y_1st_derivative_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%ibcy(:, 5), fbcy )
    ene_rhs_ypencil = ene_rhs_ypencil + accc_ypencil

    call transpose_y_to_x(ene_rhs_ypencil, accc, dm%dccc)
    tm%ene_rhs = tm%ene_rhs + accc
!==========================================================================================================
! the RHS of energy equation
! z-pencil : the RHS terms of energy (derivative) operating in the z direction
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
! z-pencil : d (gz * h_ccp) / dz 
!----------------------------------------------------------------------------------------------------------
    ene_rhs_zpencil = ZERO
    call Get_z_1st_derivative_P2C_3D( - gz_zpencil * hEnth_ccp_zpencil, accc_zpencil, dm, dm%ibcz(:, 3), dm%ftpbcz(:, :, :)%dh * dm%fbcz(:, :, :, 3 + NBC )
    ene_rhs_zpencil = ene_rhs_zpencil + accc_zpencil
!----------------------------------------------------------------------------------------------------------
! z-pencil : d (T) / dz
!----------------------------------------------------------------------------------------------------------
    call Get_z_1st_derivative_C2P_3D(tTemp_zpencil, accp_zpencil, dm, dm%ibcz(i, 5), dm%ftpbcz(:, :, :)%t )
!----------------------------------------------------------------------------------------------------------
! z-pencil : k_ccp * d (T) / dz 
!----------------------------------------------------------------------------------------------------------
    accp_zpencil = accp_zpencil * kCond_ccp_zpencil
    do n = 1, 2
      if (dm%ibcz(1, 5) == IBC_NEUMANN) then
        fbcz(:, :, :) = dm%fbcz_var(:, :, :, 5)
      else if (dm%ibcy(1, 5) == IBC_DIRICHLET) then
        fbcz(:, :, :) = dm%ftpbcz(:, :, :)%t * dm%ftpbcz(:, :, :)%k
      end if
    end if
!----------------------------------------------------------------------------------------------------------
! z-pencil : d ( k_ccp * d (T) / dz ) / dz
!----------------------------------------------------------------------------------------------------------
    call Get_z_1st_derivative_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%ibcz(:, 5), fbcz )
    ene_rhs_zpencil = ene_rhs_zpencil + accc_zpencil

    call transpose_z_to_y(ene_rhs_zpencil, ene_rhs_ypencil, dm%dccc)
    call transpose_y_to_x(ene_rhs_ypencil, accc,            dm%dccc)
    tm%ene_rhs = tm%ene_rhs + accc

!==========================================================================================================
! time approaching
!==========================================================================================================
    call Calculate_energy_fractional_step(tm%ene_rhs0, tm%ene_rhs, dm, isub)

    return
  end subroutine Compute_energy_rhs
!==========================================================================================================
!==========================================================================================================
  subroutine Solve_energy_eq(fl, tm, dm, isub)
    use udf_type_mod
    use thermo_info_mod 
    use solver_tools_mod
    implicit none
    type(t_domain), intent(in)    :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm
    integer,        intent(in)    :: isub

!----------------------------------------------------------------------------------------------------------
!   calculate rhs of energy equation
!----------------------------------------------------------------------------------------------------------
    call Compute_energy_rhs(fl, tm, dm, isub)
!----------------------------------------------------------------------------------------------------------
!   update rho * h
!----------------------------------------------------------------------------------------------------------
    tm%dh = tm%dh + tm%ene_rhs
!----------------------------------------------------------------------------------------------------------
!   update other properties from rho * h
!----------------------------------------------------------------------------------------------------------
    call Update_thermal_properties(fl, tm, dm)
!----------------------------------------------------------------------------------------------------------
!   No Need to apply b.c.
!----------------------------------------------------------------------------------------------------------
  return
  end subroutine

!==========================================================================================================
!> \brief Calculate the conservative variables from primary variable.     
!> This subroutine is called to update $\rho u_i$ from $u_i$.
!---------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    needed         specified   pubic
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     dm             domain
!> \param[in]     fm             flow
!==========================================================================================================
  subroutine Calculate_massflux_from_velocity(fl, dm)
    use udf_type_mod
    use operations
    use decomp_2d
    implicit none
    type(t_domain), intent(in )   :: dm
    type(t_flow  ), intent(inout) :: fl

    real(WP), dimension( dm%dpcc%xsz(1), &
                         dm%dpcc%xsz(2), &
                         dm%dpcc%xsz(3)) :: d_pcc
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: d_cpc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: d_ccp_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: qy_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: gy_ypencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) :: qz_ypencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) :: gz_ypencil
    real(WP), dimension( dm%dccc%ysz(1), &
                         dm%dccc%ysz(2), &
                         dm%dccc%ysz(3)) ::  d_ypencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: qz_zpencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: gz_zpencil
    real(WP), dimension( dm%dccc%zsz(1), &
                         dm%dccc%zsz(2), &
                         dm%dccc%zsz(3)) ::  d_zpencil
    
!----------------------------------------------------------------------------------------------------------
! x-pencil : u1 -> g1
!----------------------------------------------------------------------------------------------------------
    call Get_x_midp_C2P_1D (fl%dDens, d_pcc, dm, dm%ibcx(:, 5), dm%ftpbcx(:, :, :)%d)
    fl%gx = fl%qx * d_pcc
!----------------------------------------------------------------------------------------------------------
! y-pencil : u2 -> g2
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qy,    qy_ypencil, dm%dcpc)
    call transpose_x_to_y(fl%dDens,  d_ypencil, dm%dccc)

    call Get_y_midp_C2P_1D (d_ypencil, d_cpc_ypencil, dm, dm%ibcy(:, 5), dm%ftpbcy(:, :, :)%d)
    gy_ypencil = qy_ypencil * d_cpc_ypencil
    call transpose_y_to_x(gy_ypencil, fl%gy, dm%dcpc)
!----------------------------------------------------------------------------------------------------------
! Z-pencil : u3 -> g3
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_z( d_ypencil,  d_zpencil, dm%dccc)
    call transpose_x_to_y(fl%qz,      qz_ypencil, dm%dccp)
    call transpose_y_to_z(qz_ypencil, qz_zpencil, dm%dccp)

    call Get_z_midp_C2P_1D (d_zpencil, d_ccp_zpencil, dm, dm%ibcz(:, 5), dm%ftpbcz(:, :, :)%d)
    gz_zpencil = qz_zpencil * d_ccp_zpencil

    call transpose_z_to_y(gz_zpencil, gz_ypencil, dm%dccp)
    call transpose_y_to_x(gz_ypencil, fl%gz,      dm%dccp)

    return
  end subroutine Calculate_massflux_from_velocity


  !==========================================================================================================
  subroutine Calculate_velocity_from_massflux(fl, dm)
    use udf_type_mod
    use operations
    use decomp_2d
    implicit none
    type(t_domain), intent(in )   :: dm
    type(t_flow  ), intent(inout) :: fl

    real(WP), dimension( dm%dpcc%xsz(1), &
                         dm%dpcc%xsz(2), &
                         dm%dpcc%xsz(3)) :: d_pcc
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: d_cpc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: d_ccp_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: qy_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: gy_ypencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) :: qz_ypencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) :: gz_ypencil
    real(WP), dimension( dm%dccc%ysz(1), &
                         dm%dccc%ysz(2), &
                         dm%dccc%ysz(3)) ::  d_ypencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: qz_zpencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: gz_zpencil
    real(WP), dimension( dm%dccc%zsz(1), &
                         dm%dccc%zsz(2), &
                         dm%dccc%zsz(3)) ::  d_zpencil
    
!----------------------------------------------------------------------------------------------------------
! x-pencil : g1 -> u1
!----------------------------------------------------------------------------------------------------------
    call Get_x_midp_C2P_1D (fl%dDens, d_pcc, dm, dm%ibcx(:, 5), dm%ftpbcx(:, :, :)%d)
    fl%qx = fl%gx / d_pcc
!----------------------------------------------------------------------------------------------------------
! y-pencil : u2 -> g2
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%gy,    gy_ypencil, dm%dcpc)
    call transpose_x_to_y(fl%dDens,  d_ypencil, dm%dccc)

    call Get_y_midp_C2P_1D (d_ypencil, d_cpc_ypencil, dm, dm%ibcy(:, 5), dm%ftpbcy(:, :, :)%d)
    qy_ypencil = gy_ypencil / d_cpc_ypencil
    call transpose_y_to_x(qy_ypencil, fl%qy, dm%dcpc)
!----------------------------------------------------------------------------------------------------------
! Z-pencil : u3 -> g3
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_z( d_ypencil,  d_zpencil, dm%dccc)
    call transpose_x_to_y(fl%gz,      gz_ypencil, dm%dccp)
    call transpose_y_to_z(gz_ypencil, gz_zpencil, dm%dccp)

    call Get_z_midp_C2P_1D (d_zpencil, d_ccp_zpencil, dm, dm%ibcz(:, 5), dm%ftpbcz(:, :, :)%d)
    qz_zpencil = gz_zpencil / d_ccp_zpencil

    call transpose_z_to_y(qz_zpencil, qz_ypencil, dm%dccp)
    call transpose_y_to_x(qz_ypencil, fl%qz, dm%dccp)

    return
  end subroutine Calculate_velocity_from_massflux

end module eq_energy_mod
