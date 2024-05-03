module eq_energy_mod
  use operations
  use decomp_2d
  implicit none

  private :: Compute_energy_rhs
  private :: Calculate_energy_fractional_step
  public  :: Solve_energy_eq
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
    
    real(WP) :: fbcx(4, dm%np(2), dm%np(3))
    real(WP) :: fbcy(dm%np(1), 4, dm%np(3))
    real(WP) :: fbcz(dm%np(1), dm%np(2), 4)
    integer  :: n
    integer  :: nbc(1:2, 1:3)
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
    call Get_x_midp_C2P_3D(tm%hEnth,     hEnth_pcc,         dm, dm%ibcx_Th(:, IBC_CCC), dm%ftpbcx_var(:, :, :)%h ) ! for d(g_x h_pcc))/dy
    call transpose_x_to_y (tm%hEnth,     accc_ypencil, dm%dccc)                     !intermediate, accc_ypencil = hEnth_ypencil
    call Get_y_midp_C2P_3D(accc_ypencil, hEnth_cpc_ypencil, dm, dm%ibcy_Th(:, IBC_CCC), dm%ftpbcy_var(:, :, :)%h)! for d(g_y h_cpc)/dy
    call transpose_y_to_z (accc_ypencil, accc_zpencil, dm%dccc) !intermediate, accc_zpencil = hEnth_zpencil
    call Get_z_midp_C2P_3D(accc_zpencil, hEnth_ccp_zpencil, dm, dm%ibcz_Th(:, IBC_CCC), dm%ftpbcz_var(:, :, :)%h) ! for d(g_z h_ccp)/dz
!----------------------------------------------------------------------------------------------------------
!    k --> k_pcc
!      --> k_ypencil --> k_cpc_ypencil
!                    --> k_zpencil --> k_ccp_zpencil              
!----------------------------------------------------------------------------------------------------------
    call Get_x_midp_C2P_3D(tm%kCond,      kCond_pcc,         dm, dm%ibcx_Th(:, IBC_CCC), dm%ftpbcx_var(:, :, :)%k) ! for d(k_pcc * (dT/dx) )/dx
    call transpose_x_to_y (tm%kCond,      accc_ypencil, dm%dccc)  ! for k d2(T)/dy^2
    call Get_y_midp_C2P_3D(accc_ypencil,  kCond_cpc_ypencil, dm, dm%ibcy_Th(:, IBC_CCC), dm%ftpbcy_var(:, :, :)%k)
    call transpose_y_to_z (accc_ypencil,  kCond_zpencil, dm%dccc) 
    call Get_z_midp_C2P_3D(kCond_zpencil, kCond_ccp_zpencil, dm, dm%ibcz_Th(:, IBC_CCC), dm%ftpbcz_var(:, :, :)%k)
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
    call get_ibc_for_calcuation(dm%ibcx_qx(:, IBC_PCC), nbc, dm%ibcx_Th(:, IBC_PCC))
    call Get_x_1st_derivative_P2C_3D( - fl%gx * hEnth_pcc, accc, dm, nbc(:, 1)) ! accc = -d(gx * h)/dx
    tm%ene_rhs = tm%ene_rhs + accc
!----------------------------------------------------------------------------------------------------------
! x-pencil : d (T) / dx 
!----------------------------------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3D(tm%tTemp, apcc, dm, dm%ibcx_Th(:, IBC_CCC), dm%ftpbcx_var(:, :, :)%t )
!----------------------------------------------------------------------------------------------------------
! x-pencil : k_pcc * d (T) / dx 
!----------------------------------------------------------------------------------------------------------
    apcc = apcc * kCond_pcc
!----------------------------------------------------------------------------------------------------------
! x-pencil : d ( k_pcc * d (T) / dx ) dx
!----------------------------------------------------------------------------------------------------------
    call get_ibc_for_calcuation(dm%ibcx_Th(:, IBC_PCC), nbc)
    call Get_x_1st_derivative_P2C_3D(apcc, accc, dm, nbc(:, 2))

    tm%ene_rhs = tm%ene_rhs + accc
!==========================================================================================================
! the RHS of energy equation
! y-pencil : the RHS terms of energy (derivative) operating in the y direction
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
! y-pencil : d (gy * h_cpc) / dy 
!----------------------------------------------------------------------------------------------------------
    ene_rhs_ypencil = ZERO
    call get_ibc_for_calcuation(dm%ibcy_qy(:, IBC_CPC), nbc, dm%ibcy_Th(:, IBC_CPC))
    call Get_y_1st_derivative_P2C_3D( - gy_ypencil * hEnth_cpc_ypencil, accc_ypencil, dm, nbc(:, 1))
    ene_rhs_ypencil = ene_rhs_ypencil + accc_ypencil
!----------------------------------------------------------------------------------------------------------
! y-pencil : d (T) / dy
!----------------------------------------------------------------------------------------------------------
    call Get_y_1st_derivative_C2P_3D(tTemp_ypencil, acpc_ypencil, dm, dm%ibcy_Th(:, IBC_CCC), dm%ftpbcy_var(:, :, :)%t )
!----------------------------------------------------------------------------------------------------------
! y-pencil : k_cpc * d (T) / dy 
!----------------------------------------------------------------------------------------------------------
    acpc_ypencil = acpc_ypencil * kCond_cpc_ypencil
!----------------------------------------------------------------------------------------------------------
! y-pencil : d ( k_cpc * d (T) / dy ) dy
!----------------------------------------------------------------------------------------------------------
    call get_ibc_for_calcuation(dm%ibcy_Th(:, IBC_CPC), nbc)
    call Get_y_1st_derivative_P2C_3D(acpc_ypencil, accc_ypencil, dm, nbc(:, 2))
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
    call get_ibc_for_calcuation(dm%ibcz_qz(:, IBC_CCP), nbc, dm%ibcz_Th(:, IBC_CCP))
    call Get_z_1st_derivative_P2C_3D( - gz_zpencil * hEnth_ccp_zpencil, accc_zpencil, dm, nbc(:, 1))
    ene_rhs_zpencil = ene_rhs_zpencil + accc_zpencil
!----------------------------------------------------------------------------------------------------------
! z-pencil : d (T) / dz
!----------------------------------------------------------------------------------------------------------
    call Get_z_1st_derivative_C2P_3D(tTemp_zpencil, accp_zpencil, dm, dm%ibcz_Th(:, IBC_CCC), dm%ftpbcz_var(:, :, :)%t )
!----------------------------------------------------------------------------------------------------------
! z-pencil : k_ccp * d (T) / dz 
!----------------------------------------------------------------------------------------------------------
    accp_zpencil = accp_zpencil * kCond_ccp_zpencil
!----------------------------------------------------------------------------------------------------------
! z-pencil : d ( k_ccp * d (T) / dz ) / dz
!----------------------------------------------------------------------------------------------------------
    call get_ibc_for_calcuation(dm%ibcz_Th(:, IBC_CCP), nbc)
    call Get_z_1st_derivative_P2C_3D(accp_zpencil, accc_zpencil, dm, nbc(:, 2))
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

end module eq_energy_mod
