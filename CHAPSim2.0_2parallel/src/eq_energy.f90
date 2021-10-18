module eq_energy_mod
  use precision_mod, only : WP
  implicit none
contains

  subroutine Compute_momentum_rhs(f, t, d, isub)
    implicit none

    real(WP), dimension( d%ps_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) :: accc
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) :: accc_ypencil
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) :: accc_zpencil

    real(WP), dimension( d%ps_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) :: accc0
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) :: accc0_ypencil
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) :: accc0_zpencil

    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) :: ene_rhs_ypencil
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) :: ene_rhs_zpencil

    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_yzsz(3) ) :: gy_ypencil
    real(WP), dimension( d%uz_ysz(1), d%uz_ysz(2), d%uz_ysz(3) ) :: gz_ypencil ! intermediate
    real(WP), dimension( d%uz_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) :: gz_zpencil ! intermediate

    real(WP), dimension( d%ux_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) :: hEnth_xpcc
    real(WP), dimension( d%ps_ysz(1), d%uy_ysz(2), d%ps_ysz(3) ) :: hEnth_ycpc_ypencil
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%uz_zsz(3) ) :: hEnth_zccp_zpencil

    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) :: Ttemp_ypencil
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) :: Ttemp_zpencil

    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) :: kCond_ypencil
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) :: kCond_zpencil

    !bc? different from flow b.c.

    call transpose_x_to_y(f%gy,       gy_ypencil, d%dcpc) 
    call transpose_x_to_y(f%gz,       gz_ypencil, d%dccp)  ! intermediate, accp_ypencil = gz_ypencil
    call transpose_y_to_z(gz_ypencil, gz_zpencil, d%dccp) 

    call Get_x_midp_C2P_3dArray ( t%hEnth, d, hEnth_xpcc )
    call transpose_x_to_y(t%hEnth,   accc_ypencil, d%dccc)    !intermediate, accc_ypencil = hEnth_ypencil
    call Get_y_midp_C2P_3dArray ( accc_ypencil, d, hEnth_ycpc_ypencil)
    call transpose_y_to_z(accc_ypencil, accc_zpencil, d%dccc)
    call Get_z_midp_C2P_3dArray ( accc_zpencil, d, hEnth_zccp_zpencil)

    call transpose_x_to_y(t%Ttemp,       Ttemp_ypencil, d%dccc) 
    call transpose_x_to_y(Ttemp_ypencil, Ttemp_zpencil, d%dccc) 

    call transpose_x_to_y(t%kCond,       kCond_ypencil, d%dccc) 
    call transpose_x_to_y(kCond_ypencil, kCond_zpencil, d%dccc) 

!===============================================================================
! the RHS of energy equation
! x-pencil : the RHS terms of energy (derivative) operating in the x direction
!===============================================================================
    t%ene_rhs = ZERO

    call Get_x_1st_derivative_P2C_3dArray( - f%gx * hEnth_xpcc, d, accc )
    t%ene_rhs = t%ene_rhs + accc

    call Get_x_2nd_derivative_C2C_3dArray( t%tTemp,  d, accc )
    t%ene_rhs = t%ene_rhs + t%kCond * accc

    call Get_x_1st_derivative_C2C_3dArray( t%kCond, d, accc )
    call Get_x_1st_derivative_C2C_3dArray( t%tTemp, d, accc0 )
    t%ene_rhs = t%ene_rhs + accc * accc0
!===============================================================================
! the RHS of energy equation
! y-pencil : the RHS terms of energy (derivative) operating in the y direction
!===============================================================================
    ene_rhs_ypencil = ZERO
    call Get_y_1st_derivative_P2C_3dArray( - gy_ypencil * hEnth_ycpc_ypencil, d, accc_ypencil )
    ene_rhs_ypencil = ene_rhs_ypencil + accc_ypencil

    call Get_y_2nd_derivative_C2C_3dArray( tTemp_ypencil,  d, accc_ypencil )
    ene_rhs_ypencil = ene_rhs_ypencil + kCond_ypencil * accc_ypencil

    call Get_y_1st_derivative_C2C_3dArray( kCond_ypencil, d, accc_ypencil )
    call Get_y_1st_derivative_C2C_3dArray( tTemp_ypencil, d, accc0_ypencil )
    ene_rhs_ypencil = ene_rhs_ypencil + accc_ypencil * accc0_ypencil

    call transpose_y_to_x(ene_rhs_ypencil, accc, d%dccc)
    t%ene_rhs = t%ene_rhs + accc
!===============================================================================
! the RHS of energy equation
! z-pencil : the RHS terms of energy (derivative) operating in the z direction
!===============================================================================
    ene_rhs_zpencil = ZERO
    call Get_z_1st_derivative_P2C_3dArray( - gz_zpencil * hEnth_zccp_zpencil, d, accc_zpencil )
    ene_rhs_zpencil = ene_rhs_zpencil + accc_zpencil

    call Get_z_2nd_derivative_C2C_3dArray( tTemp_zpencil,  d, accc_ypencil )
    ene_rhs_zpencil = ene_rhs_zpencil + kCond_zpencil * accc_zpencil

    call Get_z_1st_derivative_C2C_3dArray( kCond_zpencil, d, accc_zpencil )
    call Get_z_1st_derivative_C2C_3dArray( tTemp_zpencil, d, accc0_zpencil )
    ene_rhs_zpencil = ene_rhs_zpencil + accc_zpencil * accc0_zpencil

    call transpose_z_to_y(ene_rhs_zpencil, ene_rhs_ypencil, d%dccc)
    call transpose_y_to_x(ene_rhs_ypencil, accc,            d%dccc)
    t%ene_rhs = t%ene_rhs + accc

  end subroutine 

end module eq_energy_mod