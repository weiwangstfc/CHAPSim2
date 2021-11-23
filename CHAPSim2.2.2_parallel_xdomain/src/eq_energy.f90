module eq_energy_mod
  use precision_mod, only : WP
  implicit none

  private :: Calculate_energy_fractional_step
  public :: Solve_energy_eq
contains
!===============================================================================
  subroutine Calculate_energy_fractional_step(rhs0, rhs1, dt, isub)
    use parameters_constant_mod
    implicit none
    real(WP), dimension(:, :, :), intent(inout) :: rhs0, rhs1
    integer(4),                   intent(in   ) :: isub
    real(WP),                     intent(in   ) :: dt

    integer(4) :: n(3)
    real(WP), allocatable :: rhs_dummy(:, :, :)

    n(1:3) = shape(rhs1)
    allocate( rhs_dummy (n(1), n(2), n(3)) )

  ! add explicit terms
    rhs_dummy(:, :, :) = rhs1(:, :, :)
    rhs1(:, :, :) = tGamma(isub) * rhs1(:, :, :) + &
                    tZeta (isub) * rhs0(:, :, :)
    rhs0(:, :, :) = rhs_dummy(:, :, :)

  ! times the time step 
    rhs1(:, :, :) = dt * rhs1(:, :, :)

    deallocate (rhs_dummy)

    return
  end subroutine
!===============================================================================
  subroutine Compute_energy_rhs(fl, tm, dm, isub)
    use operations
    use udf_type_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(in) :: fl
    type(t_thermo), intent(inout) :: tm
    integer,        intent(in) :: isub

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc0
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc0_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc0_zpencil

    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: ene_rhs_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: ene_rhs_zpencil

    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: gy_ypencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: gz_ypencil ! intermediate
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: gz_zpencil ! intermediate
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dppc%xsz(3) ) :: hEnth_xpcc
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dppc%xsz(3) ) :: kCond_xpcc
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: hEnth_ycpc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: hEnth_zccp_zpencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: Ttemp_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: Ttemp_zpencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: kCond_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: kCond_zpencil

!-------------------------------------------------------------------------------
!   preparation
!-------------------------------------------------------------------------------
    call transpose_x_to_y(fl%gy,       gy_ypencil, dm%dcpc)   ! for d(g_y h)/dy
    call transpose_x_to_y(fl%gz,       gz_ypencil, dm%dccp)   ! intermediate, accp_ypencil = gz_ypencil
    call transpose_y_to_z(gz_ypencil,  gz_zpencil, dm%dccp)   ! for d(g_z h)/dz

    do i = 1, 2
      fbc(i) = tm%tpbcx(i)%h
    end do
    call Get_x_midp_C2P_3D(dm%ibcx(5, :), fbc(:), dm, tm%hEnth, hEnth_xpcc ) ! for d(g_y h)/dy
    do i = 1, 2
      fbc(i) = tm%tpbcx(i)%k
    end do
    call Get_x_midp_C2P_3D(dm%ibcx(5, :), fbc(:), dm, tm%kCond, kCond_xpcc ) ! for d(k*(dT/dx))/dx

    call transpose_x_to_y (tm%hEnth,      accc_ypencil, dm%dccc)    !intermediate, accc_ypencil = hEnth_ypencil
    call Get_y_midp_C2P_3D(dm%ibcx(5, :), fbc(:), dm, accc_ypencil, hEnth_ycpc_ypencil)! for d(g_y h)/dy

    call transpose_y_to_z (accc_ypencil,  accc_zpencil, dm%dccc) !intermediate, accc_zpencil = hEnth_zpencil
    call Get_z_midp_C2P_3D(dm%ibcx(5, :), fbc(:), dm, accc_zpencil, hEnth_zccp_zpencil) ! for d(g_z h)/dz

    call transpose_x_to_y (tm%Ttemp,      Ttemp_ypencil, dm%dccc) ! for k d2(T)/dy^2
    call transpose_x_to_y (Ttemp_ypencil, Ttemp_zpencil, dm%dccc)  ! for k d2(T)/dz^2

    call transpose_x_to_y (tm%kCond,      kCond_ypencil, dm%dccc)  ! for k d2(T)/dy^2
    call transpose_x_to_y (kCond_ypencil, kCond_zpencil, dm%dccc) 



!-------------------------------------------------------------------------------
! the RHS of energy equation
! x-pencil : the RHS terms of energy (derivative) operating in the x direction
!-------------------------------------------------------------------------------
    tm%ene_rhs = ZERO
    do i = 1, 2
      fbc(i) = dm%fbcx(1, i) * tm%tpbcx(i)%h
    end do
    ! accc = -d(gx * h)/dx
    call Get_x_1st_derivative_P2C_3D(dm%ibcx(1, :), fbc(:),        dm, - fl%gx * hEnth_xpcc, accc )
    tm%ene_rhs = tm%ene_rhs + accc
    ! apcc = d(T)/dx; accc = d( k * dT/dx)/dx
    do i = 1, 2
      if (dm%ibcx(5, i) == IBC_NEUMANN) then
        ibc(i) = IBC_UNKNOWN
        fbc(i) = ZERO
      else
        ibc(i) = dm%ibcx(5, i)
        fbc(i) = dm%fbcx(5, i)
      end if
    end do
    call Get_x_1st_derivative_C2P_3D(ibc(:), fbc(:), dm, tm%tTemp, apcc )
    do i = 1, 2
      if (dm%ibcx(5, i) == IBC_DIRICHLET) then
        ibc(i) = dm%ibcx(5, i)
        fbc(i) = ZERO
      else
        ibc(i) = dm%ibcx(5, :)
        fbc(i) = dm%fbcx(5, :)
      end if
    end do
    call Get_x_1st_derivative_P2C_3D(dm%ibcx(5, :), dm%fbcx(5, :), dm, apcc * kCond_xpcc, accc )

    tm%ene_rhs = tm%ene_rhs + accc * accc0
!-------------------------------------------------------------------------------
! the RHS of energy equation
! y-pencil : the RHS terms of energy (derivative) operating in the y direction
!-------------------------------------------------------------------------------
    ene_rhs_ypencil = ZERO
    do i = 1, 2
      fbc(i) = dm%fbcy(2, i) * tm%tpbcy(i)%h
    end do
    call Get_y_1st_derivative_P2C_3D(dm%ibcy(2, :), fbc(:), dm, - gy_ypencil * hEnth_ycpc_ypencil, accc_ypencil )
    ene_rhs_ypencil = ene_rhs_ypencil + accc_ypencil

    call Get_y_2nd_derivative_C2C_3D(dm%ibcy(5, :), dm%fbcy(5, :), dm, tTemp_ypencil, accc_ypencil )
    ene_rhs_ypencil = ene_rhs_ypencil + kCond_ypencil * accc_ypencil

    do i = 1, 2
      fbc(i) = tm%tpbcy(i)%k
    end do
    call Get_y_1st_derivative_C2C_3D(dm%ibcy(5, :), fbc(:),        dm, kCond_ypencil, accc_ypencil )
    call Get_y_1st_derivative_C2C_3D(dm%ibcy(5, :), dm%fbcy(5, :), dm, tTemp_ypencil, accc0_ypencil )
    ene_rhs_ypencil = ene_rhs_ypencil + accc_ypencil * accc0_ypencil

    call transpose_y_to_x(ene_rhs_ypencil, accc, dm%dccc)
    tm%ene_rhs = tm%ene_rhs + accc
!-------------------------------------------------------------------------------
! the RHS of energy equation
! z-pencil : the RHS terms of energy (derivative) operating in the z direction
!-------------------------------------------------------------------------------
    ene_rhs_zpencil = ZERO
    do i = 1, 2
      fbc(i) = dm%fbcz(2, i) * tm%tpbcz(i)%h
    end do
    call Get_z_1st_derivative_P2C_3D(dm%ibcz(3, :), fbc(:), dm, - gz_zpencil * hEnth_zccp_zpencil, accc_zpencil )
    ene_rhs_zpencil = ene_rhs_zpencil + accc_zpencil

    call Get_z_2nd_derivative_C2C_3D(dm%ibcz(5, :), dm%fbcz(5, :), dm, tTemp_zpencil, accc_ypencil )
    ene_rhs_zpencil = ene_rhs_zpencil + kCond_zpencil * accc_zpencil

    do i = 1, 2
      fbc(i) = tm%tpbcz(i)%k
    end do
    call Get_z_1st_derivative_C2C_3D(dm%ibcz(5, :), fbc(5, :),     dm, kCond_zpencil, accc_zpencil )
    call Get_z_1st_derivative_C2C_3D(dm%ibcz(5, :), dm%fbcz(5, :), dm, tTemp_zpencil, accc0_zpencil )
    ene_rhs_zpencil = ene_rhs_zpencil + accc_zpencil * accc0_zpencil

    call transpose_z_to_y(ene_rhs_zpencil, ene_rhs_ypencil, dm%dccc)
    call transpose_y_to_x(ene_rhs_ypencil, accc,            dm%dccc)
    tm%ene_rhs = tm%ene_rhs + accc

    call Calculate_energy_fractional_step(tm%ene_rhs0, tm%ene_rhs, dm%dt, isub)

    return
  end subroutine Compute_energy_rhs
!===============================================================================
  subroutine Update_thermal_properties_from_dh(fl, tm, dm)
    use udf_type_mod
    use input_thermo_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm

    integer :: i, j, k
    type(thermoProperty_t) :: tp

    do k = dm%dccc%xst(3), dm%dccc%xen(3)
      do j = dm%dccc%xst(2), dm%dccc%xen(2)
        do i = dm%dccc%xst(1), dm%dccc%xen(1)
          tp%dh = tm%dh
          call tp%Refresh_thermal_properties_from_DH()
          tm%hEnth(i, j, k) = tp%h
          tm%tTemp(i, j, k) = tp%T
          tm%kCond(i, j, k) = tp%k
          fl%dDens(i, j, k) = tp%d
          fl%mVisc(i, j, k) = tp%m
        end do
      end do
    end do

  return
  end subroutine Update_thermal_properties_from_dh
!===============================================================================
  subroutine Solve_energy_eq(fl, tm, dm, isub)
    use udf_type_mod
    
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(in) :: fl
    type(t_thermo), intent(inout) :: tm
    integer,        intent(in) :: isub

!-------------------------------------------------------------------------------
!   calculate rhs of energy equation
!-------------------------------------------------------------------------------
    call Compute_energy_rhs(fl, tm, dm, isub)
!-------------------------------------------------------------------------------
!   update rho * h
!-------------------------------------------------------------------------------
    tm%dh = tm%dh + tm%ene_rhs
!-------------------------------------------------------------------------------
!   update other properties from rho * h
!-------------------------------------------------------------------------------
    call Update_thermal_properties_from_dh(fl, tm, dm)
!-------------------------------------------------------------------------------
!   No Need to apply b.c.
!-------------------------------------------------------------------------------

    !constat heat flux 

  return
  end subroutine

end module eq_energy_mod