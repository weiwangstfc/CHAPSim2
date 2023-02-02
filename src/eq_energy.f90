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
    do i = 1, 2
      fbcx(i) = tm%ftpbcx(i)%h
      fbcy(i) = tm%ftpbcy(i)%h
      fbcz(i) = tm%ftpbcz(i)%h
    end do
    call Get_x_midp_C2P_3D(tm%hEnth,     hEnth_pcc,         dm, dm%ibcx(:, 5), fbcx(:)) ! for d(g_x h_pcc))/dy
    call transpose_x_to_y (tm%hEnth,     accc_ypencil, dm%dccc)                     !intermediate, accc_ypencil = hEnth_ypencil
    call Get_y_midp_C2P_3D(accc_ypencil, hEnth_cpc_ypencil, dm, dm%ibcy(:, 5), fbcy(:))! for d(g_y h_cpc)/dy
    call transpose_y_to_z (accc_ypencil, accc_zpencil, dm%dccc) !intermediate, accc_zpencil = hEnth_zpencil
    call Get_z_midp_C2P_3D(accc_zpencil, hEnth_ccp_zpencil, dm, dm%ibcz(:, 5), fbcz(:)) ! for d(g_z h_ccp)/dz
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
    call Get_x_midp_C2P_3D(tm%kCond,      kCond_pcc,     dm, dm%ibcx(:, 5), fbcx(:) ) ! for d(k_pcc * (dT/dx) )/dx
    call transpose_x_to_y (tm%kCond,      accc_ypencil, dm%dccc)  ! for k d2(T)/dy^2
    call Get_y_midp_C2P_3D(accc_ypencil,  kCond_cpc_ypencil, dm, dm%ibcy(:, 5), fbcy(:))
    call transpose_y_to_z (accc_ypencil,  kCond_zpencil, dm%dccc) 
    call Get_z_midp_C2P_3D(kCond_zpencil, kCond_ccp_zpencil, dm, dm%ibcz(:, 5), fbcz(:))
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
    call Get_x_1st_derivative_P2C_3D( - fl%gx * hEnth_pcc, accc, dm, dm%ibcx(:, 1) ) ! accc = -d(gx * h)/dx
    tm%ene_rhs = tm%ene_rhs + accc
!----------------------------------------------------------------------------------------------------------
! x-pencil : d (T) / dx 
!----------------------------------------------------------------------------------------------------------
    do i = 1, 2
      if (dm%ibcx(i, 5) == IBC_NEUMANN) then
        ibcx(i) = IBC_INTERIOR
        fbcx(i) = ZERO
      else
        ibcx(i) = dm%ibcx(i, 5)
        fbcx(i) = dm%fbcx(i, 5)
      end if
    end do
    call Get_x_1st_derivative_C2P_3D(tm%tTemp, apcc, dm, ibcx(:), fbcx(:) )
!----------------------------------------------------------------------------------------------------------
! x-pencil : k_pcc * d (T) / dx 
!----------------------------------------------------------------------------------------------------------
    apcc = apcc * kCond_pcc
    if (dm%ibcx(1, 5) == IBC_NEUMANN) then
      apcc(1, :, :) = dm%fbcx(1, 5)
    end if
    if (dm%ibcx(2, 5) == IBC_NEUMANN) then
      apcc(dm%dpcc%xen(1), :, :) = dm%fbcx(2, 5)
    end if
!----------------------------------------------------------------------------------------------------------
! x-pencil : d ( k_pcc * d (T) / dx ) dx
!----------------------------------------------------------------------------------------------------------
    call Get_x_1st_derivative_P2C_3D(apcc, accc, dm, dm%ibcx(:, 5) )

    tm%ene_rhs = tm%ene_rhs + accc
!==========================================================================================================
! the RHS of energy equation
! y-pencil : the RHS terms of energy (derivative) operating in the y direction
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
! y-pencil : d (gy * h_cpc) / dy 
!----------------------------------------------------------------------------------------------------------
    ene_rhs_ypencil = ZERO
    call Get_y_1st_derivative_P2C_3D( - gy_ypencil * hEnth_cpc_ypencil, accc_ypencil, dm, dm%ibcy(:, 2) )
    ene_rhs_ypencil = ene_rhs_ypencil + accc_ypencil
!----------------------------------------------------------------------------------------------------------
! y-pencil : d (T) / dy
!----------------------------------------------------------------------------------------------------------
    do i = 1, 2
      if (dm%ibcy(i, 5) == IBC_NEUMANN) then
        ibcy(i) = IBC_INTERIOR
        fbcy(i) = ZERO
      else
        ibcy(i) = dm%ibcy(i, 5)
        fbcy(i) = dm%fbcy(i, 5)
      end if
    end do
    call Get_y_1st_derivative_C2P_3D(tTemp_ypencil, acpc_ypencil, dm, ibcy(:), fbcy(:) )
!----------------------------------------------------------------------------------------------------------
! y-pencil : k_cpc * d (T) / dy 
!----------------------------------------------------------------------------------------------------------
    acpc_ypencil = acpc_ypencil * kCond_cpc_ypencil
    if (dm%ibcy(1, 5) == IBC_NEUMANN) then
      acpc_ypencil(:, 1, :) = dm%fbcy(1, 5)
    end if
    if (dm%ibcy(2, 5) == IBC_NEUMANN) then
      acpc_ypencil(:, dm%dcpc%yen(2), :) = dm%fbcy(2, 5)
    end if
!----------------------------------------------------------------------------------------------------------
! y-pencil : d ( k_cpc * d (T) / dy ) dy
!----------------------------------------------------------------------------------------------------------
    call Get_y_1st_derivative_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%ibcy(:, 5) )
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
    call Get_z_1st_derivative_P2C_3D( - gz_zpencil * hEnth_ccp_zpencil, accc_zpencil, dm, dm%ibcz(:, 3) )
    ene_rhs_zpencil = ene_rhs_zpencil + accc_zpencil
!----------------------------------------------------------------------------------------------------------
! z-pencil : d (T) / dz
!----------------------------------------------------------------------------------------------------------
    do i = 1, 2
      if (dm%ibcz(i, 5) == IBC_NEUMANN) then
        ibcz(i) = IBC_INTERIOR
        fbcz(i) = ZERO
      else
        ibcz(i) = dm%ibcz(i, 5)
        fbcz(i) = dm%fbcz(i, 5)
      end if
    end do
    call Get_z_1st_derivative_C2P_3D(tTemp_zpencil, accp_zpencil, dm, ibcz(:), fbcz(:) )
!----------------------------------------------------------------------------------------------------------
! z-pencil : k_ccp * d (T) / dz 
!----------------------------------------------------------------------------------------------------------
    accp_zpencil = accp_zpencil * kCond_ccp_zpencil
    if (dm%ibcz(1, 5) == IBC_NEUMANN) then
      accp_zpencil(:, 1, :) = dm%fbcz(1, 5)
    end if
    if (dm%ibcz(2, 5) == IBC_NEUMANN) then
      accp_zpencil(:, :, dm%dccp%zen(3)) = dm%fbcz(2, 5)
    end if
!----------------------------------------------------------------------------------------------------------
! z-pencil : d ( k_ccp * d (T) / dz ) / dz
!----------------------------------------------------------------------------------------------------------
    call Get_z_1st_derivative_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%ibcz(:, 5) )
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
    
    integer :: i, j, k
    type(DECOMP_INFO) :: dtmp

    integer  :: ibc(2)
    real(WP) :: fbc(2)

!----------------------------------------------------------------------------------------------------------
! x-pencil : u1 -> g1
!----------------------------------------------------------------------------------------------------------
    dtmp = dm%dpcc
    ibc(:) = dm%ibcx(:, 5)
    fbc(:) = dm%fbc_dend(:, 1)
    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
        fix(:) = fl%dDens(:, j, k)
        call Get_x_midp_C2P_1D (fix, fox, dm, ibc, fbc)
        fl%gx(:, j, k) = fox(:) * fl%qx(:, j, k)
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! x-pencil --> y-pencil
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qy,    uy_ypencil, dm%dcpc)
    call transpose_x_to_y(fl%dDens,  d_ypencil, dm%dccc)
    call transpose_x_to_y(fl%qz,    uz_ypencil, dm%dccp)
!----------------------------------------------------------------------------------------------------------
! y-pencil : u2 -> g2
!----------------------------------------------------------------------------------------------------------
    dtmp = dm%dcpc
    ibc(:) = dm%ibcy(:, 5)
    fbc(:) = dm%fbc_dend(:, 2)
    do k = 1, dtmp%ysz(3)
      do i = 1, dtmp%ysz(1)
        fiy(:) = d_ypencil(i, :, k)
        call Get_y_midp_C2P_1D (fiy, foy, dm, ibc, fbc)
        duy_ypencil(i, :, k) = foy(:) * uy_ypencil(i, :, k)
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! y-pencil --> z-pencil
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_z( d_ypencil,  d_zpencil, dm%dccc)
    call transpose_y_to_z(uz_ypencil, uz_zpencil, dm%dccp)
!----------------------------------------------------------------------------------------------------------
! Z-pencil : u3 -> g3
!----------------------------------------------------------------------------------------------------------
    dtmp = dm%dcpc
    ibc(:) = dm%ibcz(:, 5)
    fbc(:) = dm%fbc_dend(:, 3)
    do j = 1, dtmp%zsz(2)
      do i = 1, dtmp%zsz(1)
        fiz(:) = d_zpencil(i, j, :)
        call Get_z_midp_C2P_1D (fiz, foz, dm, ibc, fbc)
        duz_zpencil(i, j, :) = foz(:) * uz_zpencil(i, j, :)
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! z-pencil --> y-pencil
!----------------------------------------------------------------------------------------------------------
    call transpose_z_to_y(duz_zpencil, duz_ypencil, dm%dccp)
!----------------------------------------------------------------------------------------------------------
! y-pencil --> x-pencil
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_x(duz_ypencil, fl%gz, dm%dccp)
    call transpose_y_to_x(duy_ypencil, fl%gy, dm%dcpc)

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
    real(WP), dimension( dm%nc(1) ) :: fix
    real(WP), dimension( dm%np(1) ) :: fox
    real(WP), dimension( dm%nc(2) ) :: fiy
    real(WP), dimension( dm%np(2) ) :: foy
    real(WP), dimension( dm%nc(3) ) :: fiz
    real(WP), dimension( dm%np(3) ) :: foz
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
                         dm%dccc%ysz(3)) :: d_ypencil
    
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: qz_zpencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: gz_zpencil
    real(WP), dimension( dm%dccc%zsz(1), &
                         dm%dccc%zsz(2), &
                         dm%dccc%zsz(3)) :: d_zpencil
    
    integer :: i, j, k
    type(DECOMP_INFO) :: dtmp

    integer  :: ibc(2)
    real(WP) :: fbc(2)

!----------------------------------------------------------------------------------------------------------
! x-pencil : g1 -> u1
!----------------------------------------------------------------------------------------------------------
    dtmp = dm%dpcc
    ibc(:) = dm%ibcx(:, 5)
    fbc(:) = dm%fbc_dend(:, 1)
    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
        fix(:) = fl%dDens(:, j, k)
        call Get_x_midp_C2P_1D (fix, fox, dm, ibc, fbc)
        fl%qx(:, j, k) = fl%gx(:, j, k) / fox(:)
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! x-pencil --> y-pencil
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%gy,    gy_ypencil, dm%dcpc)
    call transpose_x_to_y(fl%dDens,  d_ypencil, dm%dccc)
    call transpose_x_to_y(fl%gz,    gz_ypencil, dm%dccp)
!----------------------------------------------------------------------------------------------------------
! y-pencil : g2 -> u2
!----------------------------------------------------------------------------------------------------------
    dtmp = dm%dcpc
    ibc(:) = dm%ibcy(:, 5)
    fbc(:) = dm%fbc_dend(:, 2)
    do k = 1, dtmp%ysz(3)
      do i = 1, dtmp%ysz(1)
        fiy(:) = d_ypencil(i, :, k)
        call Get_y_midp_C2P_1D (fiy, foy, dm, ibc, fbc)
        qy_ypencil(i, :, k) = gy_ypencil(i, :, k) / foy(:) 
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! y-pencil --> z-pencil
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_z( d_ypencil,  d_zpencil, dm%dccc)
    call transpose_y_to_z(gz_ypencil, gz_zpencil, dm%dccp)
!----------------------------------------------------------------------------------------------------------
! Z-pencil : u3 -> g3
!----------------------------------------------------------------------------------------------------------
    dtmp = dm%dcpc
    ibc(:) = dm%ibcz(:, 5)
    fbc(:) = dm%fbc_dend(:, 3)
    do j = 1, dtmp%zsz(2)
      do i = 1, dtmp%zsz(1)
        fiz(:) = d_zpencil(i, j, :)
        call Get_z_midp_C2P_1D (fiz, foz, dm, ibc, fbc)
        qz_zpencil(i, j, :) = gz_zpencil(i, j, :) / foz(:)
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! z-pencil --> y-pencil
!----------------------------------------------------------------------------------------------------------
    call transpose_z_to_y(qz_zpencil, qz_ypencil, dm%dccp)
!----------------------------------------------------------------------------------------------------------
! y-pencil --> x-pencil
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_x(qz_ypencil, fl%qz, dm%dccp)
    call transpose_y_to_x(qy_ypencil, fl%qy, dm%dcpc)

    return
  end subroutine Calculate_velocity_from_massflux

end module eq_energy_mod
