module eq_momentum_mod
  use precision_mod, only : WP
  implicit none
  private :: Calculate_momentum_driven_source
  private :: Calculate_momentum_fractional_step
  private :: Compute_momentum_rhs
  private :: Correct_massflux
  public  :: Solve_momentum_eq

contains
!===============================================================================
!===============================================================================
!> \brief To calcuate the driven force for streamwise peridic flow.
!>
!> This subroutine is called everytime when calcuting the rhs of momentum eqs.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                         !
!______________________________________________________________________________!
!> \param[inout]  rhs          the rhs of the momentum equation
!> \param[in]     d            domain name
!> \param[in]     isub         the RK iteration to get correct Coefficient 
!_______________________________________________________________________________
  subroutine Calculate_momentum_driven_source(rhs, str, d, isub)
    use input_general_mod,       only : idriven, drvf,  IDRVF_NO, IDRVF_MASSFLUX, &
                                        IDRVF_SKINFRIC, IDRVF_PRESLOSS, &
                                        tAlpha, dt
    use operations,              only : Get_volumetric_average_3d
    use udf_type_mod,            only : t_domain
    use parameters_constant_mod, only : HALF, ZERO
    implicit none

    type(t_domain), intent(in) :: d
    real(WP),    intent(inout) :: rhs(:, :, :)
    integer(4),     intent(in) :: isub
    character(2),   intent(in) :: str
    
    real(WP) :: rhs_bulk
    logical :: is_stored_nyp

    rhs_bulk = ZERO

    if(idriven == IDRVF_MASSFLUX) then

      is_stored_nyp = .false.
      call Get_volumetric_average_3d(rhs, str, d, rhs_bulk)
      
    else if (idriven == IDRVF_SKINFRIC) then

      rhs_bulk = - HALF * drvf * tAlpha(isub) * dt

    else if (idriven == IDRVF_PRESLOSS ) then

    ! to check this part
      rhs_bulk = - HALF * drvf * tAlpha(isub) * dt

    else 
      return
    end if

    rhs(:, :, :) = rhs(:, :, :) - rhs_bulk

    return
  end subroutine 
!===============================================================================
!===============================================================================
!> \brief To calcuate the convection and diffusion terms in rhs of momentum eq.
!>
!> This subroutine is called everytime when calcuting the rhs of momentum eqs.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                         !
!______________________________________________________________________________!
!> \param[inout]  rhs0          the last iteration rhs
!> \param[inout]  rhs1          the current iteration rhs
!> \param[in]     rhs1_semi     the semi-implicit term
!> \param[in]     isub          the RK iteration to get correct Coefficient 
!_______________________________________________________________________________
  subroutine Calculate_momentum_fractional_step(rhs0, rhs1, rhs1_semi, isub)
    use input_general_mod, only : tGamma, tZeta, tAlpha, dt, &
                                  IVIS_SEMIMPLT, iviscous
    implicit none
    real(WP), dimension(:, :, :), intent(in   ) :: rhs1_semi
    real(WP), dimension(:, :, :), intent(inout) :: rhs0, rhs1
    integer(4),                   intent(in   ) :: isub

    integer(4) :: n(3)
    real(WP), allocatable :: rhs_dummy(:, :, :)

    n(1:3) = shape(rhs1)
    allocate( rhs_dummy (n(1), n(2), n(3)) )

  ! add explicit terms
    rhs_dummy(:, :, :) = rhs1(:, :, :)
    rhs1(:, :, :) = tGamma(isub) * rhs1(:, :, :) + &
                    tZeta (isub) * rhs0(:, :, :)
    rhs0(:, :, :) = rhs_dummy(:, :, :)

  ! add implicit
    rhs1(:, :, :) = rhs1(:, :, :) + &
                    tAlpha(isub) * rhs1_semi(:, :, :)
  
  ! times the time step 
    rhs1(:, :, :) = dt * rhs1(:, :, :)

    deallocate (rhs_dummy)

    return
  end subroutine
!===============================================================================
!===============================================================================
!> \brief To calcuate all rhs of momentum eq.
!>
!> This subroutine is called everytime when calcuting the rhs of momentum eqs.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                         !
!______________________________________________________________________________!
!> \param[inout]  f             flow field
!> \param[inout]  d             domain    
!> \param[in]     isub          the RK iteration to get correct Coefficient 
!_______________________________________________________________________________
  subroutine Compute_momentum_rhs(f, t, d, isub)
    use input_general_mod,       only : ithermo, igravity, idriven, &
                                        sigma1p, &   
                                        IDRVF_NO, IVIS_EXPLICIT, IVIS_SEMIMPLT
    use udf_type_mod,            only : t_flow, t_domain
    use parameters_constant_mod, only : TWO, ZERO, ONE_THIRD, TWO_THIRD
    use operations
    implicit none

    type(t_flow),   intent(inout) :: f
    type(t_domain), intent(in   ) :: d
    integer(4),     intent(in   ) :: isub
!-------------------------------------------------------------------------------
! common vars
!-------------------------------------------------------------------------------
    real(WP), dimension( d%ux_ysz(1), d%ux_ysz(2), d%ux_ysz(3) ) ::      qx_ypencil
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) ::      qy_ypencil
    real(WP), dimension( d%uz_ysz(1), d%uz_ysz(2), d%uz_ysz(3) ) ::      qz_ypencil 

    real(WP), dimension( d%ux_zsz(1), d%ux_zsz(2), d%ux_zsz(3) ) ::      qx_zpencil
    real(WP), dimension( d%uy_zsz(1), d%uy_zsz(2), d%uy_zsz(3) ) ::      qy_zpencil
    real(WP), dimension( d%uz_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) ::      qz_zpencil

    real(WP), dimension( d%ux_ysz(1), d%uy_ysz(2), d%ux_ysz(3) ) :: qx_yppc_ypencil ! <ux>^y at (xp, yp, zc)
    real(WP), dimension( d%ux_zsz(1), d%ux_zsz(2), d%uz_zsz(3) ) :: qx_zpcp_zpencil ! <ux>^z at (xp, yc, zp)
    real(WP), dimension( d%uy_zsz(1), d%uy_zsz(2), d%uz_zsz(3) ) :: qy_zcpp_zpencil ! <uy>^z at (xc, yp, zp)
    real(WP), dimension( d%ux_xsz(1), d%uz_xsz(2), d%uz_xsz(3) ) :: qz_xpcp         ! <uz>^x at (xp, yc, zp) 
    real(WP), dimension( d%uz_ysz(1), d%uy_ysz(2), d%uz_ysz(3) ) :: qz_ycpp_ypencil ! <uz>^y at (xc, yp, zp)
    
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) ::    pres_ypencil ! p
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) ::    pres_zpencil ! p

    real(WP), dimension( d%ux_ysz(1), d%ux_ysz(2), d%ux_ysz(3) ) ::  mx_rhs_ypencil ! 
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) ::  my_rhs_ypencil ! 
    real(WP), dimension( d%uz_ysz(1), d%uz_ysz(2), d%uz_ysz(3) ) ::  mz_rhs_ypencil ! 
 
    real(WP), dimension( d%ux_zsz(1), d%ux_zsz(2), d%ux_zsz(3) ) ::  mx_rhs_zpencil ! 
    real(WP), dimension( d%uy_zsz(1), d%uy_zsz(2), d%uy_zsz(3) ) ::  my_rhs_zpencil ! 
    real(WP), dimension( d%uz_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) ::  mz_rhs_zpencil ! 
 
    real(WP), dimension( d%uz_xsz(1), d%uz_xsz(2), d%uz_xsz(3) ) :: mx_rhs_implicit ! 
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) :: my_rhs_implicit_ypencil ! 
    real(WP), dimension( d%uz_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) :: mz_rhs_implicit_zpencil ! 

    real(WP), dimension( d%ux_xsz(1), d%ux_xsz(2), d%ux_xsz(3) ) :: apcc ! 
    real(WP), dimension( d%ux_ysz(1), d%ux_ysz(2), d%ux_ysz(3) ) :: apcc_ypencil ! 
    real(WP), dimension( d%ux_zsz(1), d%ux_zsz(2), d%ux_zsz(3) ) :: apcc_zpencil ! 

    real(WP), dimension( d%uy_xsz(1), d%uy_xsz(2), d%uy_xsz(3) ) :: acpc ! 
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) :: acpc_ypencil !
    real(WP), dimension( d%uy_zsz(1), d%uy_zsz(2), d%uy_zsz(3) ) :: acpc_zpencil !  

    real(WP), dimension( d%uz_xsz(1), d%uz_xsz(2), d%uz_xsz(3) ) :: accp ! 
    real(WP), dimension( d%uz_ysz(1), d%uz_ysz(2), d%uz_ysz(3) ) :: accp_ypencil !
    real(WP), dimension( d%uz_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) :: accp_zpencil ! 

    real(WP), dimension( d%ux_ysz(1), d%ux_ysz(2), d%uz_ysz(3) ) :: apcp_ypencil ! 
    
!-------------------------------------------------------------------------------
! thermal == 0 only
!-------------------------------------------------------------------------------
    real(WP), dimension( d%ux_xsz(1), d%ux_xsz(2), d%uz_xsz(3) ) :: qx_zpcp ! <ux>^z at (xp, yc, zp)
    real(WP), dimension( d%ux_xsz(1), d%uy_xsz(2), d%ux_xsz(3) ) :: qx_yppc ! <ux>^y at (xp, yp, zc)
    real(WP), dimension( d%ux_xsz(1), d%uy_xsz(2), d%uy_xsz(3) ) :: qy_xppc ! <uy>^x at (xp, yp, zc)
    real(WP), dimension( d%ux_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) :: qy_xppc_ypencil ! <uy>^x at (xp, yp, zc)
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uz_ysz(3) ) :: qy_zcpp_ypencil ! <uy>^z at (xc, yp, zp)
    real(WP), dimension( d%ux_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) :: qz_xpcp_zpencil ! <uz>^x at (xp, yc, zp)
    real(WP), dimension( d%uz_zsz(1), d%uy_zsz(2), d%uz_zsz(3) ) :: qz_ycpp_zpencil ! <uz>^y at (xc, yp, zp)
!-------------------------------------------------------------------------------
! thermal == 1 only
!-------------------------------------------------------------------------------
    real(WP), dimension( d%ps_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) :: accc
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) :: accc_ypencil
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) :: accc_zpencil
    real(WP), dimension( d%ux_xsz(1), d%ux_xsz(2), d%uz_xsz(3) ) :: apcp
    real(WP), dimension( d%ux_zsz(1), d%ux_zsz(2), d%uz_zsz(3) ) :: apcp_zpencil
    real(WP), dimension( d%ux_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) :: appc_ypencil
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uz_ysz(3) ) :: acpp_ypencil
    real(WP), dimension( d%uy_zsz(1), d%uy_zsz(2), d%uz_zsz(3) ) :: acpp_zpencil
    
    real(WP), dimension( d%ps_ysz(1), d%ux_ysz(2), d%ux_ysz(3) ) :: qx_xccc_ypencil   ! <ux>^x at (xc, yc, zc)
    real(WP), dimension( d%ps_zsz(1), d%ux_zsz(2), d%ux_zsz(3) ) :: qx_xccc_zpencil   ! <ux>^x at (xc, yc, zc)
    real(WP), dimension( d%uy_xsz(1), d%ps_xsz(2), d%uy_xsz(3) ) :: qy_yccc           ! <uy>^y at (xc, yc, zc)
    real(WP), dimension( d%uy_zsz(1), d%ps_zsz(2), d%uy_zsz(3) ) :: qy_yccc_zpencil   ! <uy>^y at (xc, yc, zc)
    real(WP), dimension( d%uz_xsz(1), d%ps_xsz(2), d%uz_xsz(3) ) :: qz_zccc           ! <uz>^z at (xc, yc, zc)
    real(WP), dimension( d%uz_ysz(1), d%ps_ysz(2), d%uz_ysz(3) ) :: qz_zccc_ypencil   ! <uz>^z at (xc, yc, zc), intermediate

    real(WP), dimension( d%ux_xsz(1), d%uy_xsz(2), d%ux_xsz(3) ) :: gx_yppc           ! <gx>^y at (xp, yp, zc)
    real(WP), dimension( d%ux_xsz(1), d%ux_xsz(2), d%uz_xsz(3) ) :: gx_zpcp           ! <gx>^z at (xp, yc, zp)
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) :: gy_ypencil
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uz_ysz(3) ) :: gy_zcpp_ypencil   ! <gy>^z at (xc, yp, zp)
    real(WP), dimension( d%ux_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) :: gz_xpcp_zpencil   ! <gz>^x at (xp, yc, zp)
    real(WP), dimension( d%uz_zsz(1), d%uy_zsz(2), d%uz_zsz(3) ) :: gz_ycpp_zpencil   ! <gz>^y at (xc, yp, zp)
    
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) ::    dDens_ypencil  ! d 
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) ::    dDens_zpencil  ! d 

    real(WP), dimension( d%ux_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) ::    m_xpcc         ! <mu>^x       at (xp, yc, zc)  
    real(WP), dimension( d%ux_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) ::    m_xpcc_ypencil ! <mu>^x       at (xp, yc, zc)
    real(WP), dimension( d%ux_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) ::    m_xpcc_zpencil ! <mu>^x       at (xp, yc, zc)
    real(WP), dimension( d%ps_xsz(1), d%uy_xsz(2), d%ps_xsz(3) ) ::    m_ycpc         ! <mu>^y       at (xc, yp, zc)
    real(WP), dimension( d%ps_ysz(1), d%uy_ysz(2), d%ps_ysz(3) ) ::    m_ycpc_ypencil ! <mu>^y       at (xc, yp, zc)
    real(WP), dimension( d%ps_zsz(1), d%uy_zsz(2), d%ps_zsz(3) ) ::    m_ycpc_zpencil ! <mu>^y       at (xc, yp, zc)
    real(WP), dimension( d%ps_xsz(1), d%ps_xsz(2), d%uz_xsz(3) ) ::    m_zccp         ! <mu>^z       at (xc, yc, zp)
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%uz_ysz(3) ) ::    m_zccp_ypencil ! <mu>^z       at (xc, yc, zp)
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%uz_zsz(3) ) ::    m_zccp_zpencil ! <mu>^z       at (xc, yc, zp)
    
    real(WP), dimension( d%ux_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) :: dmdx_xpcc ! d( mu   )/dx at (xp, yc, zc)
    real(WP), dimension( d%ps_xsz(1), d%uy_xsz(2), d%ps_xsz(3) ) :: dmdx_ycpc ! d(<mu>^y)/dx at (xc, yp, zc)
    real(WP), dimension( d%ps_xsz(1), d%ps_xsz(2), d%uz_xsz(3) ) :: dmdx_zccp ! d(<mu>^z)/dx at (xc, yc, zp)
    real(WP), dimension( d%ux_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) :: dmdy_xpcc ! d(<mu>^x)/dy at (xp, yc, zc)
    real(WP), dimension( d%ux_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) :: dmdz_xpcc ! d(<mu>^x)/dz at (xp, yc, zc)

    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%uz_zsz(3) ) :: dmdx_zccp_zpencil ! d(<mu>^z)/dx at (xc, yc, zp)
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%uz_ysz(3) ) :: dmdy_zccp_ypencil ! d(<mu>^z)/dy at (xc, yc, zp)
    real(WP), dimension( d%ux_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) :: dmdy_xpcc_ypencil ! d(<mu>^x)/dy at (xp, yc, zc)
    real(WP), dimension( d%ux_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) :: dmdz_xpcc_zpencil ! d(<mu>^x)/dz at (xp, yc, zc)
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%uz_zsz(3) ) :: dmdz_zccp_zpencil ! d( mu   )/dz at (xc, yc, zp)

    real(WP), dimension( d%ps_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) ::       div
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) ::       div_ypencil
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) ::       div_zpencil
!-------------------------------------------------------------------------------
! others
!-------------------------------------------------------------------------------
    real(WP)              :: one_third_rre, two_third_rre, two_rre
!-------------------------------------------------------------------------------
!   Ensure it is in x-pencil
!-------------------------------------------------------------------------------
    if(d%ux_xsz(1) /= d%np(1)) call Print_error_msg("Error, not X-pencil")
!===============================================================================
! variable preparation
! In the comments: 
!     I = Intermediate
!     O = no thermal fields
!     W = with thermal fields
!     WO = both, commom
!===============================================================================
!-------------------------------------------------------------------------------
!    p --> p_ypencil --> p_zpencil
!-------------------------------------------------------------------------------
    call transpose_x_to_y(f%pres, pres_ypencil, d%dccc)           ! y-pencil : y-mom, w+o   thermal
    call transpose_y_to_z(pres_ypencil, pres_zpencil, d%dccc)     ! z-pencil : z-mom, w+o   thermal
!-------------------------------------------------------------------------------
!    qx --> qx_xccc (I) --> qx_xccc_ypencil(W) --> qx_xccc_zpencil(W)
!     | --> qx_ypencil(WO) --> qx_yppc_ypencil(WO) --> qx_yppc(O)
!                       |  --> qx_zpencil(WO) --> qx_zpcp_zpencil(WO) --> qx_zpcp_ypencil(I) --> qx_zpcp(O)
!------------------------------------------------------------------------------- 
    call transpose_x_to_y(f%qx, qx_ypencil, d%dpcc)               ! y-pencil : x-mom, w+o thermal
    call Get_y_midp_C2P_3dArray ( qx_ypencil, d, qx_yppc_ypencil )! y-pencil : x-mom, w+o thermal
    if(ithermo == 0) &
    call transpose_y_to_x(qx_yppc_ypencil, qx_yppc, d%dppc)       ! x-pencil : y-mom, o   thermal
  
    call transpose_y_to_z(qx_ypencil, qx_zpencil, d%dpcc)         ! z-pencil : x-mom, w+o thermal
    call Get_z_midp_C2P_3dArray ( qx_zpencil, d, qx_zpcp_zpencil )! z-pencil : x-mom, w+o thermal
    if(ithermo == 0) then
    call transpose_z_to_y(qx_zpcp_zpencil, apcp_ypencil, d%dpcp)  ! intermediate, apcp_ypencil = qx_zpcp_ypencil
    call transpose_y_to_x(   apcp_ypencil,      qx_zpcp, d%dpcp)  ! x-pencil : z-mom,  o  thermal
    end if
    if(ithermo == 1) then
    call Get_x_midp_P2C_3dArray ( f%qx,    d, accc )              ! intermediate, accc = qx_xccc
    call transpose_x_to_y(accc, qx_xccc_ypencil, d%dccc)          ! y-pencil : y-mom, w   thermal
    call transpose_y_to_z(qx_xccc_ypencil, qx_xccc_zpencil, d%dccc) ! z-pencil : z-mom, w   thermal
    end if
!-------------------------------------------------------------------------------
!    qy--> qy_xppc(WO) --> qy_xppc_ypencil(O) 
!     |--> qy_ypencil(WO) --> qy_zpencil(WO) --> qy_zcpp_zpencil(WO) --> qy_zcpp_ypencil(O)
!                       | --> qy_yccc_ypencil(I) --> qy_yccc(W)
!                                              | --> qy_yccc_zpencil(W)              
!------------------------------------------------------------------------------- 
    call Get_x_midp_C2P_3dArray ( f%qy, d, qy_xppc )              ! xpencil : y-mom, w+o thermal
    if(ithermo == 0) &
    call transpose_x_to_y(qy_xppc, qy_xppc_ypencil, d%dppc)       ! ypencil : x-mom, o thermal

    call transpose_x_to_y(f%qy, qy_ypencil, d%dcpc)               ! y-pencil : y-mom, w+o thermal
    call transpose_y_to_z(qy_ypencil, qy_zpencil, d%dcpc)         ! z-pencil : y-mom, w+o thermal
    call Get_z_midp_C2P_3dArray ( qy_zpencil, d, qy_zcpp_zpencil )! z-pencil : y-mom, w+o thermal
    if(ithermo == 0) &
    call transpose_z_to_y(qy_zcpp_zpencil, qy_zcpp_ypencil, d%dcpp) ! y-pencil : z-mom, o thermal

    if(ithermo == 1) then
    call Get_y_midp_P2C_3dArray ( qy_ypencil, d, accc_ypencil )   ! intermediate, accc_ypencil = qy_yccc_ypencil
    call transpose_y_to_x(accc_ypencil, qy_yccc, d%dccc)          ! x-pencil : x-mom, w   thermal
    call transpose_y_to_z(accc_ypencil, qy_yccc_zpencil, d%dccc)  ! z-pencil : z-mom, w   thermal
    end if
!-------------------------------------------------------------------------------
!    qz --> qz_xpcp(WO) --> qz_xpcp_ypencil(I) --> qz_xpcp_zpencil(O)
!     | --> qz_ypencil(WO) --> qz_ycpp_ypencil(WO) --> qz_ycpp_zpencil(O)
!                       |  --> qz_zpencil(WO) --> qz_zccc_zpencil(I) --> qz_zccc_ypencil(W) --> qz_zccc(W)
!------------------------------------------------------------------------------- 
    call Get_x_midp_C2P_3dArray ( f%qz, d, qz_xpcp )                  ! x-pencil : z-mom, w+o   thermal
    if(ithermo == 0) then
    call transpose_x_to_y(qz_xpcp,         apcp_ypencil, d%dpcp)      ! intermediate, apcp_ypencil = qz_xpcp_ypencil
    call transpose_y_to_z(apcp_ypencil, qz_xpcp_zpencil, d%dpcp)      ! z-pencil : x-mom, o   thermal
    end if
  
    call transpose_x_to_y(f%qz, qz_ypencil, d%dccp)                   ! y-pencil : z-mom, w+o   thermal
    call Get_y_midp_C2P_3dArray ( qz_ypencil, d, qz_ycpp_ypencil )    ! y-pencil : z-mom, w+o   thermal
    if(ithermo == 0) &
    call transpose_y_to_z ( qz_ycpp_ypencil, qz_ycpp_zpencil, d%dcpp )! z-pencil : y-mom, o   thermal
  
    call transpose_y_to_z(qz_ypencil, qz_zpencil, d%dccp)             ! z-pencil : z-mom, w+o   thermal
    if(ithermo == 1) then
    call Get_z_midp_P2C_3dArray ( qz_zpencil, d, accc_zpencil )       ! intermediate, accc_zpencil = qz_zccc_zpencil
    call transpose_z_to_y(   accc_zpencil, qz_zccc_ypencil, d%dccc)   ! y-pencil : y-mom, w   thermal
    call transpose_y_to_x(qz_zccc_ypencil, qz_zccc,         d%dccc)   ! x-pencil : x-mom, w   thermal
    end if

    if(ithermo == 1) then
!-------------------------------------------------------------------------------
!    gx --> gx_ypencil(I) --> gx_yppc_ypencil(I)--> gx_yppc(W)
!                     |--> gx_zpencil(I) --> gx_zpcp_zpencil(I) --> gx_zpcp_ypencil(I) --> qx_zpcp
!------------------------------------------------------------------------------- 
    call transpose_x_to_y(f%gx, apcc_ypencil, d%dpcc)               ! intermediate, apcc_ypencil = gx_ypencil
    call Get_y_midp_C2P_3dArray ( apcc_ypencil, d, appc_ypencil )   ! intermediate, appc_ypencil = gx_yppc_ypencil
    call transpose_y_to_x(appc_ypencil, gx_yppc, d%dppc)            ! x-pencil : y-mom, w   thermal
  
    call transpose_y_to_z(apcc_ypencil, apcc_zpencil, d%dpcc)       ! intermediate, apcc_zpencil = gx_zpencil
    call Get_z_midp_C2P_3dArray ( apcc_zpencil, d, apcp_zpencil )   ! intermediate, apcp_zpencil = gx_zpcp_zpencil
    call transpose_z_to_y(apcp_zpencil, apcp_ypencil, d%dpcp)       ! intermediate, apcp_ypencil = gx_zpcp_ypencil
    call transpose_y_to_x(apcp_ypencil,      gx_zpcp, d%dpcp)       ! x-pencil : z-mom, wo  thermal
!-------------------------------------------------------------------------------
!    gy --> gy_ypencil(W) --> gy_zpencil(I) --> gy_zcpp_zpencil(I) --> gy_zcpp_ypencil(W)
!-------------------------------------------------------------------------------
    call transpose_x_to_y(f%gy,   gy_ypencil, d%dcpc)               ! y-pencil : y-mom, w   thermal
    call transpose_y_to_z(gy_ypencil, acpc_zpencil, d%dcpc)         ! intermediate, acpc_zpencil = gy_zpencil
    call Get_z_midp_C2P_3dArray ( acpc_zpencil, d, acpp_zpencil )   ! intermediate, acpp_zpencil = gy_zcpp_zpencil
    call transpose_z_to_y(acpp_zpencil, gy_zcpp_ypencil, d%dcpp)    ! y-pencil : z-mom, w   thermal
!-------------------------------------------------------------------------------
!    gz --> gz_xpcp(I)    --> gz_xpcp_ypencil(I) --> gz_xpcp_zpencil(W)
!     | --> gz_ypencil(I) --> gz_ycpp_ypencil(I) --> gz_ycpp_zpencil(W)
!------------------------------------------------------------------------------- 
    call Get_x_midp_C2P_3dArray ( f%gz, d, apcp )                   ! intermediate, apcp = gz_xpcp
    call transpose_x_to_y(apcp,            apcp_ypencil, d%dpcp)    ! intermediate  apcp_ypencil = gz_xpcp_ypencil
    call transpose_y_to_z(apcp_ypencil, gz_xpcp_zpencil, d%dpcp)    ! z-pencil : x-mom, w   thermal
  
    call transpose_x_to_y(f%gz,   accp, d%dccp)                     ! intermediate, accp = gz_ypencil
    call Get_y_midp_C2P_3dArray ( accp, d, acpp_ypencil )           ! intermediate, acpp_ypencil = gz_ycpp_ypencil
    call transpose_y_to_z(acpp_ypencil, gz_ycpp_zpencil, d%dcpp)    ! z-pencil : y-mom, w   thermal
!-------------------------------------------------------------------------------
!   d --> d_ypencil --> d_zpencil
!-------------------------------------------------------------------------------
    call transpose_x_to_y(f%dDens,       dDens_ypencil, d%dccc)     ! y-pencil : y-mom, w   thermal
    call transpose_y_to_z(dDens_ypencil, dDens_zpencil, d%dccc)     ! z-pencil : z-mom, w   thermal
!-------------------------------------------------------------------------------
!    m --> dmdx_xpcc
!    | --> m_xpcc -->m_xpcc_ypencil -->dmdy_xpcc_ypencil-->dmdy_xpcc
!                                 | -->m_xpcc_zpencil --> dmdz_xpcc_zpencil--> dmdz_xpcc_ypencil(I) --> dmdz_xpcc 
!    | --> m_ypencil(I) --> m_ycpc_ypencil --> m_ycpc --> dmdx_ycpc
!                                        | --> m_ycpc_zpencil
!                     | --> m_zpencil(I) --> dmdz_zccp_zpencil
!                                      | --> m_zccp_zpencil --> m_zccp_ypencil --> m_zccp --> dmdx_zccp --> dmdx_zccp_ypencil(I) --> dmdx_zccp_zpencil--> dmdy_zccp_ypencil
!------------------------------------------------------------------------------- 
    call Get_x_1st_derivative_C2P_3dArray( f%mVisc, d, dmdx_xpcc )                ! x-pencil : x-mom, w thermal
    call Get_x_midp_C2P_3dArray ( f%mVisc, d,  m_xpcc )                           ! x-pencil : x-mom, w thermal
    call transpose_x_to_y(m_xpcc, m_xpcc_ypencil, d%dpcc)                         ! y-pencil : x-mom, w thermal

    call Get_y_1st_derivative_C2C_3dArray( m_xpcc_ypencil, d, dmdy_xpcc_ypencil ) ! y-pencil : x-mom, w thermal
    call transpose_y_to_x(dmdy_xpcc_ypencil, dmdy_xpcc, d%dpcc)                   ! x-pencil : x-mom, w thermal

    call transpose_y_to_z(m_xpcc_ypencil, m_xpcc_zpencil, d%dpcc)                 ! z-pencil : x-mom, w thermal
    call Get_z_1st_derivative_C2C_3dArray( m_xpcc_zpencil, d, dmdz_xpcc_zpencil ) ! z-pencil : x-mom, w thermal
    call transpose_z_to_y(dmdz_xpcc_zpencil, apcc_ypencil, d%dpcc)                ! intermediate, apcc_ypencil = dmdz_xpcc_ypencil
    call transpose_y_to_x(apcc_ypencil, dmdz_xpcc,         d%dpcc)                ! x-pencil : x-mom, w thermal

    call transpose_x_to_y(f%mVisc, accc_ypencil, d%dccc)                          ! intermediate, accc_ypencil = m_ypencil 
    call Get_y_midp_C2P_3dArray ( accc_ypencil, d,  m_ycpc_ypencil )              ! y-pencil : y-mom, w thermal
    call transpose_y_to_z(m_ycpc_ypencil, m_ycpc_zpencil, d%dcpc)                 ! z-pencil : y-mom, w thermal
    call transpose_y_to_x(m_ycpc_ypencil, m_ycpc,         d%dcpc)                 ! x-pencil : y-mom, w thermal
    call Get_x_1st_derivative_C2C_3dArray(  m_ycpc, d, dmdx_ycpc )                ! x-pencil : y-mom, w thermal
      
    call transpose_y_to_z(accc_ypencil,  accc_zpencil, d%dccc)                    ! intermediate, accc_zpencil = m_zpencil
    call Get_z_1st_derivative_C2P_3dArray( accc_zpencil, d, dmdz_zccp_zpencil )   ! z-pencil : z-mom, w thermal
    call Get_z_midp_C2P_3dArray ( accc_zpencil, d,  m_zccp_zpencil )              ! z-pencil : z-mom, w thermal
    call transpose_z_to_y(m_zccp_zpencil, m_zccp_ypencil, d%dccp)                 ! y-pencil : z-mom, w thermal
    call transpose_y_to_x(m_zccp_ypencil, m_zccp,         d%dccp)                 ! x-pencil : z-mom, w thermal
    call Get_x_1st_derivative_C2C_3dArray(  m_zccp, d, dmdx_zccp )                ! x-pencil : z-mom, w thermal
    call transpose_x_to_y(dmdx_zccp,         accp_ypencil, d%dccp)                ! intermidate, accp_ypencil = dmdx_zccp_ypencil
    call transpose_y_to_z(accp_ypencil, dmdx_zccp_zpencil, d%dccp)                ! z-pencil : z-mom, w thermal
    call Get_y_1st_derivative_C2C_3dArray( m_zccp_ypencil, d, dmdy_zccp_ypencil ) ! y-pencil : z-mom, w thermal
!-------------------------------------------------------------------------------
! calculate div(u_vec)
!-------------------------------------------------------------------------------
    div  = ZERO 
    accc = ZERO

    call Get_x_1st_derivative_P2C_3dArray( f%qx, d, accc ) ! accc = d(qx)/d(x)_xccc
    div = div + accc ! = d(qx)/d(x)_ccc

    call Get_y_1st_derivative_P2C_3dArray( qy_ypencil, d, accc_ypencil ) ! accc_ypencil = d(qy)/(y)_yccc_ypencil
    call transpose_y_to_x(accc_ypencil,       accc,      d%dccc)       ! accc = d(qy)/d(y)_yccc
    div = div + accc ! = d(qx)/d(x)_ccc + d(qy)/d(y)_ccc

    call Get_z_1st_derivative_P2C_3dArray( qz_zpencil, d,  accc_zpencil ) ! accc_zpencil = d(qz)/(z)_zccc_zpencil
    call transpose_z_to_y(accc_zpencil, accc_ypencil, d%dccc)           ! accc_ypencil = d(qz)/(z)_zccc_ypencil
    call transpose_y_to_x(accc_ypencil, accc,         d%dccc)           ! accc = d(qz)/d(z)_zccc
    div = div + accc ! = d(qx)/d(x)_ccc + d(qy)/d(y)_ccc + d(qz)/d(z)_ccc

    call transpose_x_to_y(div,         div_ypencil, d%dccc)
    call transpose_y_to_z(div_ypencil, div_zpencil, d%dccc)
    end if
!===============================================================================
! the RHS of momentum equation
! x-pencil : the RHS terms of all 3 momentum equations (derivative) operating in the x direction
!===============================================================================
    f%mx_rhs = ZERO
    f%my_rhs = ZERO
    f%mz_rhs = ZERO
    mx_rhs_implicit  = ZERO

    apcc = ZERO
    acpc = ZERO
    accp = ZERO
!-------------------------------------------------------------------------------
! X-pencil : X-mom convection term (x-c1/3): d(gx * qx)/dx at (i', j, k)
!-------------------------------------------------------------------------------
    if(ithermo == 0) call Get_x_1st_derivative_P2P_3dArray( -f%qx * f%qx, d, apcc )
    if(ithermo == 1) call Get_x_1st_derivative_P2P_3dArray( -f%gx * f%qx, d, apcc )
    f%mx_rhs = f%mx_rhs + apcc
!-------------------------------------------------------------------------------
! X-pencil : X-mom diffusion term (x-v1-1/7), \mu^x * L11(ux) at (i', j, k)
!-------------------------------------------------------------------------------
    call Get_x_2nd_derivative_P2P_3dArray( f%qx, d, apcc )
    if(ithermo == 0) f%mx_rhs = f%mx_rhs +          f%rre * apcc
    if(ithermo == 1) f%mx_rhs = f%mx_rhs + m_xpcc * f%rre * apcc
!-------------------------------------------------------------------------------
! X-pencil : Y-mom convection term (y-c1/3), d(gx^y * qy^x)/dx at (i, j', k)
!-------------------------------------------------------------------------------
    if(ithermo == 0) call Get_x_1st_derivative_P2C_3dArray( -qx_yppc * qy_xppc, d, acpc )
    if(ithermo == 1) call Get_x_1st_derivative_P2C_3dArray( -gx_yppc * qy_xppc, d, acpc )
    f%my_rhs = f%my_rhs + acpc
!-------------------------------------------------------------------------------
! X-pencil : Y-mom diffusion term (y-v1-1/1), \mu * L11(uy) at (i, j', k)
!-------------------------------------------------------------------------------
    call Get_x_2nd_derivative_P2P_3dArray( f%qy, d, acpc )
    if(ithermo == 0) f%my_rhs = f%my_rhs +          f%rre * acpc
    if(ithermo == 1) f%my_rhs = f%my_rhs + m_ycpc * f%rre * acpc
!-------------------------------------------------------------------------------
! X-pencil : Z-mom convection term (z-c1/3), d(gx^z * qz^x)/dx at (i, j, k')
!-------------------------------------------------------------------------------
    if(ithermo == 0) call Get_x_1st_derivative_P2C_3dArray( -qx_zpcp * qz_xpcp, d, accp )
    if(ithermo == 1) call Get_x_1st_derivative_P2C_3dArray( -gx_zpcp * qz_xpcp, d, accp )
    f%mz_rhs = f%mz_rhs + accp
!-------------------------------------------------------------------------------
! X-pencil : Z-mom diffusion term (z-v1-1/1), \mu * L11(uz) at (i, j, k')
!-------------------------------------------------------------------------------
    call Get_x_2nd_derivative_P2P_3dArray( f%qz, d, accp )
    if(ithermo == 0) f%mz_rhs = f%mz_rhs +          f%rre * accp
    if(ithermo == 1) f%mz_rhs = f%mz_rhs + m_zccp * f%rre * accp
!-------------------------------------------------------------------------------
! X-pencil : X-mom, pressure gradient in x direction, sigma_1*d(p)/dx
!-------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3dArray( f%pres, d, apcc )
    mx_rhs_implicit =  mx_rhs_implicit - sigma1p * apcc

    if(ithermo == 1) then
!-------------------------------------------------------------------------------
! x-pencil : X-mom, gravity force in x direction
!-------------------------------------------------------------------------------
      if(igravity == 1 .or. igravity == -1)  then
        call Get_x_midp_C2P_3dArray( f%dDens, d, apcc )
        mx_rhs_implicit =  mx_rhs_implicit + f%fgravity * apcc
      end if
!-------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v2/7), \mu^x * 1/3 * d (div)/dx at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3dArray( div, d, apcc )
      f%mx_rhs = f%mx_rhs + one_third_rre * m_xpcc * apcc
!-------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v3/7), -2/3 * d\mu/dx * (div(u)^x) +  
!                                                    2 * d\mu/dx * du/dx
!-------------------------------------------------------------------------------
      call Get_x_midp_C2P_3dArray          ( div,  d, apcc )
      f%mx_rhs = f%mx_rhs - two_third_rre * dmdx_xpcc * apcc
      call Get_x_1st_derivative_P2P_3dArray( f%qx, d, apcc )
      f%mx_rhs = f%mx_rhs + two_rre       * dmdx_xpcc * apcc
!-------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v4/7), d(mu^x)/dy * d(qy^y)/dx at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3dArray( qy_yccc, d, apcc )
      f%mx_rhs =  f%mx_rhs + f%rre * dmdy_xpcc * apcc
!-------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v5/7), d(mu^x)/dz * d(qz^z)/dx at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3dArray( qz_zccc, d, apcc )
      f%mx_rhs =  f%mx_rhs + f%rre * dmdz_xpcc * apcc
!-------------------------------------------------------------------------------
!   X-pencil : Y-mom diffusion term (y-v6/7), d(mu^y)/dx * d(qy^x))/dx at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2C_3dArray( f%qy, d, acpc )
      f%my_rhs =  f%my_rhs + f%rre * dmdx_ycpc * acpc
!-------------------------------------------------------------------------------
!   X-pencil : Z-mom diffusion term (z-v6/7), d(mu^z)/dx * d(qz)/dx at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2C_3dArray( f%qz, d, accp )
      f%mz_rhs =  f%mz_rhs + f%rre * dmdx_zccp * accp
    end if

!===============================================================================
! the RHS of momentum equation
! Y-pencil : the RHS terms of all 3 momentum equations (derivative) operating in the y direction
!===============================================================================
    mx_rhs_ypencil = ZERO
    my_rhs_ypencil = ZERO
    mz_rhs_ypencil = ZERO
    my_rhs_implicit_ypencil = ZERO
    apcc_ypencil = ZERO
    acpc_ypencil = ZERO
    accp_ypencil = ZERO
!-------------------------------------------------------------------------------
! Y-pencil : X-mom convection term (x-c2/3): d(<gy>^x * <qx>^y)/dy at (i', j, k)
!-------------------------------------------------------------------------------
    if(ithermo == 0) call Get_y_1st_derivative_P2C_3dArray( -qy_xppc_ypencil * qx_yppc_ypencil,  d, apcc_ypencil )
    if(ithermo == 1) call Get_y_1st_derivative_P2C_3dArray( -gy_xppc_ypencil * qx_yppc_ypencil,  d, apcc_ypencil )
    mx_rhs_ypencil = mx_rhs_ypencil + apcc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : X-mom diffusion term (x-v1-2/7), \mu^x * L22(ux) at (i', j, k)
!-------------------------------------------------------------------------------
    call Get_y_2nd_derivative_P2P_3dArray( qx_ypencil, d, apcc_ypencil )
    if(ithermo == 0) mx_rhs_ypencil = mx_rhs_ypencil +                  f%rre * apcc_ypencil
    if(ithermo == 1) mx_rhs_ypencil = mx_rhs_ypencil + m_xpcc_ypencil * f%rre * apcc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom convection term (y-c2/3), d(gy * qy)/dy at (i, j', k)
!-------------------------------------------------------------------------------
    if(ithermo == 0) call Get_y_1st_derivative_P2P_3dArray( -qy_ypencil * qy_ypencil,  d, acpc_ypencil )
    if(ithermo == 1) call Get_y_1st_derivative_P2P_3dArray( -gy_ypencil * qy_ypencil,  d, acpc_ypencil )
    my_rhs_ypencil = my_rhs_ypencil + acpc_ypencil

!-------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v1-2/1), \mu * L22(uy) at (i, j', k)
!-------------------------------------------------------------------------------
    call Get_y_2nd_derivative_P2P_3dArray( qy_ypencil,  d, acpc_ypencil )
    if(ithermo == 0) my_rhs_ypencil = my_rhs_ypencil +                  f%rre * acpc_ypencil
    if(ithermo == 1) my_rhs_ypencil = my_rhs_ypencil + m_ycpc_ypencil * f%rre * acpc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Z-mom convection term (z-c2/3), d(gy^z * qz^y)/dy at (i, j, k')
!-------------------------------------------------------------------------------
    if(ithermo == 0) call Get_y_1st_derivative_P2C_3dArray( -qy_zcpp_ypencil * qz_ycpp_ypencil,  d, accp_ypencil )
    if(ithermo == 1) call Get_y_1st_derivative_P2C_3dArray( -gy_zcpp_ypencil * qz_ycpp_ypencil,  d, accp_ypencil )
    mz_rhs_ypencil = mz_rhs_ypencil + accp_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Z-mom diffusion term (z-v1-2/1), \mu * L22(uz) at (i, j, k')
!-------------------------------------------------------------------------------
    call Get_y_2nd_derivative_P2P_3dArray( qz_ypencil,  d, accp_ypencil )
    if(ithermo == 0) mz_rhs_ypencil = mz_rhs_ypencil +                  f%rre * accp_ypencil
    if(ithermo == 1) mz_rhs_ypencil = mz_rhs_ypencil + m_zccp_ypencil * f%rre * accp_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom pressure gradient in y direction, d(sigma_1 p)
!-------------------------------------------------------------------------------
    call Get_y_1st_derivative_C2P_3dArray( pres_ypencil, d, acpc_ypencil )
    my_rhs_implicit_ypencil =  my_rhs_implicit_ypencil - sigma1p * acpc_ypencil

    if(ithermo == 1) then
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom gravity force in y direction
!-------------------------------------------------------------------------------
      if( igravity == 2 .or. igravity == -2) then
        call Get_y_midp_C2P_3dArray( dDens_ypencil, d, acpc_ypencil )
        my_rhs_implicit_ypencil =  my_rhs_implicit_ypencil + f%fgravity * acpc_ypencil
      end if
!-------------------------------------------------------------------------------
! Y-pencil : X-mom diffusion term (x-v6/7), d(mu^x)/dy * d(qx)/dy at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2C_3dArray( qx_ypencil, d, apcc_ypencil )
      mx_rhs_ypencil =  mx_rhs_ypencil + f%rre * dmdy_xpcc_ypencil * apcc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v2/7), \mu^y * 1/3 * d (div)/dy at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2P_3dArray( div_ypencil, d, acpc_ypencil )
      my_rhs_ypencil = my_rhs_ypencil + one_third_rre * m_ycpc_ypencil * acpc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v3/7), 2d\mu/dy * (-1/3 * div(u)) +  2d\mu/dy * dv/dy
!-------------------------------------------------------------------------------
      call Get_y_midp_C2P_3dArray          ( div_ypencil, d, acpc_ypencil )
      my_rhs_ypencil = my_rhs_ypencil - two_third_rre * dmdy_ycpc_ypencil * acpc_ypencil
      call Get_y_1st_derivative_P2P_3dArray( qy_ypencil, d, acpc_ypencil )
      my_rhs_ypencil = my_rhs_ypencil + two_rre       * dmdy_ycpc_ypencil * acpc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v4/7), d(mu^y)/dx * d(qx^x)/dy at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2P_3dArray( qx_xccc_ypencil, d, acpc_ypencil )
      my_rhs_ypencil =  my_rhs_ypencil + f%rre * dmdx_ycpc_ypencil * acpc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v6/7), d(mu)/dz * d(qz^z)/dy at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2P_3dArray( qz_zccc_ypencil, d, acpc_ypencil )
      my_rhs_ypencil =  my_rhs_ypencil + f%rre * dmdz_ycpc_ypencil * acpc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Z-mom diffusion term (z-v7/7), d(mu^z)/dy * d(qz)/dy at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2C_3dArray( qz_ypencil, d, accp_ypencil )
      mz_rhs_ypencil =  mz_rhs_ypencil + f%rre * dmdy_zccp_ypencil * accp_ypencil
    end if
!-------------------------------------------------------------------------------
! Data from Y-pencil to X-pencil
!-------------------------------------------------------------------------------    
    call transpose_from_y_to_x(mx_rhs_ypencil, apcc, d%dpcc)
    call transpose_from_y_to_x(my_rhs_ypencil, acpc, d%dcpc)
    call transpose_from_y_to_x(mz_rhs_ypencil, accp, d%dccp)
    f%mx_rhs = f%mx_rhs + apcc
    f%my_rhs = f%my_rhs + acpc
    f%mz_rhs = f%mz_rhs + accp
    call transpose_from_y_to_x(my_rhs_implicit_ypencil, acpc, d%dcpc)
    my_rhs_implicit = my_rhs_implicit + acpc

!===============================================================================
! the RHS of momentum equation
! z-pencil : the RHS terms of all 3 momentum equations (derivative) operating in the z direction
!===============================================================================
    mx_rhs_zpencil = ZERO
    my_rhs_zpencil = ZERO
    mz_rhs_zpencil = ZERO
    mz_rhs_implicit_zpencil = ZERO
    apcc_zpencil = ZERO
    acpc_zpencil = ZERO
    accp_zpencil = ZERO
!-------------------------------------------------------------------------------
! Z-pencil : X-mom convection term (x-c3/3): d(<gz>^x * <qx>^z)/dz at (i', j, k)
!-------------------------------------------------------------------------------
    if(ithermo == 0) call Get_z_1st_derivative_P2C_3dArray( -qz_xpcp_zpencil * qx_zpcp_zpencil,  d, apcc_zpencil )
    if(ithermo == 1) call Get_z_1st_derivative_P2C_3dArray( -gz_xpcp_zpencil * qx_zpcp_zpencil,  d, apcc_zpencil )
    mx_rhs_zpencil = mx_rhs_zpencil + apcc_zpencil

!-------------------------------------------------------------------------------
! Z-pencil : X-mom diffusion term (x-v1-3/7), \mu^x * L33(ux) at (i', j, k)
!-------------------------------------------------------------------------------
    call Get_z_2nd_derivative_P2P_3dArray( qx_zpencil, d, apcc_zpencil )
    if(ithermo == 0) mx_rhs_zpencil = mx_rhs_zpencil +                  f%rre * apcc_zpencil
    if(ithermo == 1) mx_rhs_zpencil = mx_rhs_zpencil + m_xpcc_zpencil * f%rre * apcc_zpencil

!-------------------------------------------------------------------------------
! Z-pencil : Y-mom convection term (y-c3/3), d(<gz>^y * <qy>^z)/dz at (i, j', k)
!-------------------------------------------------------------------------------
    if(ithermo == 0) call Get_z_1st_derivative_P2C_3dArray( -qz_ycpp_zpencil * qy_zcpp_zpencil,  d, acpc_zpencil )
    if(ithermo == 1) call Get_z_1st_derivative_P2C_3dArray( -gz_ycpp_zpencil * qy_zcpp_zpencil,  d, acpc_zpencil )
    my_rhs_zpencil = my_rhs_zpencil + acpc_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Y-mom diffusion term (y-v1-3/7), \mu * L33(uy) at (i, j', k)
!-------------------------------------------------------------------------------
    call Get_z_2nd_derivative_P2P_3dArray( qy_zpencil,  d, acpc_zpencil )
    if(ithermo == 0) my_rhs_zpencil = my_rhs_zpencil +                  f%rre * acpc_zpencil
    if(ithermo == 1) my_rhs_zpencil = my_rhs_zpencil + m_ycpc_zpencil * f%rre * acpc_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v1/1), \mu * L33(uz) at (i, j, k')
!-------------------------------------------------------------------------------
    call Get_z_2nd_derivative_P2P_3dArray( qz_zpencil,  d, accp_zpencil )
    if(ithermo == 0) mz_rhs_zpencil = mz_rhs_zpencil +                  f%rre * accp_zpencil
    if(ithermo == 1) mz_rhs_zpencil = mz_rhs_zpencil + m_zccp_zpencil * f%rre * accp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : pressure gradient in z direction, d(sigma_1 p)
!-------------------------------------------------------------------------------
    call Get_z_1st_derivative_C2P_3dArray( pres_zpencil, d, accp_zpencil )
    mz_rhs_implicit_zpencil =  mz_rhs_implicit_zpencil - sigma1p * accp_zpencil

    if(ithermo == 1) then
!-------------------------------------------------------------------------------
! z-pencil : gravity force in z direction
!-------------------------------------------------------------------------------
      if( igravity == 3 .or. igravity == -3) then
        call Get_z_midp_C2P_3dArray( dDens_zpencil, d, accp_zpencil )
        mz_rhs_implicit_zpencil =  mz_rhs_implicit_zpencil + f%fgravity * accp_zpencil
      end if
!-------------------------------------------------------------------------------
! Z-pencil : X-mom diffusion term (x-v7/7), d(mu^x)/dz * d(qx)/dz at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2C_3dArray( qx_zpencil, d, apcc_zpencil )
      mx_rhs_zpencil = mx_rhs_zpencil + f%rre * dmdz_xpcc_zpencil * apcc_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Y-mom diffusion term (y-v7/7), d(mu^y)/dz * d(qy)/dz at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2C_3dArray( qy_zpencil, d, acpc_zpencil )
      my_rhs_zpencil =  my_rhs_zpencil + f%rre * dmdz_ycpc_zpencil * acpc_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v2/7), \mu * 1/3 * d (div)/dz at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2P_3dArray( div_zpencil, d, accp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil + one_third_rre * m_zccp_zpencil * accp_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v3/7), 2d\mu/dz * (-1/3 * div(u)) +  2d\mu/dz * dw/dz
!-------------------------------------------------------------------------------
      call Get_z_midp_C2P_3dArray          ( div_zpencil, d, accp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil - two_third_rre * dmdz_zccp_zpencil * accp_zpencil
      call Get_z_1st_derivative_P2P_3dArray( qz_zpencil, d, accp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil + two_rre * dmdz_zccp_zpencil * accp_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v4/7), d(mu^z)/dx * d(qx^x)/dz at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2P_3dArray( qx_xccc_zpencil, d, accp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil + f%rre * dmdx_zccp_zpencil * accp_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v6/7), d(mu^z)/dy * d(qy^y)/dz at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2P_3dArray( qy_yccc_zpencil, d, accp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil + f%rre * dmdy_zccp_zpencil * accp_zpencil

    end if
!-------------------------------------------------------------------------------
! Data from Z-pencil to X-pencil
!-------------------------------------------------------------------------------    
    call transpose_from_z_to_y(mx_rhs_zpencil, apcc_ypencil, d%dpcc)
    call transpose_from_z_to_y(my_rhs_zpencil, acpc_ypencil, d%dcpc)
    call transpose_from_z_to_y(mz_rhs_zpencil, accp_ypencil, d%dccp)
    
    call transpose_from_y_to_x(apcc_ypencil, apcc, d%dpcc)
    call transpose_from_y_to_x(acpc_ypencil, acpc, d%dcpc)
    call transpose_from_y_to_x(accp_ypencil, accp, d%dccp)

    f%mx_rhs = f%mx_rhs + apcc
    f%my_rhs = f%my_rhs + acpc
    f%mz_rhs = f%mz_rhs + accp

    call transpose_from_z_to_y(mz_rhs_implicit_zpencil, accp_ypencil, d%dccp)
    call transpose_from_y_to_x(accp_ypencil, accp, d%dccp)

    mz_rhs_implicit = mz_rhs_implicit + accp
!===============================================================================
! x-pencil : to build up rhs in total, in all directions
!===============================================================================
!-------------------------------------------------------------------------------
! x-pencil : x-momentum
!-------------------------------------------------------------------------------
    call Calculate_momentum_fractional_step(f%mx_rhs0, f%mx_rhs, mx_rhs_implicit, isub)
    if(idriven /= IDRVF_NO) call Calculate_momentum_driven_source(f%mx_rhs, 'ux', d, isub) 
!-------------------------------------------------------------------------------
! x-pencil : y-momentum
!-------------------------------------------------------------------------------
    call Calculate_momentum_fractional_step(f%my_rhs0, f%my_rhs, my_rhs_implicit, isub)
!-------------------------------------------------------------------------------
! x-pencil : z-momentum
!-------------------------------------------------------------------------------
    call Calculate_momentum_fractional_step(f%mz_rhs0, f%mz_rhs, mz_rhs_implicit, isub)
 
    return
  end subroutine Compute_momentum_rhs

!===============================================================================
!===============================================================================
!> \brief To update the provisional u or rho u.
!>
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                         !
!______________________________________________________________________________!
!> \param[inout]  rhs          the rhs
!> \param[inout]  u            provisional u or rho u.
!_______________________________________________________________________________
  subroutine Calculate_intermediate_mvar(rhs, u)
    implicit none
    real(WP), dimension(:, :, :), intent(inout) :: rhs, u

    u(:, :, :) = u(:, :, :) + rhs(:, :, :)

    return
  end subroutine Calculate_intermediate_mvar
!===============================================================================
!===============================================================================
  subroutine Correct_massflux(ux, uy, uz, phi, d, isub)
    use udf_type_mod,      only : t_domain
    use input_general_mod, only : tAlpha, dt, sigma2p
    use operations
    implicit none

    type(t_domain), intent(in   ) :: d
    integer(4),     intent(in   ) :: isub
    real(WP), dimension( d%np(1), d%nc(2), d%nc(3) ), intent(inout) :: ux
    real(WP), dimension( d%nc(1), d%np(2), d%nc(3) ), intent(inout) :: uy
    real(WP), dimension( d%nc(1), d%nc(2), d%np(3) ), intent(inout) :: uz
    real(WP), dimension( d%nc(1), d%nc(2), d%nc(3) ), intent(in   ) :: phi

    real(WP), dimension( d%np(1), d%nc(2), d%nc(3) ) :: dphidx
    real(WP), dimension( d%nc(1), d%np(2), d%nc(3) ) :: dphidy
    real(WP), dimension( d%nc(1), d%nc(2), d%np(3) ) :: dphidz

  
    call Get_x_1st_derivative_C2P_3dArray( phi,  dphidx )
    ux = ux - dt * tAlpha(isub) * sigma2p * dphidx

    call Get_y_1st_derivative_C2P_3dArray( phi,  dphidy )
    uy = uy - dt * tAlpha(isub) * sigma2p * dphidy

    call Get_z_1st_derivative_C2P_3dArray( phi,  dphidz )
    uz = uz - dt * tAlpha(isub) * sigma2p * dphidz

    return
  end subroutine Correct_massflux
!===============================================================================
!===============================================================================
!> \brief To update the provisional u or rho u.
!>
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                         !
!______________________________________________________________________________!
!> \param[inout]  f            flow field
!> \param[inout]  d            domain
!> \param[in]     isub         RK sub-iteration
!_______________________________________________________________________________
  subroutine Solve_momentum_eq(f, d, isub)
    use input_general_mod, only : iviscous, IVIS_SEMIMPLT, IVIS_EXPLICIT, ithermo
    use udf_type_mod,      only : t_flow, t_domain
    use typeconvert_mod
    use continuity_eq_mod
    use poisson_mod
    use boundary_conditions_mod
    implicit none

    type(t_flow), intent(inout) :: f
    type(t_domain),  intent(in) :: d
    integer(4),      intent(in) :: isub

!-------------------------------------------------------------------------------
! to calculate the rhs of the momenturn equation in stepping method
!_______________________________________________________________________________ 
    call Compute_momentum_rhs(f, d, isub)
!-------------------------------------------------------------------------------
! to update intermediate (\hat{q}) or (\hat{g})
!_______________________________________________________________________________
 
    if(iviscous == IVIS_EXPLICIT) then

      if(ithermo == 0) then 
        call Calculate_intermediate_mvar(f%mx_rhs, f%qx)
        call Calculate_intermediate_mvar(f%my_rhs, f%qy)
        call Calculate_intermediate_mvar(f%mz_rhs, f%qz)
      else
        call Calculate_intermediate_mvar(f%mx_rhs, f%gx)
        call Calculate_intermediate_mvar(f%my_rhs, f%gy)
        call Calculate_intermediate_mvar(f%mz_rhs, f%gz)
      end if

    else if(iviscous == IVIS_SEMIMPLT) then
    !in order for a high order spacial accuracy
    ! to use Alternating direction implicit method
    ! ref: Cui2013: Convergence analysis of high-order compact 
    ! alternating direction implicit schemes for the two-dimensional 
    ! time fractional equation
      stop
    else 

    end if
!-------------------------------------------------------------------------------
! to update b.c. values
!-------------------------------------------------------------------------------
    call Apply_BC_velocity (f%gx, f%gy, f%gz, d)
    call Apply_BC_velocity (f%qx, f%qy, f%qz, d)
!-------------------------------------------------------------------------------
! to calculate the provisional divergence constrains
!-------------------------------------------------------------------------------
  !call Print_debug_mid_msg("  Computing provisional divergence constrains ...")
    call Calculate_continuity_constrains(f, d, isub)
!-------------------------------------------------------------------------------
! to solve Poisson equation
!-------------------------------------------------------------------------------
  !call Print_debug_mid_msg("  Solving Poisson Equation ...")
    call Solve_poisson(f%pcor)
!-------------------------------------------------------------------------------
! to update velocity/massflux correction
!-------------------------------------------------------------------------------
  !call Print_debug_mid_msg("  Updating velocity/mass flux ...")
    if(ithermo == 0) then 
      call Correct_massflux(f%qx, f%qy, f%qz, f%pcor, d, isub)
    else
      call Correct_massflux(f%gx, f%gy, f%gz, f%pcor, d, isub)
    end if
!-------------------------------------------------------------------------------
! to update pressure
!-------------------------------------------------------------------------------
    f%pres = f%pres + f%pcor
!-------------------------------------------------------------------------------
! to update b.c. values
!-------------------------------------------------------------------------------
    call Apply_BC_velocity (f%gx, f%gy, f%gz, d)
    call Apply_BC_velocity (f%qx, f%qy, f%qz, d)

    return
  end subroutine

end module eq_momentum_mod
