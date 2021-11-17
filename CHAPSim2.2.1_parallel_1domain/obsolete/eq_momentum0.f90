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
!  mode           name          role                                           !
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
!  mode           name          role                                           !
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
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  f             flow field
!> \param[inout]  d             domain    
!> \param[in]     isub          the RK iteration to get correct Coefficient 
!_______________________________________________________________________________
  subroutine Compute_momentum_rhs(f, d, isub)
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
! ithermo == 0, x-pencil
!-------------------------------------------------------------------------------
    
    real(WP), dimension( d%ux_xsz(1), d%uy_xsz(2), d%uy_xsz(3) ) :: qy_xppc ! <uy>^x at (xp, yp, zc)
    real(WP), dimension( d%ux_xsz(1), d%uz_xsz(2), d%uz_xsz(3) ) :: qz_xpcp ! <uz>^x at (xp, yc, zp)

    real(WP), dimension( d%ux_xsz(1), d%uy_xsz(2), d%ux_xsz(3) ) :: qx_yppc ! <ux>^y at (xp, yp, zc)
    real(WP), dimension( d%ux_xsz(1), d%ux_xsz(2), d%uz_xsz(3) ) :: qx_zpcp ! <ux>^z at (xp, yc, zp)
    
    
!-------------------------------------------------------------------------------
! ithermo == 0, y-pencil
!-------------------------------------------------------------------------------
    real(WP), dimension( d%ux_ysz(1), d%ux_ysz(2), d%ux_ysz(3) ) ::      qx_ypencil
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) ::      qy_ypencil
    real(WP), dimension( d%uz_ysz(1), d%uz_ysz(2), d%uz_ysz(3) ) ::      qz_ypencil

    real(WP), dimension( d%ux_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) :: qy_xppc_ypencil ! <uy>^x at (xp, yp, zc)

    real(WP), dimension( d%ux_ysz(1), d%uy_ysz(2), d%ux_ysz(3) ) :: qx_yppc_ypencil ! <ux>^y at (xp, yp, zc)
    
    real(WP), dimension( d%uz_ysz(1), d%uy_ysz(2), d%uz_ysz(3) ) :: qz_ycpp_ypencil ! <uz>^y at (xc, yp, zp)
    
    real(WP), dimension( d%ux_ysz(1), d%ux_ysz(2), d%uz_ysz(3) ) :: qx_zpcp_ypencil ! <ux>^z at (xp, yc, zp), intermediate
    
    
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uz_ysz(3) ) :: qy_zcpp_ypencil ! <uy>^z at (xc, yp, zp)
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) ::    pres_ypencil ! p 
!-------------------------------------------------------------------------------
! ithermo == 0, z-pencil
!-------------------------------------------------------------------------------
    real(WP), dimension( d%ux_zsz(1), d%ux_zsz(2), d%ux_zsz(3) ) ::      qx_zpencil
    real(WP), dimension( d%uy_zsz(1), d%uy_zsz(2), d%uy_zsz(3) ) ::      qy_zpencil
    real(WP), dimension( d%uz_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) ::      qz_zpencil
    real(WP), dimension( d%ux_zsz(1), d%ux_zsz(2), d%uz_zsz(3) ) :: qx_zpcp_zpencil ! <ux>^z at (xp, yc, zp)
    real(WP), dimension( d%uy_zsz(1), d%uy_zsz(2), d%uz_zsz(3) ) :: qy_zcpp_zpencil ! <uy>^z at (xc, yp, zp)
    
    real(WP), dimension( d%uz_zsz(1), d%uy_zsz(2), d%uz_zsz(3) ) :: qz_ycpp_zpencil ! <uz>^y at (xc, yp, zp)
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) ::    pres_zpencil ! p 
!-------------------------------------------------------------------------------
! ithermo == 1, x-pencil
!-------------------------------------------------------------------------------
    real(WP), dimension( d%uz_xsz(1), d%ps_xsz(2), d%uz_xsz(3) ) ::   qz_zccc ! <uz>^z at (xc, yc, zc)
    real(WP), dimension( d%ps_xsz(1), d%ux_xsz(2), d%ux_xsz(3) ) ::   qx_xccc ! <ux>^x at (xc, yc, zc), intermediate
    real(WP), dimension( d%uy_xsz(1), d%ps_xsz(2), d%uy_xsz(3) ) ::   qy_yccc ! <uy>^y at (xc, yc, zc)
    real(WP), dimension( d%ps_xsz(1), d%ux_xsz(2), d%ux_xsz(3) ) ::   gx_xccc ! <gx>^x at (xc, yc, zc), g_i = rho * u_i
    real(WP), dimension( d%ux_xsz(1), d%uy_xsz(2), d%uy_xsz(3) ) ::   gy_xppc ! <gy>^x at (xp, yp, zc)
    real(WP), dimension( d%ux_xsz(1), d%uz_xsz(2), d%uz_xsz(3) ) ::   gz_xpcp ! <gz>^x at (xp, yc, zp)

    real(WP), dimension( d%ux_xsz(1), d%uy_xsz(2), d%ux_xsz(3) ) ::   gx_yppc ! <gx>^y at (xp, yp, zc)
    real(WP), dimension( d%ux_xsz(1), d%ux_xsz(2), d%uz_xsz(3) ) ::   gx_zpcp ! <gx>^z at (xp, yc, zp)
    
    real(WP), dimension( d%ux_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) ::    m_xpcc ! <mu>^x       at (xp, yc, zc)
    real(WP), dimension( d%ps_xsz(1), d%uy_xsz(2), d%ps_xsz(3) ) ::    m_ycpc ! <mu>^y       at (xc, yp, zc)
    real(WP), dimension( d%ps_xsz(1), d%ps_xsz(2), d%uz_xsz(3) ) ::    m_zccp ! <mu>^z       at (xc, yc, zp)
    real(WP), dimension( d%ux_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) :: dmdx_xpcc ! d( mu   )/dx at (xp, yc, zc)
    real(WP), dimension( d%ps_xsz(1), d%uy_xsz(2), d%ps_xsz(3) ) :: dmdx_ycpc ! d(<mu>^y)/dx at (xc, yp, zc)
    real(WP), dimension( d%ps_xsz(1), d%ps_xsz(2), d%uz_xsz(3) ) :: dmdx_zccp ! d(<mu>^z)/dx at (xc, yc, zp)
    real(WP), dimension( d%ux_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) :: dmdy_xpcc ! d(<mu>^x)/dy at (xp, yc, zc)
    real(WP), dimension( d%ux_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) :: dmdz_xpcc ! d(<mu>^x)/dz at (xp, yc, zc)

    real(WP), dimension( d%ps_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) :: div0, div ! divergence at cell centre at (xc, yc, zc)
!-------------------------------------------------------------------------------
! ithermo == 1, y-pencil
!-------------------------------------------------------------------------------
    real(WP), dimension( d%ps_ysz(1), d%ux_ysz(2), d%ux_ysz(3) ) :: qx_xccc_ypencil ! <ux>^x at (xc, yc, zc)
    real(WP), dimension( d%uy_ysz(1), d%ps_ysz(2), d%uy_ysz(3) ) :: qy_yccc_ypencil ! <uy>^y at (xc, yc, zc)
    real(WP), dimension( d%ux_ysz(1), d%uz_ysz(2), d%uz_ysz(3) ) :: qz_xpcp_ypencil ! <uz>^x at (xp, yc, zp), intermediate
    real(WP), dimension( d%uz_ysz(1), d%ps_ysz(2), d%uz_ysz(3) ) :: qz_zccc_ypencil ! <uz>^z at (xc, yc, zc), intermediate

    real(WP), dimension( d%ux_ysz(1), d%ux_ysz(2), d%ux_ysz(3) ) ::      gx_ypencil
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) ::      gy_ypencil
    real(WP), dimension( d%uz_ysz(1), d%uz_ysz(2), d%uz_ysz(3) ) ::      gz_ypencil
    real(WP), dimension( d%ux_ysz(1), d%uy_ysz(2), d%ux_ysz(3) ) :: gx_yppc_ypencil ! <gx>^y at (xp, yp, zc)
    real(WP), dimension( d%uy_ysz(1), d%ps_ysz(2), d%uy_ysz(3) ) :: gy_yccc_ypencil ! <gy>^y at (xc, yc, zc)
    real(WP), dimension( d%ux_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) :: gy_xppc_ypencil ! <gy>^x at (xp, yp, zc)
    real(WP), dimension( d%uz_ysz(1), d%uy_ysz(2), d%uz_ysz(3) ) :: gz_ycpp_ypencil ! <gz>^y at (xc, yp, zp)
    real(WP), dimension( d%ux_ysz(1), d%ux_ysz(2), d%uz_ysz(3) ) :: gx_zpcp_ypencil ! <gx>^z at (xp, yc, zp), intermediate
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uz_ysz(3) ) :: gy_zcpp_ypencil ! <gy>^z at (xc, yp, zp)
    
    real(WP), dimension( d%ux_ysz(1), d%uz_ysz(2), d%uz_ysz(3) ) :: gz_xpcp_ypencil ! <gz>^x at (xp, yc, zp), intermediate

    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) ::         m_ypencil
    real(WP), dimension( d%ux_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) ::    m_xpcc_ypencil ! <mu>^x       at (xp, yc, zc)
    real(WP), dimension( d%ps_ysz(1), d%uy_ysz(2), d%ps_ysz(3) ) ::    m_ycpc_ypencil ! <mu>^y       at (xc, yp, zc)
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%uz_ysz(3) ) ::    m_zccp_ypencil ! <mu>^z       at (xc, yc, zp)
    real(WP), dimension( d%ux_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) :: dmdy_xpcc_ypencil ! d(<mu>^x)/dy at (xp, yc, zc)
    real(WP), dimension( d%ps_ysz(1), d%uy_ysz(2), d%ps_ysz(3) ) :: dmdy_ycpc_ypencil ! d( mu   )/dy at (xc, yp, zc)
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%uz_ysz(3) ) :: dmdy_zccp_ypencil ! d(<mu>^z)/dy at (xc, yc, zp)
    real(WP), dimension( d%ux_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) :: dmdz_xpcc_ypencil ! d(<mu>^x)/dz at (xp, yc, zc), intermediate
    real(WP), dimension( d%ps_ysz(1), d%uy_ysz(2), d%ps_ysz(3) ) :: dmdx_ycpc_ypencil ! d(<mu>^y)/dx at (xc, yp, zc)
    real(WP), dimension( d%ps_ysz(1), d%uy_ysz(2), d%ps_ysz(3) ) :: dmdz_ycpc_ypencil ! d(<mu>^y)/dz at (xc, yp, zc)
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%uz_ysz(3) ) :: dmdx_zccp_ypencil ! d(<mu>^z)/dx at (xc, yc, zp)
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) ::       div_ypencil

    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3) ) ::    dDens_ypencil ! d 
!-------------------------------------------------------------------------------
! ithermo == 1, z-pencil
!-------------------------------------------------------------------------------
    real(WP), dimension( d%ps_zsz(1), d%ux_zsz(2), d%ux_zsz(3) ) ::   qx_xccc_zpencil ! <ux>^x at (xc, yc, zc)
    real(WP), dimension( d%ux_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) ::   qz_xpcp_zpencil ! <uz>^x at (xp, yc, zp)
    real(WP), dimension( d%uy_zsz(1), d%ps_zsz(2), d%uy_zsz(3) ) ::   qy_yccc_zpencil ! <uy>^y at (xc, yc, zc)
    real(WP), dimension( d%uz_zsz(1), d%ps_zsz(2), d%uz_zsz(3) ) ::   qz_zccc_zpencil ! <uz>^z at (xc, yc, zc)

    real(WP), dimension( d%ux_zsz(1), d%ux_zsz(2), d%ux_zsz(3) ) ::        gx_zpencil
    real(WP), dimension( d%uy_zsz(1), d%uy_zsz(2), d%uy_zsz(3) ) ::        gy_zpencil
    real(WP), dimension( d%uz_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) ::        gz_zpencil
    real(WP), dimension( d%ux_zsz(1), d%ux_zsz(2), d%uz_zsz(3) ) ::   gx_zpcp_zpencil ! <gx>^z at (xp, yc, zp)
    real(WP), dimension( d%uy_zsz(1), d%uy_zsz(2), d%uz_zsz(3) ) ::   gy_zcpp_zpencil ! <gy>^z at (xc, yp, zp)
    real(WP), dimension( d%uz_zsz(1), d%uz_zsz(2), d%ps_zsz(3) ) ::   gz_zccc_zpencil ! <gz>^z at (xc, yc, zc)
    real(WP), dimension( d%ux_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) ::   gz_xpcp_zpencil ! <gz>^x at (xp, yc, zp)
    real(WP), dimension( d%uz_zsz(1), d%uy_zsz(2), d%uz_zsz(3) ) ::   gz_ycpp_zpencil ! <gz>^y at (xc, yp, zp)
    
    

    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) ::         m_zpencil
    real(WP), dimension( d%ux_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) ::    m_xpcc_zpencil ! <mu>^x       at (xp, yc, zc)
    real(WP), dimension( d%ps_zsz(1), d%uy_zsz(2), d%ps_zsz(3) ) ::    m_ycpc_zpencil ! <mu>^y       at (xc, yp, zc)
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%uz_zsz(3) ) ::    m_zccp_zpencil ! <mu>^z       at (xc, yc, zp)
    real(WP), dimension( d%ux_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) :: dmdz_xpcc_zpencil ! d(<mu>^x)/dz at (xp, yc, zc)
    real(WP), dimension( d%ps_zsz(1), d%uy_zsz(2), d%ps_zsz(3) ) :: dmdz_ycpc_zpencil ! d(<mu>^y)/dz at (xc, yp, zc)
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%uz_zsz(3) ) :: dmdz_zccp_zpencil ! d( mu   )/dz at (xc, yc, zp)
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%uz_zsz(3) ) :: dmdx_zccp_zpencil ! d(<mu>^z)/dx at (xc, yc, zp)
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%uz_zsz(3) ) :: dmdy_zccp_zpencil ! d(<mu>^z)/dy at (xc, yc, zp)
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) ::       div_zpencil

    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3) ) ::    dDens_zpencil ! d 
!-------------------------------------------------------------------------------
! common vars
!-------------------------------------------------------------------------------
    ! natural position as in staggered storage
    real(WP), dimension( d%ux_xsz(1), d%ux_xsz(2), d%ux_xsz(3) ) :: mx_rhs ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( d%uy_xsz(1), d%uy_xsz(2), d%uy_xsz(3) ) :: my_rhs ! rhs for momentum-y at (xc, yp, zc)
    real(WP), dimension( d%uz_xsz(1), d%uz_xsz(2), d%uz_xsz(3) ) :: mz_rhs ! rhs for momentum-z at (xc, yc, zp)

    real(WP), dimension( d%ux_ysz(1), d%ux_ysz(2), d%ux_ysz(3) ) ::  mx_rhs_ypencil ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) ::  my_rhs_ypencil ! rhs for momentum-y at (xc, yp, zc)
    real(WP), dimension( d%uz_ysz(1), d%uz_ysz(2), d%uz_ysz(3) ) ::  mz_rhs_ypencil ! rhs for momentum-z at (xc, yc, zp)
    real(WP), dimension( d%ux_ysz(1), d%ux_ysz(2), d%ux_ysz(3) ) :: mx_temp_ypencil ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) :: my_temp_ypencil ! rhs for momentum-y at (xc, yp, zc)
    real(WP), dimension( d%uz_ysz(1), d%uz_ysz(2), d%uz_ysz(3) ) :: mz_temp_ypencil ! rhs for momentum-z at (xc, yc, zp)
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3) ) :: my_rhs_implicit_ypencil ! rhs for momentum-y at (xc, yp, zc)

    real(WP), dimension( d%ux_zsz(1), d%ux_zsz(2), d%ux_zsz(3) ) ::  mx_rhs_zpencil ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( d%uy_zsz(1), d%uy_zsz(2), d%uy_zsz(3) ) ::  my_rhs_zpencil ! rhs for momentum-y at (xc, yp, zc)
    real(WP), dimension( d%uz_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) ::  mz_rhs_zpencil ! rhs for momentum-z at (xc, yc, zp)
    real(WP), dimension( d%ux_zsz(1), d%ux_zsz(2), d%ux_zsz(3) ) :: mx_temp_zpencil ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( d%uy_zsz(1), d%uy_zsz(2), d%uy_zsz(3) ) :: my_temp_zpencil ! rhs for momentum-y at (xc, yp, zc)
    real(WP), dimension( d%uz_zsz(1), d%uz_zsz(2), d%uz_zsz(3) ) :: mz_temp_zpencil ! rhs for momentum-z at (xc, yc, zp)
    real(WP), dimension( d%uy_zsz(1), d%uy_zsz(2), d%uy_zsz(3) ) :: mz_rhs_implicit_zpencil ! rhs for momentum-y at (xc, yp, zc)

    real(WP), dimension( d%ux_xsz(1), d%ux_xsz(2), d%ux_xsz(3) ) :: mx_rhs_implicit ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( d%uy_xsz(1), d%uy_xsz(2), d%uy_xsz(3) ) :: my_rhs_implicit ! rhs for momentum-y at (xc, yp, zc)
    real(WP), dimension( d%uz_xsz(1), d%uz_xsz(2), d%uz_xsz(3) ) :: mz_rhs_implicit ! rhs for momentum-z at (xc, yc, zp)

    real(WP), dimension( d%ux_xsz(1), d%ux_xsz(2), d%ux_xsz(3) ) :: mx_rhs_implicit0 ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( d%uy_xsz(1), d%uy_xsz(2), d%uy_xsz(3) ) :: my_rhs_implicit0 ! rhs for momentum-y at (xc, yp, zc)
    real(WP), dimension( d%uz_xsz(1), d%uz_xsz(2), d%uz_xsz(3) ) :: mz_rhs_implicit0 ! rhs for momentum-z at (xc, yc, zp)

    integer(4), parameter :: II = 1, JJ = 2, KK = 3
    real(WP)              :: one_third_rre, two_third_rre, two_rre
!-------------------------------------------------------------------------------
!   Ensure it is in x-pencil
!-------------------------------------------------------------------------------
    if(d%ux_xsz(1) /= d%np(1)) call Print_error_msg("Error, not X-pencil")
!-------------------------------------------------------------------------------
! x-pencil : Initilisation
!-------------------------------------------------------------------------------
    one_third_rre = ONE_THIRD * f%rre ! = 1/(3Re)
    two_third_rre = TWO_THIRD * f%rre ! = 2/(3Re)
    two_rre       = TWO * f%rre       ! = 2/Re
!    
    f%mx_rhs(:, :, :) = ZERO
    f%my_rhs(:, :, :) = ZERO
    f%mz_rhs(:, :, :) = ZERO
!
    mx_rhs_implicit(:, :, :) = ZERO
    my_rhs_implicit(:, :, :) = ZERO
    mz_rhs_implicit(:, :, :) = ZERO
!
    if(ithermo == 1) then
      div(:, :, :) = ZERO
    end if
!_______________________________________________________________________________
! interpolation
!_______________________________________________________________________________
    call Get_x_midp_C2P_3dArray ( f%qy, d, qy_xppc) ! used in y-mom, w/o thermal 
    call Get_x_midp_C2P_3dArray ( f%qz, d, qz_xpcp) ! used in z-mom, w/o thermal
    call Get_y_midp_C2P_3dArray ( f%qx, d, qx_yppc) ! used in x-mom, w/o thermal

    call transpose_x_to_y(f%qy, qy_ypencil, d%dcpc)                 ! used in y-mom, w/o thermal
    call transpose_y_to_z(qy_ypencil, qy_zpencil, d%dcpc)           ! used in y-mom, w/o thermal
    call Get_z_midp_C2P_3dArray ( qy_zpencil, d, qy_zcpp_zpencil )  ! used in y-mon, w/o thermal
    
    call transpose_x_to_y(f%qz, qz_ypencil, d%dcpc)                   ! used in z-mom, w/o thermal
    call Get_y_midp_C2P_3dArray ( qz_ypencil, d, qz_ycpp_ypencil )    ! used in mom-z, w/o thermal

    call transpose_y_to_z(qz_ypencil, qz_zpencil, d%ccp)
    call transpose

  if(ithermo == 0) then
    call Get_y_midp_C2P_3dArray ( f%qx, d, qx_yppc) ! used in x-mom, wo thermal
    call Get_z_midp_C2P_3dArray ( f%qx, d, qx_zpcp) ! used in z-mom, wo thermal
    call transpose_z_to_y(qy_zcpp_zpencil, qy_zcpp_ypencil, d%dcpp) ! used in z-mom, wo  thermal
  end if

  if(ithermo == 1) then

  end if

  
  qz_xpcp_zpencil

!_______________________________________________________________________________
! the RHS of momentum equation
! the RHS terms of all 3 momentum equations (derivative) operating in the x direction
!_______________________________________________________________________________
!-------------------------------------------------------------------------------
!     x-pencil : for x-mom convection term (x-c1/3): d(gx * qx)/dx at (i', j, k)
!-------------------------------------------------------------------------------
  if(ithermo == 0) call Get_x_1st_derivative_P2P_3dArray( -f%qx * f%qx, d, mx_rhs )
  if(ithermo == 1) call Get_x_1st_derivative_P2P_3dArray( -f%gx * f%qx, d, mx_rhs )
  f%mx_rhs = f%mx_rhs + mx_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for x-mom diffusion term (x-v1/1), \mu * L11(ux) at (i', j, k)
!-------------------------------------------------------------------------------
  call Get_x_2nd_derivative_P2P_3dArray( f%qx, d, mx_rhs )
  f%mx_rhs = f%mx_rhs + f%rre * mx_rhs




!_______________________________________________________________________________
! the RHS of momentum equation
! the RHS terms of all 3 momentum equations (derivative) operating in the x direction
!_______________________________________________________________________________
    if(ithermo == 0 ) then
!-------------------------------------------------------------------------------
!     x-pencil : for x-mom convection term (x-c1/3): d(qx * qx)/dx at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_P2P_3dArray( -f%qx * f%qx, d, mx_rhs )
      f%mx_rhs = f%mx_rhs + mx_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for x-mom diffusion term (x-v1/1), \mu * L11(ux) at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_x_2nd_derivative_P2P_3dArray( f%qx, d, mx_rhs )
      f%mx_rhs = f%mx_rhs + f%rre * mx_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for y-mom convection term (y-c1/3), d(qx^y * qy^x)/dx at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_P2C_3dArray( -
       * qy_xppc, d, my_rhs )
      f%my_rhs = f%my_rhs + my_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for z-mom convection term (z-c1/3), d(qx^z * qz^x)/dx at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_P2C_3dArray( -qx_zpcp * qz_xpcp, d, mz_rhs )
      f%mz_rhs = f%mz_rhs + mz_rhs


    else if(ithermo == 1) then
!-------------------------------------------------------------------------------
!     x-pencil : for x-mom convection term (x-c1/3): d(gx * qx)/dx at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_P2P_3dArray( -f%gx * f%qx, d, mx_rhs )
      f%mx_rhs = f%mx_rhs + mx_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for y-mom convection term (y-c1/3), d(gx^y * qy^x)/dx at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_P2C_3dArray( -gx_yppc * qy_xppc, d, my_rhs )
      f%my_rhs = f%my_rhs + my_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for z-mom convection term (z-c1/3), d(gx^z * qz^x)/dx at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_P2C_3dArray( -gx_zpcp * qz_xpcp, d, mz_rhs )
      f%mz_rhs = f%mz_rhs + mz_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for x-mom diffusion term (x-v1/7), \mu^x * Ljj(ux) at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_x_2nd_derivative_P2P_3dArray( f%qx, d, mx_rhs )
      f%mx_rhs = f%mx_rhs + f%rre * m_xpcc * mx_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for x-mom diffusion term (x-v2/7), \mu^x * 1/3 * d (div)/dx at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3dArray( div, d, mx_rhs )
      f%mx_rhs = f%mx_rhs + one_third_rre * m_xpcc * mx_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for x-mom diffusion term (x-v3/7), -2/3 * d\mu/dx * (div(u)^x) +  
!                                                      2 * d\mu/dx * du/dx
!-------------------------------------------------------------------------------
      call Get_x_midp_C2P_3dArray          ( div, d, mx_rhs )
      f%mx_rhs = f%mx_rhs - two_third_rre * dmdx_xpcc * mx_rhs
      call Get_x_1st_derivative_P2P_3dArray( f%qx, d, mx_rhs )
      f%mx_rhs = f%mx_rhs + two_rre * dmdx_xpcc * mx_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for x-mom diffusion term (x-v4/7), d(mu^x)/dy * d(qy^y)/dx at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3dArray( qy_yccc, d, mx_rhs )
      f%mx_rhs =  f%mx_rhs + f%rre * dmdy_xpcc * mx_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for x-mom diffusion term (x-v6/7), d(mu^x)/dz * d(qz^z)/dx at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3dArray( qz_zccc, d, mx_rhs )
      f%mx_rhs =  f%mx_rhs + f%rre * dmdz_xpcc * mx_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for y-mom diffusion term (y-v5/7), d(mu^y)/dx * d(qy))/dx at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2C_3dArray( f%qy, d, my_rhs )
      f%my_rhs =  f%my_rhs + f%rre * dmdx_ycpc * my_rhs
!-------------------------------------------------------------------------------
!     x-pencil : for z-mom diffusion term (z-v5/7), d(mu^z)/dx * d(qz)/dx at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2C_3dArray( f%qz, d, mz_rhs )
      f%mz_rhs =  f%mz_rhs + f%rre * dmdx_zccp * mz_rhs
    else
    end if
!-------------------------------------------------------------------------------
!     x-pencil : for x-mom, pressure gradient in x direction, sigma_1*d(p)/dx
!-------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3dArray( f%pres, d, mx_rhs_implicit0 )
    mx_rhs_implicit =  mx_rhs_implicit - sigma1p * mx_rhs_implicit0
!-------------------------------------------------------------------------------
!     x-pencil : gravity force in x direction
!-------------------------------------------------------------------------------
    if( ithermo == 1 .and. ( igravity == 1 .or. igravity == -1) ) then
      call Get_x_midp_C2P_3dArray( f%dDens, d, mx_rhs_implicit0 )
      mx_rhs_implicit =  mx_rhs_implicit + f%fgravity * mx_rhs_implicit0
    end if
!_______________________________________________________________________________
! the RHS of momentum equation
! the RHS terms of all 3 momentum equations (derivative) operating in the y direction
!_______________________________________________________________________________
    mx_rhs_ypencil = ZERO
    my_rhs_ypencil = ZERO
    mz_rhs_ypencil = ZERO

    mx_temp_ypencil = ZERO
    my_temp_ypencil = ZERO
    mz_temp_ypencil = ZERO

    if(ithermo == 0 ) then
!-------------------------------------------------------------------------------
! y-pencil : for x-mom convection term (x-c2/3): d(qy^x * qx^y)/dy at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_P2C_3dArray( -qy_xppc_ypencil * qx_yppc_ypencil,  d, mx_temp_ypencil )
      mx_rhs_ypencil = mx_rhs_ypencil + mx_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for x-mom diffusion term (x-v1/2), \mu * L22(ux) at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_y_2nd_derivative_P2P_3dArray( qx_ypencil,  d, mx_temp_ypencil )
      mx_rhs_ypencil = mx_rhs_ypencil + f%rre * mx_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for y-mom convection term (y-c2/3), d(qy * qy)/dy at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_P2P_3dArray( -qy_ypencil * qy_ypencil,  d, my_temp_ypencil )
      my_rhs_ypencil = my_rhs_ypencil + my_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for z-mom convection term (z-c2/3), d(qy * qz)/dy at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_P2C_3dArray( -qy_zcpp_ypencil * qz_ycpp_ypencil,  d, mz_temp_ypencil )
      mz_rhs_ypencil = mz_rhs_ypencil + mz_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for y-mom diffusion term (y-v1/1), \mu * Ljj(uy) at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_y_2nd_derivative_P2P_3dArray( qy_ypencil,  d, my_temp_ypencil )
      my_rhs_ypencil = my_rhs_ypencil + f%rre * my_temp_ypencil

    else if (ithermo == 1) then
!-------------------------------------------------------------------------------
! y-pencil : for x-mom convection term (x-c2/3): d(gy * qx)/dy at (i', j, k)
!-------------------------------------------------------------------------------
      call transpose_x_to_y(qx_yppc, qx_yppc_ypencil, d%dppc)
      call Get_y_1st_derivative_P2C_3dArray( -gy_xppc_ypencil * qx_yppc_ypencil,  d, mx_temp_ypencil )
      mx_rhs_ypencil = mx_rhs_ypencil + mx_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for y-mom convection term (y-c2/3), d(gy * qy)/dy at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_P2P_3dArray( -gy_ypencil * qy_ypencil,  d, my_temp_ypencil )
      my_rhs_ypencil = my_rhs_ypencil + my_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for z-mom convection term (z-c2/3), d(gy * qz)/dy at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_P2C_3dArray( -gy_zcpp_ypencil * qz_ycpp_ypencil,  d, mz_temp_ypencil )
      mz_rhs_ypencil = mz_rhs_ypencil + mz_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for y-mom diffusion term (y-v1/7), \mu * Ljj(uy) at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_y_2nd_derivative_P2P_3dArray( qy_ypencil,  d, my_temp_ypencil )
      my_rhs_ypencil = my_rhs_ypencil + f%rre * m_ycpc_ypencil * my_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for y-mom diffusion term (y-v2/7), \mu * 1/3 * d (div)/dy at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2P_3dArray( div_ypencil, d, my_temp_ypencil )
      my_rhs_ypencil = my_rhs_ypencil + one_third_rre * m_ycpc_ypencil * my_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for y-mom diffusion term (y-v3/7), 2d\mu/dy * (-1/3 * div(u)) +  2d\mu/dy * dv/dy
!-------------------------------------------------------------------------------
      call Get_y_midp_C2P_3dArray          ( div_ypencil, d, my_temp_ypencil )
      my_rhs_ypencil = my_rhs_ypencil - two_third_rre * dmdy_ycpc_ypencil * my_temp_ypencil
      call Get_y_1st_derivative_P2P_3dArray( qy_ypencil, d, my_temp_ypencil )
      my_rhs_ypencil = my_rhs_ypencil + two_rre       * dmdy_ycpc_ypencil * my_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for y-mom diffusion term (y-v4/7), d(mu)/dx * d(qx)/dy at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2P_3dArray( qx_xccc_ypencil, d, my_temp_ypencil )
      my_rhs_ypencil =  my_rhs_ypencil + f%rre * dmdx_ycpc_ypencil * my_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for y-mom diffusion term (y-v6/7), d(mu)/dz * d(qz)/dy at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2P_3dArray( qz_zccc_ypencil, d, my_temp_ypencil )
      my_rhs_ypencil =  my_rhs_ypencil + f%rre * dmdz_ycpc_ypencil * my_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for x-mom diffusion term (x-v5/7), d(mu)/dy * d(qx)/dy at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2C_3dArray( qx_ypencil, d, mx_temp_ypencil )
      mx_rhs_ypencil =  mx_rhs_ypencil + f%rre * dmdy_xpcc_ypencil * mx_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for z-mom diffusion term (z-v7/7), d(mu)/dy * d(qz)/dy at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2C_3dArray( qz_ypencil, d, mz_temp_ypencil )
      mz_rhs_ypencil =  mz_rhs_ypencil + f%rre * dmdy_zccp_ypencil * mz_temp_ypencil

    else
    end if
!-------------------------------------------------------------------------------
! y-pencil : for y-mom pressure gradient in y direction, d(sigma_1 p)
!-------------------------------------------------------------------------------
    call transpose_x_to_y(f%pres, pres_ypencil, d%dccc)
    call Get_y_1st_derivative_C2P_3dArray( pres_ypencil, d, my_temp_ypencil )
    my_rhs_implicit_ypencil =  my_rhs_implicit_ypencil - sigma1p * my_temp_ypencil
!-------------------------------------------------------------------------------
! y-pencil : for y-mom gravity force in y direction
!-------------------------------------------------------------------------------
    if( ithermo == 1 .and. ( igravity == 2 .or. igravity == -2) ) then
      call tranpose_x_to_y(f%dDens, dDens_ypencil, d%dccc)
      call Get_y_midp_C2P_3dArray( dDens_ypencil, d, my_temp_ypencil )
      my_rhs_implicit_ypencil =  my_rhs_implicit_ypencil + f%fgravity * my_temp_ypencil
    end if
!-------------------------------------------------------------------------------
! Data from Y-pencil to Z-pencil
!-------------------------------------------------------------------------------    
    call transpose_from_y_to_x(mx_rhs_ypencil, mx_rhs, d%dpcc)
    call transpose_from_y_to_x(my_rhs_ypencil, my_rhs, d%dcpc)
    call transpose_from_y_to_x(mz_rhs_ypencil, mz_rhs, d%dccp)
    f%mx_rhs = f%mx_rhs + mx_rhs
    f%my_rhs = f%my_rhs + my_rhs
    f%mz_rhs = f%mz_rhs + mz_rhs
    call transpose_from_y_to_x(my_rhs_implicit_ypencil, my_rhs_implicit0, d%dcpc)
    my_rhs_implicit = my_rhs_implicit + my_rhs_implicit0
!_______________________________________________________________________________
! the RHS of momentum equation
! the RHS terms of all 3 momentum equations (derivative) operating in the z direction
!_______________________________________________________________________________
     mx_rhs_zpencil = ZERO
     my_rhs_zpencil = ZERO
     mz_rhs_zpencil = ZERO

    mx_temp_zpencil = ZERO
    my_temp_zpencil = ZERO
    mz_temp_zpencil = ZERO
    if(ithermo == 0 ) then
!-------------------------------------------------------------------------------
! z-pencil : for x-mom convection term (x-c3/3): d(qz * qx)/dz at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_P2C_3dArray( -qz_xpcp_zpencil * qx_zpcp_zpencil,  d, mx_temp_zpencil )
      mx_rhs_zpencil = mx_rhs_zpencil + mx_temp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : for y-mom convection term (y-c3/3), d(qz * qy)/dz at (i, j', k)
!-------------------------------------------------------------------------------
      call transpose_y_to_z(qz_ycpp_ypencil, qz_ycpp_zpencil, d%dcpp)
      call Get_z_1st_derivative_P2C_3dArray( -qz_ycpp_zpencil * qy_zcpp_zpencil,  d, my_temp_zpencil )
      my_rhs_zpencil = my_rhs_zpencil + my_temp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : for z-mom convection term (z-c3/3), d(qz * qz)/dz at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_P2P_3dArray( -qz_zpencil * qz_zpencil,  d, mz_temp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil + mz_temp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : for z-mom diffusion term (z-v1/1), \mu * Ljj(uz) at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_z_2nd_derivative_P2P_3dArray( qz_zpencil,  d, mz_temp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil + f%rre * mz_temp_zpencil

    else if (ithermo == 1) then
!-------------------------------------------------------------------------------
! z-pencil : for x-mom convection term (x-c3/3): d(gz * qx)/dz at (i', j, k)
!-------------------------------------------------------------------------------
      call transpose_x_to_y(gz_xpcp,         gz_xpcp_ypencil, d%dpcp)
      call transpose_y_to_z(gz_xpcp_ypencil, gz_xpcp_zpencil, d%dpcp)
      call Get_z_1st_derivative_P2C_3dArray( -gz_xpcp_zpencil * qx_zpcp_zpencil,  d, mx_temp_zpencil )
      mx_rhs_zpencil = mx_rhs_zpencil + mx_temp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : for y-mom convection term (y-c3/3), d(gz * qy)/dz at (i, j', k)
!-------------------------------------------------------------------------------
      call transpose_y_to_z(gz_ycpp_ypencil, gz_ycpp_zpencil, d%dcpp)
      call Get_z_1st_derivative_P2C_3dArray( -gz_ycpp_zpencil * qy_zcpp_zpencil,  d, my_temp_zpencil )
      my_rhs_zpencil = my_rhs_zpencil + my_temp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : for z-mom convection term (z-c3/3), d(gz * qz)/dz at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_P2P_3dArray( -gz_zpencil * qz_zpencil,  d, mz_temp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil + mz_temp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : for z-mom diffusion term (z-v1/7), \mu * Ljj(uz) at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_z_2nd_derivative_P2P_3dArray( qz_zpencil,  d, mz_temp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil + f%rre * m_zccp_zpencil * mz_temp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : for z-mom diffusion term (z-v2/7), \mu * 1/3 * d (div)/dz at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2P_3dArray( div_zpencil, d, mz_temp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil + one_third_rre * m_zccp_zpencil * mz_temp_zpencil
!-------------------------------------------------------------------------------
      ! z-pencil : for z-mom diffusion term (z-v3/7), 2d\mu/dz * (-1/3 * div(u)) +  2d\mu/dz * dw/dz
!-------------------------------------------------------------------------------
      call Get_z_midp_C2P_3dArray          ( div_zpencil, d, mz_temp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil - two_third_rre * dmdz_zccp_zpencil * mz_temp_zpencil
      call Get_z_1st_derivative_P2P_3dArray( f%qz_zpencil, d, mz_temp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil + two_rre * dmdz_zccp_zpencil * mz_temp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : for z-mom diffusion term (z-v4/7), d(mu)/dx * d(qx)/dz at (i, j, k')
!-------------------------------------------------------------------------------
      call transpose_y_to_z(qx_xccc_ypencil, qx_xccc_zpencil, d%dccc)
      call transpose_x_to_y(dmdx_zccp,         dmdx_zccp_ypencil, d%dccp)
      call transpose_y_to_z(dmdx_zccp_ypencil, dmdx_zccp_zpencil, d%dccp)
      call Get_z_1st_derivative_C2P_3dArray( qx_xccc_zpencil, d, mz_temp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil + f%rre * dmdx_zccp_zpencil * mz_temp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : for z-mom diffusion term (z-v6/7), d(mu)/dy * d(qy)/dz at (i, j, k')
!-------------------------------------------------------------------------------
      call transpose_y_to_z(qy_yccc_ypencil, qy_yccc_zpencil, d%dccc)
      call transpose_y_to_z(dmdy_zccp_ypencil, dmdy_zccp_zpencil, d%dccc)
      call Get_z_1st_derivative_C2P_3dArray( qy_yccc_zpencil, d, mz_temp_zpencil )
      mz_rhs_zpencil = mz_rhs_zpencil + f%rre * dmdy_zccp_zpencil * mz_temp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : for x-mom diffusion term (x-v7/7), d(mu)/dz * d(qx)/dz at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2C_3dArray( qx_zpencil, d, mx_temp_zpencil )
      mx_rhs_zpencil = mx_rhs_zpencil + f%rre * dmdz_xpcc_zpencil * mx_temp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : for y-mom diffusion term (y-v7/7), d(mu)/dz * d(qy)/dz at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2C_3dArray( qy_zpencil, d, my_temp_zpencil )
      my_rhs_zpencil =  my_rhs_zpencil + f%rre * dmdz_ycpc_zpencil * my_temp_zpencil

    else
    end if
!-------------------------------------------------------------------------------
! z-pencil : pressure gradient in z direction, d(sigma_1 p)
!-------------------------------------------------------------------------------
    call transpose_y_to_z(pres_ypencil, pres_zpencil, d%dccc)
    call Get_z_1st_derivative_C2P_3dArray( pres_zpencil, d, mz_temp_zpencil )
    mz_rhs_implicit_zpencil =  mz_rhs_implicit_zpencil - sigma1p * mz_temp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : gravity force in z direction
!-------------------------------------------------------------------------------
    if( ithermo == 1 .and. ( igravity == 3 .or. igravity == -3) ) then
      call tranpose_y_to_z(dDens_ypencil, dDens_zpencil, d%dccc)
      call Get_z_midp_C2P_3dArray( dDens_zpencil, d, mz_temp_zpencil )
      mz_rhs_implicit_zpencil =  mz_rhs_implicit_zpencil + f%fgravity * mz_temp_zpencil
    end if
!-------------------------------------------------------------------------------
! Data from Z-pencil to X-pencil
!-------------------------------------------------------------------------------    
    call transpose_from_z_to_y(mx_rhs_zpencil, mx_rhs_ypencil, d%dpcc)
    call transpose_from_z_to_y(my_rhs_zpencil, my_rhs_ypencil, d%dcpc)
    call transpose_from_z_to_y(mz_rhs_zpencil, mz_rhs_ypencil, d%dccp)
    
    call transpose_from_y_to_x(mx_rhs_ypencil, mx_rhs, d%dpcc)
    call transpose_from_y_to_x(my_rhs_ypencil, my_rhs, d%dcpc)
    call transpose_from_y_to_x(mz_rhs_ypencil, mz_rhs, d%dccp)

    f%mx_rhs = f%mx_rhs + mx_rhs
    f%my_rhs = f%my_rhs + my_rhs
    f%mz_rhs = f%mz_rhs + mz_rhs

    call transpose_from_z_to_y(my_rhs_implicit_zpencil, my_rhs_implicit_ypencil, d%dccp)
    call transpose_from_y_to_x(my_rhs_implicit_ypencil, my_rhs_implicit0, d%dccp)

    my_rhs_implicit = my_rhs_implicit + my_rhs_implicit0
!-------------------------------------------------------------------------------
! x-pencil : to build up rhs in total, in all directions
!-------------------------------------------------------------------------------
    ! x-momentum
    call Calculate_momentum_fractional_step(f%mx_rhs0, f%mx_rhs, mx_rhs_implicit, isub)
    if(idriven /= IDRVF_NO) call Calculate_momentum_driven_source(f%mx_rhs, 'ux', d, isub) 

    ! y-momentum
    call Calculate_momentum_fractional_step(f%my_rhs0, f%my_rhs, my_rhs_implicit, isub)

    ! z-momentum
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
!  mode           name          role                                           !
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
!  mode           name          role                                           !
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
