module eq_momentum_mod
  use precision_mod, only : WP
  implicit none
  private :: Calculate_xmomentum_driven_source
  private :: Calculate_momentum_fractional_step
  private :: Compute_momentum_rhs
  private :: Correct_massflux
  public  :: Solve_momentum_eq

contains
!===============================================================================
!> \brief To calcuate the driven force for streamwise peridic flow.
!>
!> This subroutine is called everytime when calcuting the rhs of momentum eqs.
!> check
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                         
!______________________________________________________________________________!
!> \param[inout]  rhs          the rhs of the momentum equation
!> \param[in]     d            domain name
!> \param[in]     isub         the RK iteration to get correct Coefficient 
!===============================================================================
  subroutine Calculate_xmomentum_driven_source(isub, idriven, drvf, dm, rhs)
    use parameters_constant_mod
    implicit none

    type(t_domain), intent(in) :: dm
    real(WP),    intent(inout) :: rhs(:, :, :)
    integer,     intent(in) :: isub
    integer,     intent(in) :: idriven
    real(WP),       intent(in) :: drvf
    
    real(WP) :: rhs_bulk

    rhs_bulk = ZERO

    if(idriven == IDRVF_MASSFLUX) then

      call Get_volumetric_average_3d(.false., dm%ibcy(1,:), dm%fbcy(1,:), &
          dm, dm%dpcc, rhs, rhs_bulk)
      
    else if (idriven == IDRVF_SKINFRIC) then

      rhs_bulk = - HALF * drvf * dm%tAlpha(isub) * dt

    else if (idriven == IDRVF_PRESLOSS ) then

    ! to check this part
      rhs_bulk = - HALF * drvf * dm%tAlpha(isub) * dt

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
  subroutine Calculate_momentum_fractional_step(rhs0, rhs1, rhs1_semi, dm, isub)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    real(WP), dimension(:, :, :), intent(in   ) :: rhs1_semi
    real(WP), dimension(:, :, :), intent(inout) :: rhs0, rhs1
    integer,                   intent(in   ) :: isub
    type(t_domain), intent(in) :: dm
    

    integer :: n(3)
    real(WP), allocatable :: rhs_dummy(:, :, :)

    n(1:3) = shape(rhs1)
    allocate( rhs_dummy (n(1), n(2), n(3)) )

  ! add explicit terms
    rhs_dummy(:, :, :) = rhs1(:, :, :)
    rhs1(:, :, :) = dm%tGamma(isub) * rhs1(:, :, :) + &
                    dm%tZeta (isub) * rhs0(:, :, :)
    rhs0(:, :, :) = rhs_dummy(:, :, :)

  ! add implicit
    rhs1(:, :, :) = rhs1(:, :, :) + &
                    dm%tAlpha(isub) * rhs1_semi(:, :, :)
  
  ! times the time step 
    rhs1(:, :, :) = dm%dt * rhs1(:, :, :)

    deallocate (rhs_dummy)

    return
  end subroutine
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
  subroutine Compute_momentum_rhs(fl, dm, isub)
    use parameters_constant_mod
    use operations
    implicit none

    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(in   ) :: dm
    integer,     intent(in   ) :: isub
!-------------------------------------------------------------------------------
! common vars
!-------------------------------------------------------------------------------
    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: qx_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: qy_ypencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: qz_ypencil
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: qx_zpencil
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: qy_zpencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: qz_zpencil

    real(WP), dimension( dm%dpcc%ysz(1), dm%dcpc%ysz(2), dm%dpcc%ysz(3) ) :: qx_yppc_ypencil ! <ux>^y at (xp, yp, zc)
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dccp%zsz(3) ) :: qx_zpcp_zpencil ! <ux>^z at (xp, yc, zp)
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dccp%zsz(3) ) :: qy_zcpp_zpencil ! <uy>^z at (xc, yp, zp)
    real(WP), dimension( dm%dpcc%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: qz_xpcp         ! <uz>^x at (xp, yc, zp) 
    real(WP), dimension( dm%dccp%ysz(1), dm%dcpc%ysz(2), dm%dccp%ysz(3) ) :: qz_ycpp_ypencil ! <uz>^y at (xc, yp, zp)
    
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: pres_ypencil ! p
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: pres_zpencil ! p

    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: mx_rhs_ypencil ! 
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: my_rhs_ypencil ! 
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: mz_rhs_ypencil ! 
 
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: mx_rhs_zpencil ! 
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: my_rhs_zpencil ! 
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: mz_rhs_zpencil ! 
 
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: mx_rhs_implicit ! 
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: my_rhs_implicit_ypencil ! 
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: mz_rhs_implicit_zpencil ! 

    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: apcc ! 
    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: apcc_ypencil ! 
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: apcc_zpencil ! 

    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: acpc ! 
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil !
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: acpc_zpencil !  

    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: accp ! 
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil !
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil ! 

    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dccp%ysz(3) ) :: apcp_ypencil ! 
    
!-------------------------------------------------------------------------------
! thermal == 0 only
!-------------------------------------------------------------------------------
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dccp%xsz(3) ) :: qx_zpcp ! <ux>^z at (xp, yc, zp)
    real(WP), dimension( dm%dpcc%xsz(1), dm%dcpc%xsz(2), dm%dpcc%xsz(3) ) :: qx_yppc ! <ux>^y at (xp, yp, zc)
    real(WP), dimension( dm%dpcc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: qy_xppc ! <uy>^x at (xp, yp, zc)
    real(WP), dimension( dm%dpcc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: qy_xppc_ypencil ! <uy>^x at (xp, yp, zc)
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dccp%ysz(3) ) :: qy_zcpp_ypencil ! <uy>^z at (xc, yp, zp)
    real(WP), dimension( dm%dpcc%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: qz_xpcp_zpencil ! <uz>^x at (xp, yc, zp)
    real(WP), dimension( dm%dccp%zsz(1), dm%dcpc%zsz(2), dm%dccp%zsz(3) ) :: qz_ycpp_zpencil ! <uz>^y at (xc, yp, zp)
!-------------------------------------------------------------------------------
! thermal == 1 only
!-------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dccp%xsz(3) ) :: apcp
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dccp%zsz(3) ) :: apcp_zpencil
    real(WP), dimension( dm%dpcc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: appc_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dccp%ysz(3) ) :: acpp_ypencil
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dccp%zsz(3) ) :: acpp_zpencil
    
    real(WP), dimension( dm%dccc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: qx_xccc_ypencil   ! <ux>^x at (xc, yc, zc)
    real(WP), dimension( dm%dccc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: qx_xccc_zpencil   ! <ux>^x at (xc, yc, zc)
    real(WP), dimension( dm%dcpc%xsz(1), dm%dccc%xsz(2), dm%dcpc%xsz(3) ) :: qy_yccc           ! <uy>^y at (xc, yc, zc)
    real(WP), dimension( dm%dcpc%zsz(1), dm%dccc%zsz(2), dm%dcpc%zsz(3) ) :: qy_yccc_zpencil   ! <uy>^y at (xc, yc, zc)
    real(WP), dimension( dm%dccp%xsz(1), dm%dccc%xsz(2), dm%dccp%xsz(3) ) :: qz_zccc           ! <uz>^z at (xc, yc, zc)
    real(WP), dimension( dm%dccp%ysz(1), dm%dccc%ysz(2), dm%dccp%ysz(3) ) :: qz_zccc_ypencil   ! <uz>^z at (xc, yc, zc), intermediate

    real(WP), dimension( dm%dpcc%xsz(1), dm%dcpc%xsz(2), dm%dpcc%xsz(3) ) :: gx_yppc           ! <gx>^y at (xp, yp, zc)
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dccp%xsz(3) ) :: gx_zpcp           ! <gx>^z at (xp, yc, zp)
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: gy_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dccp%ysz(3) ) :: gy_zcpp_ypencil   ! <gy>^z at (xc, yp, zp)
    real(WP), dimension( dm%dpcc%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: gz_xpcp_zpencil   ! <gz>^x at (xp, yc, zp)
    real(WP), dimension( dm%dccp%zsz(1), dm%dcpc%zsz(2), dm%dccp%zsz(3) ) :: gz_ycpp_zpencil   ! <gz>^y at (xc, yp, zp)
    
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: dDens_ypencil  ! d 
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: dDens_zpencil  ! d 

    real(WP), dimension( dm%dpcc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: m_xpcc         ! <mu>^x       at (xp, yc, zc)  
    real(WP), dimension( dm%dpcc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: m_xpcc_ypencil ! <mu>^x       at (xp, yc, zc)
    real(WP), dimension( dm%dpcc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: m_xpcc_zpencil ! <mu>^x       at (xp, yc, zc)
    real(WP), dimension( dm%dccc%xsz(1), dm%dcpc%xsz(2), dm%dccc%xsz(3) ) :: m_ycpc         ! <mu>^y       at (xc, yp, zc)
    real(WP), dimension( dm%dccc%ysz(1), dm%dcpc%ysz(2), dm%dccc%ysz(3) ) :: m_ycpc_ypencil ! <mu>^y       at (xc, yp, zc)
    real(WP), dimension( dm%dccc%zsz(1), dm%dcpc%zsz(2), dm%dccc%zsz(3) ) :: m_ycpc_zpencil ! <mu>^y       at (xc, yp, zc)
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccp%xsz(3) ) :: m_zccp         ! <mu>^z       at (xc, yc, zp)
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccp%ysz(3) ) :: m_zccp_ypencil ! <mu>^z       at (xc, yc, zp)
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccp%zsz(3) ) :: m_zccp_zpencil ! <mu>^z       at (xc, yc, zp)
    
    real(WP), dimension( dm%dpcc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: dmdx_xpcc ! d( mu   )/dx at (xp, yc, zc)
    real(WP), dimension( dm%dccc%xsz(1), dm%dcpc%xsz(2), dm%dccc%xsz(3) ) :: dmdx_ycpc ! d(<mu>^y)/dx at (xc, yp, zc)
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccp%xsz(3) ) :: dmdx_zccp ! d(<mu>^z)/dx at (xc, yc, zp)
    real(WP), dimension( dm%dpcc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: dmdy_xpcc ! d(<mu>^x)/dy at (xp, yc, zc)
    real(WP), dimension( dm%dpcc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: dmdz_xpcc ! d(<mu>^x)/dz at (xp, yc, zc)

    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccp%zsz(3) ) :: dmdx_zccp_zpencil ! d(<mu>^z)/dx at (xc, yc, zp)
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccp%ysz(3) ) :: dmdy_zccp_ypencil ! d(<mu>^z)/dy at (xc, yc, zp)
    real(WP), dimension( dm%dpcc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: dmdy_xpcc_ypencil ! d(<mu>^x)/dy at (xp, yc, zc)
    real(WP), dimension( dm%dpcc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: dmdz_xpcc_zpencil ! d(<mu>^x)/dz at (xp, yc, zc)
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccp%zsz(3) ) :: dmdz_zccp_zpencil ! d( mu   )/dz at (xc, yc, zp)

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: div
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: div_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: div_zpencil
!-------------------------------------------------------------------------------
! others
!-------------------------------------------------------------------------------
    real(WP) :: one_third_rre, two_third_rre, two_rre
    real(WP) :: fbc(2)
    integer  :: ibc(2)
    integer  :: i
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
    call transpose_x_to_y (fl%pres,      pres_ypencil, dm%dccc)      ! y-pencil : y-mom, w+o   thermal
    call transpose_y_to_z (pres_ypencil, pres_zpencil, dm%dccc)      ! z-pencil : z-mom, w+o   thermal
!-------------------------------------------------------------------------------
!    qx --> qx_xccc (I) --> qx_xccc_ypencil(W) --> qx_xccc_zpencil(W)
!     | --> qx_ypencil(WO) --> qx_yppc_ypencil(WO) --> qx_yppc(O)
!                       |  --> qx_zpencil(WO) --> qx_zpcp_zpencil(WO) --> qx_zpcp_ypencil(I) --> qx_zpcp(O)
!-------------------------------------------------------------------------------
    i = 1 
    call transpose_x_to_y (fl%qx, qx_ypencil, dm%dpcc)  
    call Get_y_midp_C2P_3D(qx_ypencil, qx_yppc_ypencil, dm, dm%ibcy(i,:), dm%fbcy(i, :)) ! y-pencil : x-mom, w+o thermal
    if(dm%ithermo == 0) &
    call transpose_y_to_x (qx_yppc_ypencil, qx_yppc, dm%dppc)                    ! x-pencil : y-mom, o   thermal
  
    call transpose_y_to_z (qx_ypencil, qx_zpencil, dm%dpcc)                    ! z-pencil : x-mom, w+o thermal
    call Get_z_midp_C2P_3D(qx_zpencil, qx_zpcp_zpencil, dm, dm%ibcz(i,:), dm%fbcz(i, :)) ! z-pencil : x-mom, w+o thermal
    if(dm%ithermo == 0) then
    call transpose_z_to_y (qx_zpcp_zpencil, apcp_ypencil,    dm%dpcp)                    ! intermediate, apcp_ypencil = qx_zpcp_ypencil
    call transpose_y_to_x (apcp_ypencil,    qx_zpcp,         dm%dpcp)                    ! x-pencil : z-mom,  o  thermal
    end if
    if(dm%ithermo == 1) then
    call Get_x_midp_P2C_3D( fl%qx,          accc,       dm, dm%ibcx(i,:) ) ! intermediate, accc = qx_xccc
    call transpose_x_to_y (accc,            qx_xccc_ypencil, dm%dccc)                    ! y-pencil : y-mom, w   thermal
    call transpose_y_to_z (qx_xccc_ypencil, qx_xccc_zpencil, dm%dccc)                    ! z-pencil : z-mom, w   thermal
    end if
!-------------------------------------------------------------------------------
!    qy--> qy_xppc(WO) --> qy_xppc_ypencil(O) 
!     |--> qy_ypencil(WO) --> qy_zpencil(WO) --> qy_zcpp_zpencil(WO) --> qy_zcpp_ypencil(O)
!                       | --> qy_yccc_ypencil(I) --> qy_yccc(W)
!                                              | --> qy_yccc_zpencil(W)              
!-------------------------------------------------------------------------------
    i = 2 
    call Get_x_midp_C2P_3D(fl%qy, qy_xppc, dm, dm%ibcx(i,:), dm%fbcx(i, :)) ! xpencil : y-mom, w+o thermal
    if(dm%ithermo == 0) &
    call transpose_x_to_y (qy_xppc,         qy_xppc_ypencil, dm%dppc)                    ! ypencil : x-mom, o thermal

    call transpose_x_to_y (fl%qy,           qy_ypencil,      dm%dcpc)                    ! y-pencil : y-mom, w+o thermal
    call transpose_y_to_z (qy_ypencil,      qy_zpencil,      dm%dcpc)                    ! z-pencil : y-mom, w+o thermal
    call Get_z_midp_C2P_3D(qy_zpencil, qy_zcpp_zpencil, dm, dm%ibcz(i,:), dm%fbcz(i, :)) ! z-pencil : y-mom, w+o thermal
    if(dm%ithermo == 0) &
    call transpose_z_to_y (qy_zcpp_zpencil, qy_zcpp_ypencil, dm%dcpp)                    ! y-pencil : z-mom, o thermal

    if ( dm%ithermo == 1) then
    call Get_y_midp_P2C_3D(qy_ypencil, accc_ypencil, dm, dm%ibcy(i,:)   ) ! intermediate, accc_ypencil = qy_yccc_ypencil
    call transpose_y_to_x (accc_ypencil,    qy_yccc,         dm%dccc)                    ! x-pencil : x-mom, w   thermal
    call transpose_y_to_z (accc_ypencil,    qy_yccc_zpencil, dm%dccc)                    ! z-pencil : z-mom, w   thermal
    end if
!-------------------------------------------------------------------------------
!    qz --> qz_xpcp(WO) --> qz_xpcp_ypencil(I) --> qz_xpcp_zpencil(O)
!     | --> qz_ypencil(WO) --> qz_ycpp_ypencil(WO) --> qz_ycpp_zpencil(O)
!                       |  --> qz_zpencil(WO) --> qz_zccc_zpencil(I) --> qz_zccc_ypencil(W) --> qz_zccc(W)
!-------------------------------------------------------------------------------
    i = 3 
    call Get_x_midp_C2P_3D(fl%qz, qz_xpcp, dm, dm%ibcx(i,:), dm%fbcx(i, :) ) ! x-pencil : z-mom, w+o   thermal
    if ( dm%ithermo == 0) then
    call transpose_x_to_y (qz_xpcp,         apcp_ypencil,    dm%dpcp)                    ! intermediate, apcp_ypencil = qz_xpcp_ypencil
    call transpose_y_to_z (apcp_ypencil,    qz_xpcp_zpencil, dm%dpcp)                    ! z-pencil : x-mom, o   thermal
    end if
  
    call transpose_x_to_y (fl%qz,           qz_ypencil,      dm%dccp)                    ! y-pencil : z-mom, w+o   thermal
    call Get_y_midp_C2P_3D(qz_ypencil, qz_ycpp_ypencil, dm, dm%ibcy(i,:), dm%fbcy(i, :)) ! y-pencil : z-mom, w+o   thermal
    if ( dm%ithermo == 0) &
    call transpose_y_to_z (qz_ycpp_ypencil, qz_ycpp_zpencil, dm%dcpp)                    ! z-pencil : y-mom, o   thermal
  
    call transpose_y_to_z (qz_ypencil,      qz_zpencil,      dm%dccp)                    ! z-pencil : z-mom, w+o   thermal
    if ( dm%ithermo == 1) then
    call Get_z_midp_P2C_3D(qz_zpencil, accc_zpencil, dm, dm%ibcz(i,:)   ) ! intermediate, accc_zpencil = qz_zccc_zpencil
    call transpose_z_to_y (accc_zpencil,    qz_zccc_ypencil, dm%dccc)                    ! y-pencil : y-mom, w   thermal
    call transpose_y_to_x (qz_zccc_ypencil, qz_zccc,         dm%dccc)                    ! x-pencil : x-mom, w   thermal
    end if

    if ( dm%ithermo == 1) then
!-------------------------------------------------------------------------------
!    gx --> gx_ypencil(I) --> gx_yppc_ypencil(I)--> gx_yppc(W)
!                     |--> gx_zpencil(I) --> gx_zpcp_zpencil(I) --> gx_zpcp_ypencil(I) --> qx_zpcp
!-------------------------------------------------------------------------------
    i = 1 
    call transpose_x_to_y (fl%gx,        apcc_ypencil,       dm%dpcc)                   ! intermediate, apcc_ypencil = gx_ypencil
    call Get_y_midp_C2P_3D(apcc_ypencil, appc_ypencil, dm, dm%ibcy(i,:), dm%fbcy(i, :) ) ! intermediate, appc_ypencil = gx_yppc_ypencil
    call transpose_y_to_x (appc_ypencil, gx_yppc,            dm%dppc)                   ! x-pencil : y-mom, w   thermal
  
    call transpose_y_to_z (apcc_ypencil, apcc_zpencil,       dm%dpcc)                   ! intermediate, apcc_zpencil = gx_zpencil
    call Get_z_midp_C2P_3D(apcc_zpencil, apcp_zpencil, dm, dm%ibcz(i,:), dm%fbcz(i, :) ) ! intermediate, apcp_zpencil = gx_zpcp_zpencil
    call transpose_z_to_y (apcp_zpencil, apcp_ypencil,       dm%dpcp)                   ! intermediate, apcp_ypencil = gx_zpcp_ypencil
    call transpose_y_to_x (apcp_ypencil, gx_zpcp,            dm%dpcp)                   ! x-pencil : z-mom, wo  thermal
!-------------------------------------------------------------------------------
!    gy --> gy_ypencil(W) --> gy_zpencil(I) --> gy_zcpp_zpencil(I) --> gy_zcpp_ypencil(W)
!-------------------------------------------------------------------------------
    i = 2
    call transpose_x_to_y (fl%gy,        gy_ypencil,        dm%dcpc)                    ! y-pencil : y-mom, w   thermal
    call transpose_y_to_z (gy_ypencil,   acpc_zpencil,      dm%dcpc)                    ! intermediate, acpc_zpencil = gy_zpencil
    call Get_z_midp_C2P_3D(acpc_zpencil, acpp_zpencil, dm, dm%ibcz(i,:), dm%fbcz(i, :) ) ! intermediate, acpp_zpencil = gy_zcpp_zpencil
    call transpose_z_to_y (acpp_zpencil, gy_zcpp_ypencil,   dm%dcpp)                    ! y-pencil : z-mom, w   thermal
!-------------------------------------------------------------------------------
!    gz --> gz_xpcp(I)    --> gz_xpcp_ypencil(I) --> gz_xpcp_zpencil(W)
!     | --> gz_ypencil(I) --> gz_ycpp_ypencil(I) --> gz_ycpp_zpencil(W)
!-------------------------------------------------------------------------------
    i = 3 
    call Get_x_midp_C2P_3D(fl%gz, apcp, dm, dm%ibcx(i,:), dm%fbcx(i, :) ) ! intermediate, apcp = gz_xpcp
    call transpose_x_to_y (apcp,         apcp_ypencil,      dm%dpcp)                    ! intermediate  apcp_ypencil = gz_xpcp_ypencil
    call transpose_y_to_z (apcp_ypencil, gz_xpcp_zpencil,   dm%dpcp)                    ! z-pencil : x-mom, w   thermal
  
    call transpose_x_to_y (fl%gz,        accp,              dm%dccp)                    ! intermediate, accp = gz_ypencil
    call Get_y_midp_C2P_3D(accp,         acpp_ypencil, dm, dm%ibcy(i,:), dm%fbcy(i, :) ) ! intermediate, acpp_ypencil = gz_ycpp_ypencil
    call transpose_y_to_z (acpp_ypencil, gz_ycpp_zpencil,   dm%dcpp)                    ! z-pencil : y-mom, w   thermal
!-------------------------------------------------------------------------------
!   d --> d_ypencil --> d_zpencil
!-------------------------------------------------------------------------------
    call transpose_x_to_y (fl%dDens,       dDens_ypencil,    dm%dccc)                    ! y-pencil : y-mom, w   thermal
    call transpose_y_to_z (dDens_ypencil,  dDens_zpencil,    dm%dccc)                    ! z-pencil : z-mom, w   thermal
!-------------------------------------------------------------------------------
!    m --> dmdx_xpcc
!    | --> m_xpcc -->m_xpcc_ypencil -->dmdy_xpcc_ypencil-->dmdy_xpcc
!                                 | -->m_xpcc_zpencil --> dmdz_xpcc_zpencil--> dmdz_xpcc_ypencil(I) --> dmdz_xpcc 
!    | --> m_ypencil(I) --> m_ycpc_ypencil --> m_ycpc --> dmdx_ycpc
!                                        | --> m_ycpc_zpencil
!                     | --> m_zpencil(I) --> dmdz_zccp_zpencil
!                                      | --> m_zccp_zpencil --> m_zccp_ypencil --> m_zccp --> dmdx_zccp --> dmdx_zccp_ypencil(I) --> dmdx_zccp_zpencil--> dmdy_zccp_ypencil
!-------------------------------------------------------------------------------
    i = 5
    call Get_x_1st_derivative_C2P_3D(fl%mVisc, dmdx_xpcc, dm, dm%ibcx(i,:), dm%fbc_vism(1, :) )  ! x-pencil : x-mom, w thermal
    call Get_x_midp_C2P_3D          (fl%mVisc, m_xpcc,    dm, dm%ibcx(i,:), dm%fbc_vism(1, :) )  ! x-pencil : x-mom, w thermal
    call transpose_x_to_y (m_xpcc, m_xpcc_ypencil, dm%dpcc)                            ! y-pencil : x-mom, w thermal
    
    call Get_y_1st_derivative_C2C_3D(m_xpcc_ypencil, dmdy_xpcc_ypencil, dm, dm%ibcy(i,:), dm%fbc_vism(2, :) )  ! y-pencil : x-mom, w thermal
    call transpose_y_to_x (dmdy_xpcc_ypencil, dmdy_xpcc,         dm%dpcc)                            ! x-pencil : x-mom, w thermal
    
    call transpose_y_to_z (m_xpcc_ypencil,    m_xpcc_zpencil,    dm%dpcc)                            ! z-pencil : x-mom, w thermal
    call Get_z_1st_derivative_C2C_3D(m_xpcc_zpencil, dmdz_xpcc_zpencil, dm, dm%ibcz(i,:), dm%fbc_vism(3, :)  )  ! z-pencil : x-mom, w thermal
    call transpose_z_to_y (dmdz_xpcc_zpencil, apcc_ypencil,      dm%dpcc)                            ! intermediate, apcc_ypencil = dmdz_xpcc_ypencil
    call transpose_y_to_x (apcc_ypencil,      dmdz_xpcc,         dm%dpcc)                            ! x-pencil : x-mom, w thermal

    call transpose_x_to_y (fl%mVisc, accc_ypencil, dm%dccc)                                          ! intermediate, accc_ypencil = m_ypencil 
    call Get_y_midp_C2P_3D(accc_ypencil, m_ycpc_ypencil, dm, dm%ibcy(i,:), dm%fbc_vism(2, :) )              ! y-pencil : y-mom, w thermal
    call transpose_y_to_z (m_ycpc_ypencil,    m_ycpc_zpencil,    dm%dcpc)                            ! z-pencil : y-mom, w thermal
    call transpose_y_to_x (m_ycpc_ypencil,    m_ycpc,            dm%dcpc)                            ! x-pencil : y-mom, w thermal
    call Get_x_1st_derivative_C2C_3D(m_ycpc,  dmdx_ycpc, dm, dm%ibcx(i,:), dm%fbcx(i,:))                ! x-pencil : y-mom, w thermal
      
    call transpose_y_to_z (accc_ypencil,      accc_zpencil,      dm%dccc)                             ! intermediate, accc_zpencil = m_zpencil
    call Get_z_1st_derivative_C2P_3D(accc_zpencil, dmdz_zccp_zpencil, dm, dm%ibcz(i,:), dm%fbc_vism(3, :) )   ! z-pencil : z-mom, w thermal
    call Get_z_midp_C2P_3D(accc_zpencil, m_zccp_zpencil, dm, dm%ibcz(i,:), dm%fbc_vism(3, :) )              ! z-pencil : z-mom, w thermal
    call transpose_z_to_y (m_zccp_zpencil,    m_zccp_ypencil,    dm%dccp)                             ! y-pencil : z-mom, w thermal
    call transpose_y_to_x (m_zccp_ypencil,    m_zccp,            dm%dccp)                             ! x-pencil : z-mom, w thermal
    call Get_x_1st_derivative_C2C_3D(m_zccp,  dmdx_zccp, dm, dm%ibcx(i,:), dm%fbcx(i,:) )                ! x-pencil : z-mom, w thermal
    call transpose_x_to_y (dmdx_zccp,         accp_ypencil,      dm%dccp)                             ! intermidate, accp_ypencil = dmdx_zccp_ypencil
    call transpose_y_to_z (accp_ypencil,      dmdx_zccp_zpencil, dm%dccp)                             ! z-pencil : z-mom, w thermal
    call Get_y_1st_derivative_C2C_3D(m_zccp_ypencil, dmdy_zccp_ypencil, dm, dm%ibcy(i,:), dm%fbc_vism(2, :) )   ! y-pencil : z-mom, w thermal
!-------------------------------------------------------------------------------
! calculate div(u_vec)
!-------------------------------------------------------------------------------
    div  = ZERO 
    accc = ZERO

    call Get_x_1st_derivative_P2C_3D(fl%qx,      accc,         dm, dm%ibcx(1, :) ) ! accc = d(qx)/d(x)_xccc
    div = div + accc ! = d(qx)/d(x)_ccc

    call Get_y_1st_derivative_P2C_3D(qy_ypencil, accc_ypencil, dm, dm%ibcy(2, :) ) ! accc_ypencil = d(qy)/(y)_yccc_ypencil
    call transpose_y_to_x (accc_ypencil, accc,         dm%dccc)                                   ! accc = d(qy)/d(y)_yccc
    div = div + accc ! = d(qx)/d(x)_ccc + d(qy)/d(y)_ccc

    call Get_z_1st_derivative_P2C_3D(qz_zpencil, accc_zpencil, dm, dm%ibcz(3, :) ) ! accc_zpencil = d(qz)/(z)_zccc_zpencil
    call transpose_z_to_y (accc_zpencil, accc_ypencil, dm%dccc)           ! accc_ypencil = d(qz)/(z)_zccc_ypencil
    call transpose_y_to_x (accc_ypencil, accc,         dm%dccc)           ! accc = d(qz)/d(z)_zccc
    div = div + accc ! = d(qx)/d(x)_ccc + d(qy)/d(y)_ccc + d(qz)/d(z)_ccc

    call transpose_x_to_y (div,          div_ypencil,  dm%dccc)
    call transpose_y_to_z (div_ypencil,  div_zpencil,  dm%dccc)
    end if
!===============================================================================
! the RHS of momentum equation
! x-pencil : the RHS terms of all 3 momentum equations (derivative) operating in the x direction
!===============================================================================
    fl%mx_rhs = ZERO
    fl%my_rhs = ZERO
    fl%mz_rhs = ZERO
    mx_rhs_implicit  = ZERO

    apcc = ZERO
    acpc = ZERO
    accp = ZERO
!-------------------------------------------------------------------------------
! X-pencil : X-mom convection term (x-c1/3): d(gx * qx)/dx at (i', j, k)
!-------------------------------------------------------------------------------
    
    if ( dm%ithermo == 0) then
      fbc(1:2) = dm%fbcx(1, 1:2) * dm%fbcx(1, 1:2)
      call Get_x_1st_derivative_P2P_3D(-fl%qx * fl%qx, apcc, dm, dm%ibc(1, :), fbc(:) )
    end if
    if ( dm%ithermo == 1) then
      fbc(1:2) = dm%fbcx(1, 1:2) * dm%fbcx(1, 1:2) * dm%fbc_dend(1, 1:2)
      call Get_x_1st_derivative_P2P_3D(-fl%gx * fl%qx, apcc, dm, dm%ibcx(1, :), fbc(:) )
    end if
    fl%mx_rhs = fl%mx_rhs + apcc
!-------------------------------------------------------------------------------
! X-pencil : X-mom diffusion term (x-v1-1/7), \mu^x * L11(ux) at (i', j, k)
!-------------------------------------------------------------------------------
    call Get_x_2nd_derivative_P2P_3D(fl%qx, apcc, dm, dm%ibcx(1, :) )
    if ( dm%ithermo == 0) fl%mx_rhs = fl%mx_rhs +          fl%rre * apcc
    if ( dm%ithermo == 1) fl%mx_rhs = fl%mx_rhs + m_xpcc * fl%rre * apcc
!-------------------------------------------------------------------------------
! X-pencil : Y-mom convection term (y-c1/3), d(gx^y * qy^x)/dx at (i, j', k)
!-------------------------------------------------------------------------------
    if ( dm%ithermo == 0) then
      fbc(1:2) = dm%fbcx(1, 1:2) * dm%fbcx(2, 1:2)
      call Get_x_1st_derivative_P2C_3D( -qx_yppc * qy_xppc, acpc, dm, dm%ibcx(1, :) )
    end if
    if ( dm%ithermo == 1) then
      fbc(1:2) = dm%fbcx(2, 1:2) * dm%fbcx(1, 1:2) * dm%fbc_dend(1, 1:2)
      call Get_x_1st_derivative_P2C_3D( -gx_yppc * qy_xppc, acpc, dm, dm%ibcx(1, :) )
    end if
    fl%my_rhs = fl%my_rhs + acpc
!-------------------------------------------------------------------------------
! X-pencil : Y-mom diffusion term (y-v1-1/1), \mu * L11(uy) at (i, j', k)
!-------------------------------------------------------------------------------
    call Get_x_2nd_derivative_P2P_3D(fl%qy, acpc, dm, dm%ibcx(2, :) )
    if ( dm%ithermo == 0) fl%my_rhs = fl%my_rhs +          fl%rre * acpc
    if ( dm%ithermo == 1) fl%my_rhs = fl%my_rhs + m_ycpc * fl%rre * acpc
!-------------------------------------------------------------------------------
! X-pencil : Z-mom convection term (z-c1/3), d(gx^z * qz^x)/dx at (i, j, k')
!-------------------------------------------------------------------------------
    if ( dm%ithermo == 0) then
      fbc(1:2) = dm%fbcx(1, 1:2) * dm%fbcx(3, 1:2)
      call Get_x_1st_derivative_P2C_3D( -qx_zpcp * qz_xpcp, accp, dm, dm%ibcx(1, :) )
    end if
    if ( dm%ithermo == 1) then
      fbc(1:2) = dm%fbcx(1, 1:2) * dm%fbcx(3, 1:2) * dm%fbc_dend(1, 1:2)
      call Get_x_1st_derivative_P2C_3D( -gx_zpcp * qz_xpcp, accp, dm, dm%ibcx(1, :) )
    end if
    fl%mz_rhs = fl%mz_rhs + accp
!-------------------------------------------------------------------------------
! X-pencil : Z-mom diffusion term (z-v1-1/1), \mu * L11(uz) at (i, j, k')
!-------------------------------------------------------------------------------
    call Get_x_2nd_derivative_P2P_3D(fl%qz, accp, dm, dm%ibcx(3, :) )
    if ( dm%ithermo == 0) fl%mz_rhs = fl%mz_rhs +          fl%rre * accp
    if ( dm%ithermo == 1) fl%mz_rhs = fl%mz_rhs + m_zccp * fl%rre * accp
!-------------------------------------------------------------------------------
! X-pencil : X-mom, pressure gradient in x direction, sigma_1*d(p)/dx
!-------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3D(fl%pres, apcc, dm, dm%ibcx(4, :), dm%fbcx(4, :) )
    mx_rhs_implicit =  mx_rhs_implicit - sigma1p * apcc

    if ( dm%ithermo == 1) then
!-------------------------------------------------------------------------------
! x-pencil : X-mom, gravity force in x direction
!-------------------------------------------------------------------------------
      if(igravity == 1 .or. igravity == -1)  then
        call Get_x_midp_C2P_3D(fl%dDens, apcc, dm, dm%ibcx(5, :), dm%fbc_dend(1, :) )
        mx_rhs_implicit =  mx_rhs_implicit + fl%fgravity * apcc
      end if
!-------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v2/7), \mu^x * 1/3 * d (div)/dx at (i', j, k)
!-------------------------------------------------------------------------------
      ! check boundary of div for inlet b.c.
      call Get_x_1st_derivative_C2P_3D(div, apcc, dm, dm%ibcx(1, :), dm%fbcx(1, :) )
      fl%mx_rhs = fl%mx_rhs + one_third_rre * m_xpcc * apcc
!-------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v3/7), -2/3 * d\mu/dx * (div(u)^x) +  
!                                                    2 * d\mu/dx * du/dx
!-------------------------------------------------------------------------------
      call Get_x_midp_C2P_3D          (div, apcc, dm, dm%ibcx(1, :), dm%fbcx(1, :) )
      fl%mx_rhs = fl%mx_rhs - two_third_rre * dmdx_xpcc * apcc
      call Get_x_1st_derivative_P2P_3D(fl%qx, apcc, dm, dm%ibcx(1, :), dm%fbcx(1, :) )
      fl%mx_rhs = fl%mx_rhs + two_rre       * dmdx_xpcc * apcc
!-------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v4/7), d(mu^x)/dy * d(qy^y)/dx at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3D(qy_yccc, apcc, dm, dm%ibcx(2, :), dm%fbcx(2, :) )
      fl%mx_rhs =  fl%mx_rhs + fl%rre * dmdy_xpcc * apcc
!-------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v5/7), d(mu^x)/dz * d(qz^z)/dx at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3D(qz_zccc, apcc, dm, dm%ibcx(3, :), dm%fbcx(3, :) )
      fl%mx_rhs =  fl%mx_rhs + fl%rre * dmdz_xpcc * apcc
!-------------------------------------------------------------------------------
!   X-pencil : Y-mom diffusion term (y-v6/7), d(mu^y)/dx * d(qy^x))/dx at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2C_3D(fl%qy, acpc, dm, dm%ibcx(2, :), dm%fbcx(2, :) )
      fl%my_rhs =  fl%my_rhs + fl%rre * dmdx_ycpc * acpc
!-------------------------------------------------------------------------------
!   X-pencil : Z-mom diffusion term (z-v6/7), d(mu^z)/dx * d(qz)/dx at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2C_3D(fl%qz, accp, dm, dm%ibcx(3, :), dm%fbcx(3, :) )
      fl%mz_rhs =  fl%mz_rhs + fl%rre * dmdx_zccp * accp
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
    if ( dm%ithermo == 0) then
      fbc(1:2) = dm%fbcy(1, 1:2) * dm%fbcy(2, 1:2)
      call Get_y_1st_derivative_P2C_3D( -qy_xppc_ypencil * qx_yppc_ypencil, apcc_ypencil, dm, dm%ibcy(2, :) )
    end if
    if ( dm%ithermo == 1) then
      fbc(1:2) = dm%fbcy(1, 1:2) * dm%fbcy(2, 1:2) * dm%fbc_dend(2, 1:2)
      call Get_y_1st_derivative_P2C_3D( -gy_xppc_ypencil * qx_yppc_ypencil, apcc_ypencil, dm, dm%ibcy(2, :) )
    end if
    mx_rhs_ypencil = mx_rhs_ypencil + apcc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : X-mom diffusion term (x-v1-2/7), \mu^x * L22(ux) at (i', j, k)
!-------------------------------------------------------------------------------
    call Get_y_2nd_derivative_P2P_3D(qx_ypencil, apcc_ypencil, dm, dm%ibcy(1, :), dm%fbcy(1, :) )
    if ( dm%ithermo == 0) mx_rhs_ypencil = mx_rhs_ypencil +                  fl%rre * apcc_ypencil
    if ( dm%ithermo == 1) mx_rhs_ypencil = mx_rhs_ypencil + m_xpcc_ypencil * fl%rre * apcc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom convection term (y-c2/3), d(gy * qy)/dy at (i, j', k)
!-------------------------------------------------------------------------------
    if ( dm%ithermo == 0) then
      fbc(1:2) = dm%fbcy(2, 1:2) * dm%fbcy(2, 1:2)
      call Get_y_1st_derivative_P2P_3D(-qy_ypencil * qy_ypencil, acpc_ypencil, dm, dm%ibcy(2, :), fbc(:) )
    end if
    if ( dm%ithermo == 1) then
      fbc(1:2) = dm%fbcy(2, 1:2) * dm%fbcy(2, 1:2) * dm%fbc_dend(2, 1:2)
      call Get_y_1st_derivative_P2P_3D(-gy_ypencil * qy_ypencil, acpc_ypencil, dm, dm%ibcy(2, :), fbc(:) )
    end if
    my_rhs_ypencil = my_rhs_ypencil + acpc_ypencil

!-------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v1-2/1), \mu * L22(uy) at (i, j', k)
!-------------------------------------------------------------------------------
    call Get_y_2nd_derivative_P2P_3D(qy_ypencil, acpc_ypencil, dm, dm%ibcy(2, :), dm%fbcy(2, :) )
    if ( dm%ithermo == 0) my_rhs_ypencil = my_rhs_ypencil +                  fl%rre * acpc_ypencil
    if ( dm%ithermo == 1) my_rhs_ypencil = my_rhs_ypencil + m_ycpc_ypencil * fl%rre * acpc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Z-mom convection term (z-c2/3), d(gy^z * qz^y)/dy at (i, j, k')
!-------------------------------------------------------------------------------
    if ( dm%ithermo == 0) then
      fbc(1:2) = dm%fbcy(2, 1:2) * dm%fbcy(3, 1:2)
      call Get_y_1st_derivative_P2C_3D( -qy_zcpp_ypencil * qz_ycpp_ypencil, accp_ypencil, dm, dm%ibcy(2, :) )
    end if
    if ( dm%ithermo == 1) then
      fbc(1:2) = dm%fbcy(3, 1:2) * dm%fbcy(2, 1:2) * dm%fbc_dend(2, 1:2)
      call Get_y_1st_derivative_P2C_3D( -gy_zcpp_ypencil * qz_ycpp_ypencil, accp_ypencil, dm, dm%ibcy(2, :) )
    end if
    mz_rhs_ypencil = mz_rhs_ypencil + accp_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Z-mom diffusion term (z-v1-2/1), \mu * L22(uz) at (i, j, k')
!-------------------------------------------------------------------------------
    call Get_y_2nd_derivative_P2P_3D( qz_ypencil, accp_ypencil, dm, dm%ibcy(3, :), dm%fbcy(3, :) )
    if ( dm%ithermo == 0) mz_rhs_ypencil = mz_rhs_ypencil +                  fl%rre * accp_ypencil
    if ( dm%ithermo == 1) mz_rhs_ypencil = mz_rhs_ypencil + m_zccp_ypencil * fl%rre * accp_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom pressure gradient in y direction, d(sigma_1 p)
!-------------------------------------------------------------------------------
    call Get_y_1st_derivative_C2P_3D(pres_ypencil, acpc_ypencil, dm, dm%ibcy(4, :), dm%fbcy(4, :) )
    my_rhs_implicit_ypencil =  my_rhs_implicit_ypencil - sigma1p * acpc_ypencil

    if ( dm%ithermo == 1) then
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom gravity force in y direction
!-------------------------------------------------------------------------------
      if( igravity == 2 .or. igravity == -2) then
        call Get_y_midp_C2P_3D(dDens_ypencil, acpc_ypencil, dm, dm%ibcy(5, :), dm%fbc_dend(2, 1:2) )
        my_rhs_implicit_ypencil =  my_rhs_implicit_ypencil + fl%fgravity * acpc_ypencil
      end if
!-------------------------------------------------------------------------------
! Y-pencil : X-mom diffusion term (x-v6/7), d(mu^x)/dy * d(qx)/dy at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2C_3D(qx_ypencil, apcc_ypencil, dm, dm%ibcy(1, :), dm%fbcy(1, :) )
      mx_rhs_ypencil =  mx_rhs_ypencil + fl%rre * dmdy_xpcc_ypencil * apcc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v2/7), \mu^y * 1/3 * d (div)/dy at (i, j', k)
!-------------------------------------------------------------------------------
      ! check b.c for div, for inlet/outlet
      call Get_y_1st_derivative_C2P_3D(div_ypencil, acpc_ypencil, dm, dm%ibcy(2, :), dm%fbcy(2, :) )
      my_rhs_ypencil = my_rhs_ypencil + one_third_rre * m_ycpc_ypencil * acpc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v3/7), 2d\mu/dy * (-1/3 * div(u)) +  2d\mu/dy * dv/dy
!-------------------------------------------------------------------------------
      call Get_y_midp_C2P_3D          (div_ypencil, acpc_ypencil, dm, dm%ibcy(2, :), dm%fbcy(2, :) )
      my_rhs_ypencil = my_rhs_ypencil - two_third_rre * dmdy_ycpc_ypencil * acpc_ypencil
      call Get_y_1st_derivative_P2P_3D(qy_ypencil,  acpc_ypencil, dm, dm%ibcy(2, :), dm%fbcy(2, :) )
      my_rhs_ypencil = my_rhs_ypencil + two_rre       * dmdy_ycpc_ypencil * acpc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v4/7), d(mu^y)/dx * d(qx^x)/dy at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2P_3D(qx_xccc_ypencil, acpc_ypencil, dm, dm%ibcy(1, :), dm%fbcy(1, :) )
      my_rhs_ypencil =  my_rhs_ypencil + fl%rre * dmdx_ycpc_ypencil * acpc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v6/7), d(mu)/dz * d(qz^z)/dy at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2P_3D(qz_zccc_ypencil, acpc_ypencil, dm, dm%ibcy(3, :), dm%fbcy(3, :) )
      my_rhs_ypencil =  my_rhs_ypencil + fl%rre * dmdz_ycpc_ypencil * acpc_ypencil
!-------------------------------------------------------------------------------
! Y-pencil : Z-mom diffusion term (z-v7/7), d(mu^z)/dy * d(qz)/dy at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2C_3D(qz_ypencil, accp_ypencil, dm, dm%ibcy(3, :), dm%fbcy(3, :) )
      mz_rhs_ypencil =  mz_rhs_ypencil + fl%rre * dmdy_zccp_ypencil * accp_ypencil
    end if
!-------------------------------------------------------------------------------
! Data from Y-pencil to X-pencil
!-------------------------------------------------------------------------------    
    call transpose_from_y_to_x (mx_rhs_ypencil, apcc, dm%dpcc)
    call transpose_from_y_to_x (my_rhs_ypencil, acpc, dm%dcpc)
    call transpose_from_y_to_x (mz_rhs_ypencil, accp, dm%dccp)
    fl%mx_rhs = fl%mx_rhs + apcc
    fl%my_rhs = fl%my_rhs + acpc
    fl%mz_rhs = fl%mz_rhs + accp
    call transpose_from_y_to_x (my_rhs_implicit_ypencil, acpc, dm%dcpc)
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
    if ( dm%ithermo == 0) then
      fbc(1:2) = dm%fbcz(1, 1:2) * dm%fbcz(3, 1:2)
      call Get_z_1st_derivative_P2C_3D( -qz_xpcp_zpencil * qx_zpcp_zpencil, apcc_zpencil, dm, dm%ibcz(1, :) )
    end if
    if ( dm%ithermo == 1) then
      fbc(1:2) = dm%fbcz(1, 1:2) * dm%fbcz(3, 1:2) * dm%fbc_dend(3, 1:2)
      call Get_z_1st_derivative_P2C_3D( -gz_xpcp_zpencil * qx_zpcp_zpencil, apcc_zpencil, dm, dm%ibcz(1, :) )
    end if
    mx_rhs_zpencil = mx_rhs_zpencil + apcc_zpencil

!-------------------------------------------------------------------------------
! Z-pencil : X-mom diffusion term (x-v1-3/7), \mu^x * L33(ux) at (i', j, k)
!-------------------------------------------------------------------------------
    call Get_z_2nd_derivative_P2P_3D(qx_zpencil, apcc_zpencil, dm, dm%ibcz(1, :), dm%fbcz(1, :) )
    if ( dm%ithermo == 0) mx_rhs_zpencil = mx_rhs_zpencil +                  fl%rre * apcc_zpencil
    if ( dm%ithermo == 1) mx_rhs_zpencil = mx_rhs_zpencil + m_xpcc_zpencil * fl%rre * apcc_zpencil

!-------------------------------------------------------------------------------
! Z-pencil : Y-mom convection term (y-c3/3), d(<gz>^y * <qy>^z)/dz at (i, j', k)
!-------------------------------------------------------------------------------
    if ( dm%ithermo == 0) then
      fbc(1:2) = dm%fbcz(2, 1:2) * dm%fbcz(3, 1:2)
      call Get_z_1st_derivative_P2C_3D( -qz_ycpp_zpencil * qy_zcpp_zpencil, acpc_zpencil, dm, dm%ibcz(3, :) )
    end if
    if ( dm%ithermo == 1) then
      fbc(1:2) = dm%fbcz(2, 1:2) * dm%fbcz(3, 1:2) * dm%fbc_dend(3, 1:2)
      call Get_z_1st_derivative_P2C_3D( -gz_ycpp_zpencil * qy_zcpp_zpencil, acpc_zpencil, dm, dm%ibcz(3, :) )
    end if
    my_rhs_zpencil = my_rhs_zpencil + acpc_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Y-mom diffusion term (y-v1-3/7), \mu * L33(uy) at (i, j', k)
!-------------------------------------------------------------------------------
    call Get_z_2nd_derivative_P2P_3D(qy_zpencil, acpc_zpencil, dm, dm%ibcz(2, :), dm%fbcz(2, :) )
    if ( dm%ithermo == 0) my_rhs_zpencil = my_rhs_zpencil +                  fl%rre * acpc_zpencil
    if ( dm%ithermo == 1) my_rhs_zpencil = my_rhs_zpencil + m_ycpc_zpencil * fl%rre * acpc_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v1/1), \mu * L33(uz) at (i, j, k')
!-------------------------------------------------------------------------------
    call Get_z_2nd_derivative_P2P_3D(qz_zpencil, accp_zpencil, dm, dm%ibcz(3, :), dm%fbcz(3, :) )
    if ( dm%ithermo == 0) mz_rhs_zpencil = mz_rhs_zpencil +                  fl%rre * accp_zpencil
    if ( dm%ithermo == 1) mz_rhs_zpencil = mz_rhs_zpencil + m_zccp_zpencil * fl%rre * accp_zpencil
!-------------------------------------------------------------------------------
! z-pencil : pressure gradient in z direction, d(sigma_1 p)
!-------------------------------------------------------------------------------
    call Get_z_1st_derivative_C2P_3D(pres_zpencil, accp_zpencil, dm, dm%ibcz(4, :), dm%fbcz(4, :) )
    mz_rhs_implicit_zpencil =  mz_rhs_implicit_zpencil - sigma1p * accp_zpencil

    if ( dm%ithermo == 1) then
!-------------------------------------------------------------------------------
! z-pencil : gravity force in z direction
!-------------------------------------------------------------------------------
      if( igravity == 3 .or. igravity == -3) then
        call Get_z_midp_C2P_3D(dDens_zpencil, accp_zpencil, dm, dm%ibcz(5, :), dm%fbc_dend(3, :) )
        mz_rhs_implicit_zpencil =  mz_rhs_implicit_zpencil + fl%fgravity * accp_zpencil
      end if
!-------------------------------------------------------------------------------
! Z-pencil : X-mom diffusion term (x-v7/7), d(mu^x)/dz * d(qx)/dz at (i', j, k)
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2C_3D(qx_zpencil, apcc_zpencil, dm, dm%ibcz(1, :), dm%fbcz(1, :) )
      mx_rhs_zpencil = mx_rhs_zpencil + fl%rre * dmdz_xpcc_zpencil * apcc_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Y-mom diffusion term (y-v7/7), d(mu^y)/dz * d(qy)/dz at (i, j', k)
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2C_3D(qy_zpencil, acpc_zpencil, dm, dm%ibcz(2, :), dm%fbcz(2, :) )
      my_rhs_zpencil =  my_rhs_zpencil + fl%rre * dmdz_ycpc_zpencil * acpc_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v2/7), \mu * 1/3 * d (div)/dz at (i, j, k')
!-------------------------------------------------------------------------------
      ! check div b.c. for inlet/outlet
      call Get_z_1st_derivative_C2P_3D(div_zpencil, accp_zpencil, dm, dm%ibcz(3, :), dm%fbcz(3, :) )
      mz_rhs_zpencil = mz_rhs_zpencil + one_third_rre * m_zccp_zpencil * accp_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v3/7), 2d\mu/dz * (-1/3 * div(u)) +  2d\mu/dz * dw/dz
!-------------------------------------------------------------------------------
      call Get_z_midp_C2P_3D(div_zpencil, accp_zpencil, dm,  dm%ibcz(3, :), dm%fbcz(3, :) )
      mz_rhs_zpencil = mz_rhs_zpencil - two_third_rre * dmdz_zccp_zpencil * accp_zpencil
      call Get_z_1st_derivative_P2P_3D( qz_zpencil,  accp_zpencil, dm, dm%ibcz(3, :), dm%fbcz(3, :) )
      mz_rhs_zpencil = mz_rhs_zpencil + two_rre * dmdz_zccp_zpencil * accp_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v4/7), d(mu^z)/dx * d(qx^x)/dz at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2P_3D(qx_xccc_zpencil, accp_zpencil, dm, dm%ibcz(1, :), dm%fbcz(1, :) )
      mz_rhs_zpencil = mz_rhs_zpencil + fl%rre * dmdx_zccp_zpencil * accp_zpencil
!-------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v6/7), d(mu^z)/dy * d(qy^y)/dz at (i, j, k')
!-------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2P_3D(qy_yccc_zpencil, accp_zpencil, dm, dm%ibcz(2, :), dm%fbcz(2, :) )
      mz_rhs_zpencil = mz_rhs_zpencil + fl%rre * dmdy_zccp_zpencil * accp_zpencil

    end if
!-------------------------------------------------------------------------------
! Data from Z-pencil to X-pencil
!-------------------------------------------------------------------------------    
    call transpose_from_z_to_y (mx_rhs_zpencil, apcc_ypencil, dm%dpcc)
    call transpose_from_z_to_y (my_rhs_zpencil, acpc_ypencil, dm%dcpc)
    call transpose_from_z_to_y (mz_rhs_zpencil, accp_ypencil, dm%dccp)
    
    call transpose_from_y_to_x (apcc_ypencil, apcc, dm%dpcc)
    call transpose_from_y_to_x (acpc_ypencil, acpc, dm%dcpc)
    call transpose_from_y_to_x (accp_ypencil, accp, dm%dccp)

    fl%mx_rhs = fl%mx_rhs + apcc
    fl%my_rhs = fl%my_rhs + acpc
    fl%mz_rhs = fl%mz_rhs + accp

    call transpose_from_z_to_y (mz_rhs_implicit_zpencil, accp_ypencil, dm%dccp)
    call transpose_from_y_to_x (accp_ypencil, accp, dm%dccp)

    mz_rhs_implicit = mz_rhs_implicit + accp
!===============================================================================
! x-pencil : to build up rhs in total, in all directions
!===============================================================================
!-------------------------------------------------------------------------------
! x-pencil : x-momentum
!-------------------------------------------------------------------------------
    if(fl%idriven /= IDRVF_NO) &
    call Calculate_xmomentum_driven_source(isub, fl%idriven, fl%drvfc, dm, fl%mx_rhs) 
    call Calculate_momentum_fractional_step(fl%mx_rhs0, fl%mx_rhs, mx_rhs_implicit, dm, isub)
!-------------------------------------------------------------------------------
! x-pencil : y-momentum
!-------------------------------------------------------------------------------
    call Calculate_momentum_fractional_step(fl%my_rhs0, fl%my_rhs, my_rhs_implicit, dm, isub)
!-------------------------------------------------------------------------------
! x-pencil : z-momentum
!-------------------------------------------------------------------------------
    call Calculate_momentum_fractional_step(fl%mz_rhs0, fl%mz_rhs, mz_rhs_implicit, dm, isub)
 
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
  subroutine Correct_massflux(ux, uy, uz, phi, dm, isub)
    use udf_type_mod,      only : t_domain
    use input_general_mod, only : tAlpha, dt, sigma2p
    use operations
    implicit none

    type(t_domain), intent(in   ) :: dm
    integer,     intent(in   ) :: isub
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ), intent(inout) :: ux
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ), intent(inout) :: uy
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ), intent(inout) :: uz
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ), intent(in   ) :: phi

    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: dphidx
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: dphidy
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: dphidz

    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: phi_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: dphidy_ypencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: dphidz_ypencil
    
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: phi_zpencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: dphidz_zpencil
    
!-------------------------------------------------------------------------------
!   x-pencil, ux = ux - dt * alpha * d(phi)/dx
!_______________________________________________________________________________
    call Get_x_1st_derivative_C2P_3D(phi,  dphidx, dm, dm%ibcx(4, :), dm%fbcx(4, :) )
    ux = ux - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidx
!-------------------------------------------------------------------------------
!   y-pencil, uy = uy - dt * alpha * d(phi)/dy
!_______________________________________________________________________________
    call transpose_from_x_to_y (phi, phi_ypencil, dm%dccc)
    call Get_y_1st_derivative_C2P_3D(phi_ypencil, dphidy_ypencil, dm, dm%ibcy(4, :), dm%fbcy(4, :) )
    call transpose_from_y_to_x (dphidy_ypencil, dphidy, dm%dcpc)
    uy = uy - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidy
!-------------------------------------------------------------------------------
!   z-pencil, uz = uz - dt * alpha * d(phi)/dz
!_______________________________________________________________________________
    call transpose_from_y_to_z (phi_ypencil, phi_zpencil, dm%dccc)
    call Get_z_1st_derivative_C2P_3D(phi_zpencil, dphidz_zpencil, dm, dm%ibcz(4, :), dm%fbcz(4, :) )
    call transpose_from_z_to_y (dphidz_zpencil, dphidz_ypencil, dm%dccp)
    call transpose_from_y_to_x (dphidz_ypencil, dphidz,         dm%dccp)
    uz = uz - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidz

    return
  end subroutine Correct_massflux
!===============================================================================
!> \brief To update the provisional u or rho u.
!>
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                         !
!______________________________________________________________________________!
!> \param[inout]  fl            flow field
!> \param[inout]  dm            domain
!> \param[in]     isub         RK sub-iteration
!===============================================================================
  subroutine Solve_momentum_eq(fl, dm, isub)
    use udf_type_mod,      only : t_flow, t_domain
    use typeconvert_mod
    use continuity_eq_mod
    use poisson_mod
    use boundary_conditions_mod
    implicit none

    type(t_flow), intent(inout) :: fl
    type(t_domain),  intent(in) :: dm
    integer,      intent(in) :: isub

!-------------------------------------------------------------------------------
! to calculate the rhs of the momenturn equation in stepping method
!_______________________________________________________________________________ 
    call Compute_momentum_rhs(fl, dm, isub)
!-------------------------------------------------------------------------------
! to update intermediate (\hat{q}) or (\hat{g})
!_______________________________________________________________________________
 
    if(iviscous == IVIS_EXPLICIT) then

      if ( dm%ithermo == 0) then 
        call Calculate_intermediate_mvar(fl%mx_rhs, fl%qx)
        call Calculate_intermediate_mvar(fl%my_rhs, fl%qy)
        call Calculate_intermediate_mvar(fl%mz_rhs, fl%qz)
      else
        call Calculate_intermediate_mvar(fl%mx_rhs, fl%gx)
        call Calculate_intermediate_mvar(fl%my_rhs, fl%gy)
        call Calculate_intermediate_mvar(fl%mz_rhs, fl%gz)
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
    call Apply_BC_velocity (dm, fl%qx, fl%qy, fl%qz)
    if(dm%ithermo ==  1) call Apply_BC_velocity (dm, fl%gx, fl%gy, fl%gz)
!-------------------------------------------------------------------------------
! to calculate the provisional divergence constrains
!-------------------------------------------------------------------------------
  !call Print_debug_mid_msg("  Computing provisional divergence constrains ...")
    call Calculate_continuity_constrains(fl, dm, isub)
!-------------------------------------------------------------------------------
! to solve Poisson equation
!-------------------------------------------------------------------------------
  !call Print_debug_mid_msg("  Solving Poisson Equation ...")
    call Solve_poisson(fl%pcor)
!-------------------------------------------------------------------------------
! to update velocity/massflux correction
!-------------------------------------------------------------------------------
  !call Print_debug_mid_msg("  Updating velocity/mass flux ...")
    if ( dm%ithermo == 0) then 
      call Correct_massflux(fl%qx, fl%qy, fl%qz, fl%pcor, dm, isub)
    else
      call Correct_massflux(fl%gx, fl%gy, fl%gz, fl%pcor, dm, isub)
    end if
!-------------------------------------------------------------------------------
! to update pressure
!-------------------------------------------------------------------------------
    fl%pres = fl%pres + fl%pcor
!-------------------------------------------------------------------------------
! to update b.c. values
!-------------------------------------------------------------------------------
    call Apply_BC_velocity (dm, fl%gx, fl%gy, fl%gz)
    call Apply_BC_velocity (dm, fl%qx, fl%qy, fl%qz)

    return
  end subroutine Solve_momentum_eq

end module eq_momentum_mod
