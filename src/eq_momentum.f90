module eq_momentum_mod
  use operations
  use precision_mod
  use decomp_2d
  implicit none

  private :: Calculate_momentum_fractional_step
  private :: Compute_momentum_rhs
  private :: Correct_massflux
  private :: solve_poisson
  private :: solve_poisson_x2z
  
  public  :: Solve_momentum_eq
  

contains
!==========================================================================================================
!==========================================================================================================
!> \brief To calcuate the convection and diffusion terms in rhs of momentum eq.
!>
!> This subroutine is called everytime when calcuting the rhs of momentum eqs.
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                         !
!______________________________________________________________________________!
!> \param[inout]  rhs0          the last iteration rhs
!> \param[inout]  rhs1          the current iteration rhs
!> \param[in]     rhs1_semi     the semi-implicit term
!> \param[in]     isub          the RK iteration to get correct Coefficient 
!_______________________________________________________________________________
  subroutine Calculate_momentum_fractional_step(rhs0, rhs1, rhs1_pfc, dtmp, dm, isub)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    type(t_domain), intent(in) :: dm
    real(WP), dimension(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3)), intent(in)    :: rhs1_pfc
    real(WP), dimension(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3)), intent(inout) :: rhs0, rhs1
    integer,  intent(in) :: isub
    
    real(WP) :: rhs_explicit_current, rhs_explicit_last, rhs_total
    integer :: i, j, k 

    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
        do i = 1, dtmp%xsz(1)

      ! add explicit terms : convection+viscous rhs
          rhs_explicit_current = rhs1(i, j, k) 
          rhs_explicit_last    = rhs0(i, j, k)
          rhs_total = dm%tGamma(isub) * rhs_explicit_current + &
                      dm%tZeta (isub) * rhs_explicit_last
          rhs0(i, j, k) = rhs_explicit_current

      ! add pressure gradient
          rhs_total = rhs_total + &
                      dm%tAlpha(isub) * rhs1_pfc(i, j, k)

      ! times the time step 
          rhs1(i, j, k) = dm%dt * rhs_total

        end do
      end do
    end do

    return
  end subroutine
!==========================================================================================================
!> \brief To calcuate all rhs of momentum eq.
!>
!> This subroutine is called everytime when calcuting the rhs of momentum eqs.
!>
!----------------------------------------------------------------------------------------------------------
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
    use udf_type_mod
    use operations
    use solver_tools_mod
    implicit none

    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(in ) :: dm
    integer,     intent(in ) :: isub
!----------------------------------------------------------------------------------------------------------
! common vars
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: qx_ypencil
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: qx_zpencil
    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: mx_rhs_ypencil
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: mx_rhs_zpencil
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: mx_rhs_pfc

    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: qy_ypencil
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: qy_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: my_rhs_ypencil ! 
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: my_rhs_zpencil ! 
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: my_rhs_pfc
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: my_rhs_pfc_ypencil

    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: qz_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: qz_zpencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: mz_rhs_ypencil !
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: mz_rhs_zpencil ! 
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: mz_rhs_pfc
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: mz_rhs_pfc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: mz_rhs_pfc_zpencil

    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: acpc
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: accp ! 

    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: apcc_ypencil
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: apcc_zpencil
    
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil !
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: acpc_zpencil !
    
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil !
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil ! 
    
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: pres_ypencil ! p
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: pres_zpencil ! p
    
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: qx_ppc_ypencil ! <ux>^y at (xp, yp, zc)
    
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qy_cpp_zpencil ! <uy>^z at (xc, yp, zp)
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qz_cpp_ypencil ! <uz>^y at (xc, yp, zp)

    real(WP), dimension( dm%dpcp%ysz(1), dm%dpcp%ysz(2), dm%dpcp%ysz(3) ) :: apcp_ypencil ! 
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: qz_pcp         ! <uz>^x at (xp, yc, zp) 
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: qx_pcp_zpencil ! <ux>^z at (xp, yc, zp)
    
!----------------------------------------------------------------------------------------------------------
! thermal == 0 only
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: qx_pcp ! <ux>^z at (xp, yc, zp)
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: qx_ppc ! <ux>^y at (xp, yp, zc)
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: qy_ppc ! <uy>^x at (xp, yp, zc)
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: qy_ppc_ypencil ! <uy>^x at (xp, yp, zc)
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qy_cpp_ypencil ! <uy>^z at (xc, yp, zp)
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: qz_pcp_zpencil ! <uz>^x at (xp, yc, zp)
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qz_cpp_zpencil ! <uz>^y at (xc, yp, zp)
!----------------------------------------------------------------------------------------------------------
! thermal == 1 only
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: apcc
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: apcp
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: appc
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: apcp_zpencil
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: appc_ypencil
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: acpp_ypencil
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: acpp_zpencil
    
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: qx_ccc           ! <ux>^x at (xc, yc, zc)
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qx_ccc_ypencil   ! <ux>^x at (xc, yc, zc)
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: qx_ccc_zpencil   ! <ux>^x at (xc, yc, zc)
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: qy_ccc           ! <uy>^y at (xc, yc, zc)
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qy_ccc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: qy_ccc_zpencil   ! <uy>^y at (xc, yc, zc)
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: qz_ccc           ! <uz>^z at (xc, yc, zc)
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qz_ccc_ypencil   ! <uz>^z at (xc, yc, zc), intermediate
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: qz_ccc_zpencil

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: gx_ccc           ! <gx>^x at (xc, yc, zc)
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: gx_ppc           ! <gx>^y at (xp, yp, zc)
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: gx_pcp           ! <gx>^z at (xp, yc, zp)
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: gy_ypencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: gy_ccc_ypencil
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: gy_cpp_ypencil   ! <gy>^z at (xc, yp, zp)
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: gy_ppc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: gz_zpencil   ! 
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: gz_pcp_zpencil   ! <gz>^x at (xp, yc, zp)
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: gz_cpp_zpencil   ! <gz>^y at (xc, yp, zp)
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: gz_ccc_zpencil   ! <gz>^y at (xc, yc, zc)
    
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: dDens_ypencil  ! d 
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: dDens_zpencil  ! d 

    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: m_pcc         ! <mu>^x       at (xp, yc, zc)  
    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: m_pcc_ypencil ! <mu>^x       at (xp, yc, zc)
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: m_pcc_zpencil ! <mu>^x       at (xp, yc, zc)
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: m_cpc         ! <mu>^y       at (xc, yp, zc)
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: m_cpc_ypencil ! <mu>^y       at (xc, yp, zc)
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: m_cpc_zpencil ! <mu>^y       at (xc, yp, zc)
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: m_ccp         ! <mu>^z       at (xc, yc, zp)
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: m_ccp_ypencil ! <mu>^z       at (xc, yc, zp)
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: m_ccp_zpencil ! <mu>^z       at (xc, yc, zp)
    
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: dmdx_pcc ! d( mu   )/dx at (xp, yc, zc)
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: dmdx_cpc ! d(<mu>^y)/dx at (xc, yp, zc)
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: dmdx_ccp ! d(<mu>^z)/dx at (xc, yc, zp)
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: dmdy_pcc ! d(<mu>^x)/dy at (xp, yc, zc)
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: dmdz_pcc ! d(<mu>^x)/dz at (xp, yc, zc)

    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: dmdx_cpc_ypencil ! d( mu   )/dx at (xc, yp, zc)
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: dmdx_ccp_zpencil ! d(<mu>^z)/dx at (xc, yc, zp)
    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: dmdy_pcc_ypencil ! d(<mu>^x)/dy at (xp, yc, zc)
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: dmdy_cpc_ypencil ! d( mu   )/dy at (xc, yp, zc)
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: dmdy_ccp_ypencil ! d(<mu>^z)/dy at (xc, yc, zp)
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: dmdy_ccp_zpencil ! d(<mu>^z)/dx at (xc, yc, zp)
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: dmdz_pcc_zpencil ! d(<mu>^x)/dz at (xp, yc, zc)
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: dmdz_cpc_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: dmdz_cpc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: dmdz_ccp_zpencil ! d( mu   )/dz at (xc, yc, zp)
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: div
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: div_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: div_zpencil
!----------------------------------------------------------------------------------------------------------
! others
!----------------------------------------------------------------------------------------------------------
    real(WP) :: one_third_rre, two_third_rre, two_rre
    real(WP) :: fbc(2)
    integer  :: i, j, k, m, jj
    real(WP) :: rhsx_bulk, rhsy_bulk, rhsz_bulk

#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_start_msg("Compute_momentum_rhs ...")
#endif

    one_third_rre = ONE_THIRD * fl%rre
    two_third_rre = TWO_THIRD * fl%rre
          two_rre = TWO * fl%rre
!==========================================================================================================
! variable preparation
! In the comments: 
!     I = Intermediate
!     O = no thermal fields
!     W = with thermal fields
!     WO = both, commom
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
!    p --> p_ypencil (y-mom, w+o) --> p_zpencil (z-mom, w+o)
!----------------------------------------------------------------------------------------------------------
    !write(*,*) nrank,  'test-0'
    call transpose_x_to_y (fl%pres,      pres_ypencil, dm%dccc) 
    call transpose_y_to_z (pres_ypencil, pres_zpencil, dm%dccc)
!----------------------------------------------------------------------------------------------------------
!    qx --> qx_ccc (OW) --> qx_ccc_ypencil(W) --> qx_ccc_zpencil(W)
!     | --> qx_ypencil(WO) --> qx_ppc_ypencil(WO) --> qx_ppc(O)
!                       |  --> qx_zpencil(WO) --> qx_pcp_zpencil(WO) --> qx_pcp_ypencil(I) --> qx_pcp(O)
!----------------------------------------------------------------------------------------------------------
    !write(*,*) nrank,  'test-1'
    i = 1 

    call Get_x_midp_P2C_3D(fl%qx, qx_ccc, dm, dm%ibcx(:, i) )

    call transpose_x_to_y (fl%qx, qx_ypencil, dm%dpcc)
    call transpose_y_to_z (qx_ypencil, qx_zpencil, dm%dpcc) ! qx_zpencil : x-mom, w+o thermal

    call Get_y_midp_C2P_3D(qx_ypencil, qx_ppc_ypencil, dm, dm%ibcy(:, i), dm%fbcy(:, i)) ! qx_ppc_ypencil : x-mom, w+o thermal

    call transpose_y_to_x (qx_ppc_ypencil, qx_ppc, dm%dppc) !

    
    call Get_z_midp_C2P_3D(qx_zpencil, qx_pcp_zpencil, dm, dm%ibcz(:, i), dm%fbcz(:, i)) ! qx_pcp_zpencil : x-mom, w+o thermal

    if(.not. dm%is_thermo) then
      call transpose_z_to_y (qx_pcp_zpencil, apcp_ypencil, dm%dpcp) ! intermediate, apcp_ypencil = qx_pcp_ypencil
      call transpose_y_to_x (apcp_ypencil,    qx_pcp,      dm%dpcp) ! qx_pcp : z-mom,  o  thermal
    end if 
    if(dm%is_thermo) then
      call transpose_x_to_y (qx_ccc,         qx_ccc_ypencil, dm%dccc) ! qx_ccc_ypencil : y-mom, w   thermal
      call transpose_y_to_z (qx_ccc_ypencil, qx_ccc_zpencil, dm%dccc) ! qx_ccc_zpencil : z-mom, w   thermal
    end if
!----------------------------------------------------------------------------------------------------------
!    qy--> qy_ppc(WO) --> qy_ppc_ypencil(O) 
!     |--> qy_ypencil(WO) --> qy_zpencil(WO) --> qy_cpp_zpencil(WO) --> qy_cpp_ypencil(O)
!                       | --> qy_ccc_ypencil(WO) --> qy_ccc(W)
!                                              | --> qy_ccc_zpencil(W)              
!----------------------------------------------------------------------------------------------------------
    !write(*,*) nrank,  'test-2'
    i = 2 
    call Get_x_midp_C2P_3D(fl%qy, qy_ppc, dm, dm%ibcx(:,i), dm%fbcx(:,i)) ! qy_ppc : y-mom, w+o thermal
    if(.not. dm%is_thermo) then
      call transpose_x_to_y (qy_ppc, qy_ppc_ypencil, dm%dppc) ! qy_ppc_ypencil : x-mom, o thermal
    end if

    call transpose_x_to_y (fl%qy,      qy_ypencil, dm%dcpc) ! qy_ypencil : y-mom, w+o thermal
    call transpose_y_to_z (qy_ypencil, qy_zpencil, dm%dcpc) ! qy_zpencil : y-mom, w+o thermal
    call Get_z_midp_C2P_3D(qy_zpencil, qy_cpp_zpencil, dm, dm%ibcz(:,i), dm%fbcz(:,i))  ! qy_cpp_zpencil : y-mom, w+o thermal

    call Get_y_midp_P2C_3D(qy_ypencil, qy_ccc_ypencil, dm, dm%ibcy(:,i)) !

    if(.not. dm%is_thermo) then
      call transpose_z_to_y (qy_cpp_zpencil, qy_cpp_ypencil, dm%dcpp) ! qy_cpp_ypencil : z-mom, o thermal
    end if

    if ( dm%is_thermo) then 
      call transpose_y_to_x (qy_ccc_ypencil, qy_ccc,         dm%dccc)! qy_ccc: x-mom, w   thermal
      call transpose_y_to_z (qy_ccc_ypencil, qy_ccc_zpencil, dm%dccc)! qy_ccc_zpencil : z-mom, w   thermal
    end if
!----------------------------------------------------------------------------------------------------------
!    qz --> qz_pcp(WO) --> qz_pcp_ypencil(I) --> qz_pcp_zpencil(O)
!     | --> qz_ypencil(WO) --> qz_cpp_ypencil(WO) --> qz_cpp_zpencil(O)
!                       |  --> qz_zpencil(WO) --> qz_ccc_zpencil(WO) --> qz_ccc_ypencil(W) --> qz_ccc(W)
!----------------------------------------------------------------------------------------------------------
    !write(*,*) nrank,  'test-3'
    i = 3 
    call Get_x_midp_C2P_3D(fl%qz, qz_pcp, dm, dm%ibcx(:,i), dm%fbcx(:,i)) ! x-pencil : z-mom, w+o   thermal
    if ( .not. dm%is_thermo) then
      call transpose_x_to_y (qz_pcp,      apcp_ypencil,    dm%dpcp) ! intermediate, apcp_ypencil = qz_pcp_ypencil
      call transpose_y_to_z (apcp_ypencil, qz_pcp_zpencil, dm%dpcp) ! qz_pcp_zpencil : x-mom, o   thermal
    end if
  
    call transpose_x_to_y (fl%qz, qz_ypencil, dm%dccp) ! qz_ypencil : z-mom, w+o   thermal
    call Get_y_midp_C2P_3D(qz_ypencil, qz_cpp_ypencil, dm, dm%ibcy(:,i), dm%fbcy(:,i)) ! qz_cpp_ypencil : z-mom, w+o   thermal
    if ( .not. dm%is_thermo) then
      call transpose_y_to_z (qz_cpp_ypencil, qz_cpp_zpencil, dm%dcpp) ! z-pencil : y-mom, o   thermal
    end if

    call transpose_y_to_z (qz_ypencil, qz_zpencil, dm%dccp) ! z-pencil : z-mom, w+o   thermal
    call Get_z_midp_P2C_3D(qz_zpencil, qz_ccc_zpencil, dm, dm%ibcz(:,i)) ! intermediate, accc_zpencil = qz_ccc_zpencil
    if ( dm%is_thermo) then
      call transpose_z_to_y (qz_ccc_zpencil, qz_ccc_ypencil, dm%dccc) ! y-pencil : y-mom, w   thermal
      call transpose_y_to_x (qz_ccc_ypencil, qz_ccc,         dm%dccc) ! x-pencil : x-mom, w   thermal
    end if

    if ( dm%is_thermo) then
!----------------------------------------------------------------------------------------------------------
!    gx --> gx_ccc
!      |--> gx_ypencil(I) --> gx_ppc_ypencil(I)--> gx_ppc(W)
!                     |--> gx_zpencil(I) --> gx_pcp_zpencil(I) --> gx_pcp_ypencil(I) --> qx_pcp
!----------------------------------------------------------------------------------------------------------
      i = 1 
      call Get_x_midp_P2C_3D(fl%gx, gx_ccc, dm, dm%ibcx(:, i) ) ! 

      call transpose_x_to_y (fl%gx,        apcc_ypencil,       dm%dpcc)                   ! intermediate, apcc_ypencil = gx_ypencil
      do m = 1, 2
        fbc(m) = dm%fbcy(m, 1) * dm%fbc_dend(m, 2)
      end do 
      call Get_y_midp_C2P_3D(apcc_ypencil, appc_ypencil, dm, dm%ibcy(:,i), fbc ) ! intermediate, appc_ypencil = gx_ppc_ypencil
      call transpose_y_to_x (appc_ypencil, gx_ppc,            dm%dppc)                   ! gx_ppc : y-mom, w   thermal
    
      call transpose_y_to_z (apcc_ypencil, apcc_zpencil,       dm%dpcc)                   ! intermediate, apcc_zpencil = gx_zpencil
      do m = 1, 2
        fbc(m) = dm%fbcz(m, 1) * dm%fbc_dend(m, 3)
      end do 
      call Get_z_midp_C2P_3D(apcc_zpencil, apcp_zpencil, dm, dm%ibcz(:,i), fbc ) ! intermediate, apcp_zpencil = gx_pcp_zpencil
      call transpose_z_to_y (apcp_zpencil, apcp_ypencil,       dm%dpcp)                   ! intermediate, apcp_ypencil = gx_pcp_ypencil
      call transpose_y_to_x (apcp_ypencil, gx_pcp,            dm%dpcp)                   ! x-pencil : z-mom, wo  thermal
!----------------------------------------------------------------------------------------------------------
!    gy --> gy_ypencil(W) --> gy_zpencil(I) --> gy_cpp_zpencil(I) --> gy_cpp_ypencil(W)
!                         --> gy_ppc(I) --> gy_ppc_ypencil(W)
!                         --> gy_ccc_ypencil
!----------------------------------------------------------------------------------------------------------
      i = 2
      call transpose_x_to_y (fl%gy, gy_ypencil, dm%dcpc)                    ! y-pencil : y-mom, w   thermal
      do m = 1, 2
        fbc(m) = dm%fbcx(m, 2) * dm%fbc_dend(m, 1) ! uy * d
      end do 
      call Get_y_midp_C2P_3D(gy_ypencil, gy_ccc_ypencil, dm, dm%ibcy(:,i), fbc )     ! 
      call Get_x_midp_C2P_3D(fl%gy,      appc,           dm, dm%ibcx(:,i), fbc )     ! 
      call transpose_x_to_y (appc,       gy_ppc_ypencil, dm%dppc)                    ! 
      do m = 1, 2
        fbc(m) = dm%fbcz(m, 2) * dm%fbc_dend(m, 3) ! uy * d
      end do 
      call transpose_y_to_z (gy_ypencil, acpc_zpencil, dm%dcpc)                    ! y-pencil : y-mom, w   thermal
      call Get_z_midp_C2P_3D(acpc_zpencil, acpp_zpencil,      dm, dm%ibcz(:,i), fbc ) ! intermediate, acpp_zpencil = gy_cpp_zpencil
      call transpose_z_to_y (acpp_zpencil, gy_cpp_ypencil,   dm%dcpp)                    ! y-pencil : z-mom, w   thermal
!----------------------------------------------------------------------------------------------------------
!    gz --> gz_pcp(I)    --> gz_pcp_ypencil(I) --> gz_pcp_zpencil(W)
!     | --> gz_ypencil(I) --> gz_cpp_ypencil(I) --> gz_cpp_zpencil(W)
!                | --> gz_zpencil | --> gz_ccc_zpencil
!----------------------------------------------------------------------------------------------------------
      i = 3
      do m = 1, 2
        fbc(m) = dm%fbcx(m, 3) * dm%fbc_dend(m, 1)
      end do 
      call Get_x_midp_C2P_3D(fl%gz, apcp, dm, dm%ibcx(:,i), fbc(:) ) ! intermediate, apcp = gz_pcp
      call transpose_x_to_y (apcp,         apcp_ypencil,      dm%dpcp)                    ! intermediate  apcp_ypencil = gz_pcp_ypencil
      call transpose_y_to_z (apcp_ypencil, gz_pcp_zpencil,   dm%dpcp)                    ! z-pencil : x-mom, w   thermal
    
      call transpose_x_to_y (fl%gz,        accp_ypencil,   dm%dccp)                    ! intermediate, accp = gz_ypencil
      call transpose_y_to_z (accp_ypencil, gz_zpencil,     dm%dccp)                    ! intermediate, accp = gz_ypencil
      call Get_z_midp_P2C_3D(gz_zpencil, gz_ccc_zpencil,   dm, dm%ibcz(:,i) ) ! intermediate, acpp_ypencil = gz_cpp_ypencil
      do m = 1, 2
        fbc(m) = dm%fbcy(m, 3) * dm%fbc_dend(m, 2)
      end do 
      call Get_y_midp_C2P_3D(accp_ypencil, acpp_ypencil,   dm, dm%ibcy(:,i), fbc ) ! intermediate, acpp_ypencil = gz_cpp_ypencil
      call transpose_y_to_z (acpp_ypencil, gz_cpp_zpencil,   dm%dcpp)                    ! z-pencil : y-mom, w   thermal
!----------------------------------------------------------------------------------------------------------
!   d --> d_ypencil --> d_zpencil
!----------------------------------------------------------------------------------------------------------
      call transpose_x_to_y (fl%dDens,       dDens_ypencil,    dm%dccc)                    ! y-pencil : y-mom, w   thermal
      call transpose_y_to_z (dDens_ypencil,  dDens_zpencil,    dm%dccc)                    ! z-pencil : z-mom, w   thermal
!----------------------------------------------------------------------------------------------------------
!    m --> dmdx_pcc
!    | --> m_pcc -->m_pcc_ypencil -->dmdy_pcc_ypencil-->dmdy_pcc
!                                 | -->m_pcc_zpencil --> dmdz_pcc_zpencil--> dmdz_pcc_ypencil(I) --> dmdz_pcc 
!    | --> m_ypencil(I) --> dmdy_cpc_ypencil
!                     | --> m_cpc_ypencil --> m_cpc --> dmdx_cpc --> dmdx_cpc_ypencil
!                                       | --> m_cpc_zpencil --> dmdz_cpc_zpencil --> dmdz_cpc_ypencil 
!                     | --> m_zpencil(I) --> dmdz_ccp_zpencil
!                                      | --> m_ccp_zpencil --> m_ccp_ypencil --> m_ccp --> dmdx_ccp --> dmdx_ccp_ypencil(I) --> dmdx_ccp_zpencil
!                                                                     |        --> dmdy_ccp_ypencil --> dmdy_ccp_zpencil
!----------------------------------------------------------------------------------------------------------
      i = 5
      call Get_x_1st_derivative_C2P_3D(fl%mVisc, dmdx_pcc, dm, dm%ibcx(:,i), dm%fbc_vism(:, 1) )  ! x-pencil : x-mom, w thermal
  
      call Get_x_midp_C2P_3D          (fl%mVisc, m_pcc,    dm, dm%ibcx(:,i), dm%fbc_vism(:, 1) )  ! x-pencil : x-mom, w thermal
      call transpose_x_to_y (m_pcc, m_pcc_ypencil, dm%dpcc)                            ! y-pencil : x-mom, w thermal
      
      call Get_y_1st_derivative_C2C_3D(m_pcc_ypencil, dmdy_pcc_ypencil, dm, dm%ibcy(:,i), dm%fbc_vism(:, 2) )  ! y-pencil : x-mom, w thermal
      call transpose_y_to_x (dmdy_pcc_ypencil, dmdy_pcc,         dm%dpcc)                            ! x-pencil : x-mom, w thermal
      
      call transpose_y_to_z (m_pcc_ypencil,    m_pcc_zpencil,    dm%dpcc)                            ! z-pencil : x-mom, w thermal
      call Get_z_1st_derivative_C2C_3D(m_pcc_zpencil, dmdz_pcc_zpencil, dm, dm%ibcz(:,i), dm%fbc_vism(:, 3)  )  ! z-pencil : x-mom, w thermal
      call transpose_z_to_y (dmdz_pcc_zpencil, apcc_ypencil,      dm%dpcc)                            ! intermediate, apcc_ypencil = dmdz_pcc_ypencil
      call transpose_y_to_x (apcc_ypencil,      dmdz_pcc,         dm%dpcc)                            ! x-pencil : x-mom, w thermal

      call transpose_x_to_y (fl%mVisc, accc_ypencil, dm%dccc)       
      call Get_y_1st_derivative_C2P_3D(accc_ypencil,   dmdy_cpc_ypencil, dm, dm%ibcy(:,i), dm%fbc_vism(:, 2))                ! x-pencil : y-mom, w thermal                                   ! intermediate, accc_ypencil = m_ypencil 
      call Get_y_midp_C2P_3D(accc_ypencil, m_cpc_ypencil, dm, dm%ibcy(:,i), dm%fbc_vism(:, 2) )              ! y-pencil : y-mom, w thermal
      call transpose_y_to_z (m_cpc_ypencil,    m_cpc_zpencil,       dm%dcpc)                            ! z-pencil : y-mom, w thermal
      call transpose_y_to_x (m_cpc_ypencil,    m_cpc,               dm%dcpc)                            ! x-pencil : y-mom, w thermal
      
      call Get_z_1st_derivative_C2C_3D (m_cpc_zpencil, dmdz_cpc_zpencil, dm, dm%ibcy(:,i), dm%fbc_vism(:, 3))
      call transpose_z_to_y (dmdz_cpc_zpencil,    dmdz_cpc_ypencil,       dm%dcpc) 
      
      call Get_x_1st_derivative_C2C_3D(m_cpc,  dmdx_cpc, dm, dm%ibcx(:,i), dm%fbc_vism(:, 1))                ! x-pencil : y-mom, w thermal
      call transpose_x_to_y (dmdx_cpc, dmdx_cpc_ypencil,    dm%dcpc)

      call transpose_y_to_z (accc_ypencil,      accc_zpencil,      dm%dccc)                             ! intermediate, accc_zpencil = m_zpencil
      call Get_z_1st_derivative_C2P_3D(accc_zpencil, dmdz_ccp_zpencil, dm, dm%ibcz(:,i), dm%fbc_vism(:, 3) )   ! z-pencil : z-mom, w thermal
      call Get_z_midp_C2P_3D(accc_zpencil, m_ccp_zpencil, dm, dm%ibcz(:,i), dm%fbc_vism(:, 3) )              ! z-pencil : z-mom, w thermal
      call transpose_z_to_y (m_ccp_zpencil,    m_ccp_ypencil,    dm%dccp)                             ! y-pencil : z-mom, w thermal
      call transpose_y_to_x (m_ccp_ypencil,    m_ccp,            dm%dccp)                             ! x-pencil : z-mom, w thermal
      call Get_x_1st_derivative_C2C_3D(m_ccp,  dmdx_ccp, dm, dm%ibcx(:,i), dm%fbc_vism(:, 1) )                ! x-pencil : z-mom, w thermal
      call transpose_x_to_y (dmdx_ccp,         accp_ypencil,      dm%dccp)                             ! intermidate, accp_ypencil = dmdx_ccp_ypencil
      call transpose_y_to_z (accp_ypencil,      dmdx_ccp_zpencil, dm%dccp)                             ! z-pencil : z-mom, w thermal
      call Get_y_1st_derivative_C2C_3D(m_ccp_ypencil, dmdy_ccp_ypencil, dm, dm%ibcy(:,i), dm%fbc_vism(:, 2) )   ! y-pencil : z-mom, w thermal
      call transpose_y_to_z (dmdy_ccp_ypencil, dmdy_ccp_zpencil, dm%dccp)
!----------------------------------------------------------------------------------------------------------
! calculate div(u_vec)
!----------------------------------------------------------------------------------------------------------
      div  = ZERO 
      accc = ZERO

      call Get_x_1st_derivative_P2C_3D(fl%qx,      accc,         dm, dm%ibcx(:, 1) ) ! accc = d(qx)/d(x)_ccc
      div = div + accc ! = d(qx)/d(x)_ccc

      call Get_y_1st_derivative_P2C_3D(qy_ypencil, accc_ypencil, dm, dm%ibcy(:, 2) ) ! accc_ypencil = d(qy)/(y)_ccc_ypencil
      call transpose_y_to_x (accc_ypencil, accc,         dm%dccc)                                   ! accc = d(qy)/d(y)_ccc
      div = div + accc ! = d(qx)/d(x)_ccc + d(qy)/d(y)_ccc

      call Get_z_1st_derivative_P2C_3D(qz_zpencil, accc_zpencil, dm, dm%ibcz(:, 3) ) ! accc_zpencil = d(qz)/(z)_ccc_zpencil
      call transpose_z_to_y (accc_zpencil, accc_ypencil, dm%dccc)           ! accc_ypencil = d(qz)/(z)_ccc_ypencil
      call transpose_y_to_x (accc_ypencil, accc,         dm%dccc)           ! accc = d(qz)/d(z)_ccc
      div = div + accc ! = d(qx)/d(x)_ccc + d(qy)/d(y)_ccc + d(qz)/d(z)_ccc

      call transpose_x_to_y (div,          div_ypencil,  dm%dccc)
      call transpose_y_to_z (div_ypencil,  div_zpencil,  dm%dccc)
    end if
!write(*,*) nrank,  'test-4'
!==========================================================================================================
! the RHS of x-momentum equation
!==========================================================================================================
    fl%mx_rhs = ZERO
    mx_rhs_ypencil = ZERO
    mx_rhs_zpencil = ZERO

    mx_rhs_pfc  = ZERO

    apcc = ZERO
    apcc_ypencil = ZERO
    apcc_zpencil = ZERO
!----------------------------------------------------------------------------------------------------------
! X-pencil : X-mom convection term (x-c1/3): d(gx * qx)/dx at (i', j, k)
!----------------------------------------------------------------------------------------------------------    
    if ( .not. dm%is_thermo) then
      do m = 1, 2
        fbc(m) = dm%fbcx(m, 1) * dm%fbcx(m, 1)
      end do
      call Get_x_1st_derivative_C2P_3D(-qx_ccc * qx_ccc, apcc, dm, dm%ibcx(:, 1), fbc(:) )
    end if

    if ( dm%is_thermo) then
      do m = 1, 2
        fbc(m) = dm%fbcx(m, 1) * dm%fbcx(m, 1) * dm%fbc_dend(m, 1)
      end do
      call Get_x_1st_derivative_C2P_3D(-gx_ccc * qx_ccc, apcc, dm, dm%ibcx(:, 1), fbc(:) )
    end if
    fl%mx_rhs = fl%mx_rhs + apcc
!----------------------------------------------------------------------------------------------------------
! Y-pencil : X-mom convection term (x-c2/3): d(<gy>^x * <qx>^y)/dy at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_y_1st_derivative_P2C_3D( -qy_ppc_ypencil * qx_ppc_ypencil, apcc_ypencil, dm, dm%ibcy(:, 1) )
    end if
    if ( dm%is_thermo) then
      call Get_y_1st_derivative_P2C_3D( -gy_ppc_ypencil * qx_ppc_ypencil, apcc_ypencil, dm, dm%ibcy(:, 1) )
    end if
    mx_rhs_ypencil = mx_rhs_ypencil + apcc_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : X-mom convection term (x-c3/3): d(<gz>^x * <qx>^z)/dz at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_z_1st_derivative_P2C_3D( -qz_pcp_zpencil * qx_pcp_zpencil, apcc_zpencil, dm, dm%ibcz(:, 1) )
    end if
    if ( dm%is_thermo) then
      call Get_z_1st_derivative_P2C_3D( -gz_pcp_zpencil * qx_pcp_zpencil, apcc_zpencil, dm, dm%ibcz(:, 1) )
    end if
    mx_rhs_zpencil = mx_rhs_zpencil + apcc_zpencil
!----------------------------------------------------------------------------------------------------------
! X-pencil : X-mom pressure gradient in x direction, d(sigma_1 p)
!----------------------------------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3D( fl%pres, apcc, dm, dm%ibcx(:, 4), dm%fbcx(:, 4) )
    mx_rhs_pfc=  mx_rhs_pfc - apcc
!----------------------------------------------------------------------------------------------------------
! X-pencil : X-mom diffusion term (x-v1-1/7), \mu^x * LL1(ux) at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    !call Get_x_2nd_derivative_P2P_3D(fl%qx, apcc, dm, dm%ibcx(:, 1)) ! check
    accc = ZERO
    apcc = ZERO
    call Get_x_1st_derivative_P2C_3D(fl%qx, accc, dm, dm%ibcx(:, 1))
    call Get_x_1st_derivative_C2P_3D(accc,  apcc, dm, dm%ibcx(:, 1), dm%fbcx(:, 1))

    if ( .not. dm%is_thermo) fl%mx_rhs = fl%mx_rhs +         fl%rre * apcc
    if ( dm%is_thermo) fl%mx_rhs = fl%mx_rhs + m_pcc * fl%rre * apcc
!----------------------------------------------------------------------------------------------------------
! Y-pencil : X-mom diffusion term (x-v1-2/7), \mu^x * LL2(ux) at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    !call Get_y_2nd_derivative_C2C_3D(qx_ypencil, apcc_ypencil, dm, dm%ibcy(:, 1), dm%fbcy(:, 1))
    appc_ypencil = ZERO
    apcc_ypencil = ZERO
    call Get_y_1st_derivative_C2P_3D(qx_ypencil,   appc_ypencil, dm, dm%ibcy(:, 1), dm%fbcy(:, 1))
    call Get_y_1st_derivative_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%ibcy(:, 1))

    if ( .not. dm%is_thermo ) mx_rhs_ypencil = mx_rhs_ypencil +                 fl%rre * apcc_ypencil
    if ( dm%is_thermo) mx_rhs_ypencil = mx_rhs_ypencil + m_pcc_ypencil * fl%rre * apcc_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : X-mom diffusion term (x-v1-3/7), \mu^x * LL3(ux) at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    !call Get_z_2nd_derivative_C2C_3D(qx_zpencil, apcc_zpencil, dm, dm%ibcz(:, 1), dm%fbcz(:, 1))
    apcp_zpencil = ZERO 
    apcc_zpencil = ZERO
    call Get_z_1st_derivative_C2P_3D(qx_zpencil,   apcp_zpencil, dm, dm%ibcz(:, 1), dm%fbcz(:, 1))
    call Get_z_1st_derivative_P2C_3D(apcp_zpencil, apcc_zpencil, dm, dm%ibcz(:, 1))

    if ( .not. dm%is_thermo) mx_rhs_zpencil = mx_rhs_zpencil +                 fl%rre * apcc_zpencil
    if ( dm%is_thermo) mx_rhs_zpencil = mx_rhs_zpencil + m_pcc_zpencil * fl%rre * apcc_zpencil

    if ( dm%is_thermo) then
!----------------------------------------------------------------------------------------------------------
! x-pencil : X-mom, gravity force in x direction
!----------------------------------------------------------------------------------------------------------
      if(fl%igravity == 1 .or. fl%igravity == -1)  then
        call Get_x_midp_C2P_3D(fl%dDens, apcc, dm, dm%ibcx(:, 5), dm%fbc_dend(:, 1) )
        mx_rhs_pfc =  mx_rhs_pfc + fl%fgravity(1) * apcc
      end if
!----------------------------------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v2/7), \mu^x * 1/3 * d (div)/dx at (i', j, k)
!----------------------------------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3D(div, apcc, dm, dm%ibcx(:, 1) ) ! apcc = d(div)/dx_xpencil
      fl%mx_rhs = fl%mx_rhs + one_third_rre * m_pcc * apcc
!----------------------------------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v3/7), -2/3 * d\mu/dx * (div(u)^x) +  
!                                                    2 * d\mu/dx * du/dx
!----------------------------------------------------------------------------------------------------------
      fbc(:) = ZERO ! check
      call Get_x_midp_C2P_3D (div, apcc, dm, dm%ibcx(:, 1), fbc ) ! apcc = div_pcc
      fl%mx_rhs = fl%mx_rhs - two_third_rre * dmdx_pcc * apcc
      call Get_x_1st_derivative_P2P_3D(fl%qx, apcc, dm, dm%ibcx(:, 1), dm%fbcx(:, 1) ) ! apcc = d(qx)/dx_pcc
      fl%mx_rhs = fl%mx_rhs + two_rre       * dmdx_pcc * apcc
!----------------------------------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v4/7), d(mu^x)/dy * d(qy^y)/dx at (i', j, k)
!----------------------------------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3D(qy_ccc, apcc, dm, dm%ibcx(:, 2), dm%fbcx(:, 2) ) !apcc = d(qy)/dx
      fl%mx_rhs =  fl%mx_rhs + fl%rre * dmdy_pcc * apcc
!----------------------------------------------------------------------------------------------------------
!   Y-pencil : X-mom diffusion term (x-v5/7), d(mu^x)/dy * d(qx)/dy at (i', j, k)
!----------------------------------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2C_3D(qx_ypencil, apcc_ypencil, dm, dm%ibcy(:, 1), dm%fbcy(:, 1) ) !apcc_ypencil = d(qx)/dy_ypencil
      mx_rhs_ypencil =  mx_rhs_ypencil + fl%rre * dmdy_pcc_ypencil * apcc_ypencil
!----------------------------------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v6/7), d(mu^x)/dz * d(qz^z)/dx at (i', j, k)
!----------------------------------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3D(qz_ccc, apcc, dm, dm%ibcx(:, 3), dm%fbcx(:, 3) ) ! apcc = d(qz)/dx
      fl%mx_rhs =  fl%mx_rhs + fl%rre * dmdz_pcc * apcc
!----------------------------------------------------------------------------------------------------------
!   Z-pencil : X-mom diffusion term (x-v7/7), d(mu^x)/dz * d(qx)/dz at (i', j, k)
!----------------------------------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2C_3D(qx_zpencil, apcc_zpencil, dm, dm%ibcz(:, 1), dm%fbcz(:, 1) ) ! apcc_zpencil = d(qx)/dz
      mx_rhs_zpencil = mx_rhs_zpencil + fl%rre * dmdz_pcc_zpencil * apcc_zpencil
    end if   
!----------------------------------------------------------------------------------------------------------
! x-mom: convert all terms to rhs
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_x (mx_rhs_ypencil, apcc, dm%dpcc)
    fl%mx_rhs =  fl%mx_rhs + apcc
    call transpose_z_to_y (mx_rhs_zpencil, apcc_ypencil, dm%dpcc)
    call transpose_y_to_x (apcc_ypencil, apcc, dm%dpcc)
    fl%mx_rhs =  fl%mx_rhs + apcc 
!write(*,*) nrank,  'test-5'
!==========================================================================================================
! the RHS of Y momentum equation
!==========================================================================================================
    fl%my_rhs = ZERO
    my_rhs_ypencil = ZERO
    my_rhs_zpencil = ZERO
    
    my_rhs_pfc  = ZERO
    my_rhs_pfc_ypencil = ZERO

    acpc = ZERO
    acpc_ypencil = ZERO
    acpc_zpencil = ZERO
!----------------------------------------------------------------------------------------------------------
! X-pencil : Y-mom convection term (y-c1/3), d(gx^y * qy^x)/dx at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_x_1st_derivative_P2C_3D( -qx_ppc * qy_ppc, acpc, dm, dm%ibcx(:, 2) )
    end if
    if ( dm%is_thermo) then
      call Get_x_1st_derivative_P2C_3D( -gx_ppc * qy_ppc, acpc, dm, dm%ibcx(:, 2) )
    end if
    fl%my_rhs = fl%my_rhs + acpc
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom convection term (y-c2/3), d(gy * qy)/dy at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      do m = 1, 2
        fbc(m) = dm%fbcy(m, 2) * dm%fbcy(m, 2)
      end do
      call Get_y_1st_derivative_C2P_3D(-qy_ccc_ypencil * qy_ccc_ypencil, acpc_ypencil, dm, dm%ibcy(:, 2), fbc(:) )
    end if
    if ( dm%is_thermo) then
      do m = 1, 2
        fbc(m) = dm%fbcy(m, 2) * dm%fbcy(m, 2) * dm%fbc_dend(m, 2)
      end do
      call Get_y_1st_derivative_C2P_3D(-gy_ccc_ypencil * qy_ccc_ypencil, acpc_ypencil, dm, dm%ibcy(:, 2), fbc(:) )
    end if
    my_rhs_ypencil = my_rhs_ypencil + acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Y-mom convection term (y-c3/3), d(<gz>^y * <qy>^z)/dz at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_z_1st_derivative_P2C_3D( -qz_cpp_zpencil * qy_cpp_zpencil, acpc_zpencil, dm, dm%ibcz(:, 2) )
    end if
    if ( dm%is_thermo) then
      call Get_z_1st_derivative_P2C_3D( -gz_cpp_zpencil * qy_cpp_zpencil, acpc_zpencil, dm, dm%ibcz(:, 2) )
    end if
    my_rhs_zpencil = my_rhs_zpencil + acpc_zpencil
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom pressure gradient in y direction, d(sigma_1 p)
!----------------------------------------------------------------------------------------------------------
    call Get_y_1st_derivative_C2P_3D(pres_ypencil, acpc_ypencil, dm, dm%ibcy(:, 4), dm%fbcy(:, 4) )
    my_rhs_pfc_ypencil =  my_rhs_pfc_ypencil - acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! X-pencil : Y-mom diffusion term (y-v1-1/7), \mu * LL1(uy) at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    !call Get_x_2nd_derivative_C2C_3D(fl%qy, acpc, dm, dm%ibcx(:, 2) )
    appc = ZERO
    acpc = ZERO
    call Get_x_1st_derivative_C2P_3D(fl%qy, appc, dm, dm%ibcx(:, 2), dm%fbcx(:, 2) )
    call Get_x_1st_derivative_P2C_3D(appc,  acpc, dm, dm%ibcx(:, 2) )
    if (.not. dm%is_thermo) fl%my_rhs = fl%my_rhs +         fl%rre * acpc
    if ( dm%is_thermo) fl%my_rhs = fl%my_rhs + m_cpc * fl%rre * acpc
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v1-2/7), \mu * LL2(uy) at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    !call Get_y_2nd_derivative_P2P_3D(qy_ypencil, acpc_ypencil, dm, dm%ibcy(:, 2))
    call Get_y_1st_derivative_P2C_3D(qy_ypencil,   accc_ypencil, dm, dm%ibcy(:, 2))
    call Get_y_1st_derivative_C2P_3D(accc_ypencil, acpc_ypencil, dm, dm%ibcy(:, 2), dm%fbcy(:, 2))
    if ( .not. dm%is_thermo ) my_rhs_ypencil = my_rhs_ypencil +                  fl%rre * acpc_ypencil
    if ( dm%is_thermo) my_rhs_ypencil = my_rhs_ypencil + m_cpc_ypencil * fl%rre * acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Y-mom diffusion term (y-v1-3/7), \mu * LL3(uy) at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    !call Get_z_2nd_derivative_C2C_3D(qy_zpencil, acpc_zpencil, dm, dm%ibcz(:, 2))
    call Get_z_1st_derivative_C2P_3D(qy_zpencil,   acpp_zpencil, dm, dm%ibcz(:, 2), dm%fbcz(:, 2))
    call Get_z_1st_derivative_P2C_3D(acpp_zpencil, acpc_zpencil, dm, dm%ibcz(:, 2))
    if ( .not. dm%is_thermo ) my_rhs_zpencil = my_rhs_zpencil +                  fl%rre * acpc_zpencil
    if ( dm%is_thermo) my_rhs_zpencil = my_rhs_zpencil + m_cpc_zpencil * fl%rre * acpc_zpencil

    if ( dm%is_thermo) then
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom gravity force in y direction
!----------------------------------------------------------------------------------------------------------
      if(fl%igravity == 2 .or. fl%igravity == -2) then
        call Get_y_midp_C2P_3D(dDens_ypencil, acpc_ypencil, dm, dm%ibcy(:, 5), dm%fbc_dend(:, 2) )
        my_rhs_pfc_ypencil =  my_rhs_pfc_ypencil + fl%fgravity(2) * acpc_ypencil
      end if
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v2/7), \mu^y * 1/3 * d (div)/dy at (i, j', k)
!----------------------------------------------------------------------------------------------------------
      fbc = ZERO ! check
      call Get_y_1st_derivative_C2P_3D(div_ypencil, acpc_ypencil, dm, dm%ibcy(:, 2), fbc)
      my_rhs_ypencil = my_rhs_ypencil + one_third_rre * m_cpc_ypencil * acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v3/7), 2d\mu/dy * (-1/3 * div(u)) +  2d\mu/dy * dv/dy
!----------------------------------------------------------------------------------------------------------
      fbc = ZERO ! check
      call Get_y_midp_C2P_3D (div_ypencil, acpc_ypencil, dm, dm%ibcy(:, 2), fbc )
      my_rhs_ypencil = my_rhs_ypencil - two_third_rre * dmdy_cpc_ypencil * acpc_ypencil
      call Get_y_1st_derivative_P2P_3D(qy_ypencil,  acpc_ypencil, dm, dm%ibcy(:, 2), dm%fbcy(:, 2) )
      my_rhs_ypencil = my_rhs_ypencil + two_rre       * dmdy_cpc_ypencil * acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v4/7), d(mu^y)/dx * d(qx^x)/dy at (i, j', k)
!----------------------------------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2P_3D(qx_ccc_ypencil, acpc_ypencil, dm, dm%ibcy(:, 1), dm%fbcy(:, 1) )
      my_rhs_ypencil =  my_rhs_ypencil + fl%rre * dmdx_cpc_ypencil * acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! X-pencil : Y-mom diffusion term (y-v5/7), d(mu^y)/dx * d(qy^x))/dx at (i, j', k)
!----------------------------------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2C_3D(fl%qy, acpc, dm, dm%ibcx(:, 2), dm%fbcx(:, 2) )
      fl%my_rhs =  fl%my_rhs + fl%rre * dmdx_cpc * acpc
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v6/7), d(mu^y)/dz * d(qz^z)/dy at (i, j', k)
!----------------------------------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2P_3D(qz_ccc_ypencil, acpc_ypencil, dm, dm%ibcy(:, 3), dm%fbcy(:, 3) )
      my_rhs_ypencil =  my_rhs_ypencil + fl%rre * dmdz_cpc_ypencil * acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Y-mom diffusion term (y-v7/7), d(mu^y)/dz * d(qy)/dz at (i, j', k)
!----------------------------------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2C_3D(qy_zpencil, acpc_zpencil, dm, dm%ibcz(:, 2), dm%fbcz(:, 2) )
      my_rhs_zpencil =  my_rhs_zpencil + fl%rre * dmdz_cpc_zpencil * acpc_zpencil
    end if
!----------------------------------------------------------------------------------------------------------
! y-mom: convert all terms to x-pencil
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_x (my_rhs_ypencil, acpc, dm%dcpc)
    fl%my_rhs =  fl%my_rhs + acpc
    call transpose_z_to_y (my_rhs_zpencil, acpc_ypencil, dm%dcpc)
    call transpose_y_to_x (acpc_ypencil, acpc, dm%dcpc)
    fl%my_rhs =  fl%my_rhs + acpc

    call transpose_y_to_x (my_rhs_pfc_ypencil,  my_rhs_pfc,  dm%dcpc)

!write(*,*) nrank,  'test-6'
!==========================================================================================================
! the RHS of Z momentum equation
!==========================================================================================================
    fl%mz_rhs = ZERO
    mz_rhs_ypencil = ZERO
    mz_rhs_zpencil = ZERO

    mz_rhs_pfc  = ZERO
    mz_rhs_pfc_ypencil = ZERO
    mz_rhs_pfc_zpencil = ZERO

    accp = ZERO
    accp_ypencil = ZERO
    accp_zpencil = ZERO
!----------------------------------------------------------------------------------------------------------
! X-pencil : Z-mom convection term (z-c1/3), d(gx^z * qz^x)/dx at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_x_1st_derivative_P2C_3D( -qx_pcp * qz_pcp, accp, dm, dm%ibcx(:, 3) )
    end if
    if ( dm%is_thermo) then
      call Get_x_1st_derivative_P2C_3D( -gx_pcp * qz_pcp, accp, dm, dm%ibcx(:, 3) )
    end if
    fl%mz_rhs = fl%mz_rhs + accp
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Z-mom convection term (z-c2/3), d(gy^z * qz^y)/dy at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_y_1st_derivative_P2C_3D( -qy_cpp_ypencil * qz_cpp_ypencil, accp_ypencil, dm, dm%ibcy(:, 3) )
    end if
    if ( dm%is_thermo) then
      call Get_y_1st_derivative_P2C_3D( -gy_cpp_ypencil * qz_cpp_ypencil, accp_ypencil, dm, dm%ibcy(:, 3) )
    end if
    mz_rhs_ypencil = mz_rhs_ypencil + accp_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Z-mom convection term (y-c3/3), d(gz * qz)/dz at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      do m = 1, 2
        fbc(m) = dm%fbcz(m, 3) * dm%fbcz(m, 3)
      end do
      call Get_z_1st_derivative_C2P_3D(-qz_ccc_zpencil * qz_ccc_zpencil, accp_zpencil, dm, dm%ibcz(:, 3), fbc(:) )
    end if
    if ( dm%is_thermo) then
      do m = 1, 2
        fbc(m) = dm%fbcz(m, 3) * dm%fbcz(m, 3) * dm%fbc_dend(m, 3)
      end do
      call Get_z_1st_derivative_C2P_3D(-gz_ccc_zpencil * qz_ccc_zpencil, accp_zpencil, dm, dm%ibcz(:, 3), fbc(:) )
    end if
    mz_rhs_zpencil = mz_rhs_zpencil + accp_zpencil

!----------------------------------------------------------------------------------------------------------
! z-pencil : pressure gradient in z direction, d(sigma_1 p)
!----------------------------------------------------------------------------------------------------------
    call Get_z_1st_derivative_C2P_3D(pres_zpencil, accp_zpencil, dm, dm%ibcz(:, 4), dm%fbcz(:, 4) )
    mz_rhs_pfc_zpencil =  mz_rhs_pfc_zpencil - accp_zpencil
!----------------------------------------------------------------------------------------------------------
! X-pencil : Z-mom diffusion term (z-v1-1/7), \mu * L11(uz) at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    !call Get_x_2nd_derivative_C2C_3D(fl%qz, accp, dm, dm%ibcx(:, 3) )
    call Get_x_1st_derivative_C2P_3D(fl%qz, apcp, dm, dm%ibcx(:, 3),  dm%fbcx(:, 3))
    call Get_x_1st_derivative_P2C_3D(apcp,  accp, dm, dm%ibcx(:, 3) )
    if ( .not. dm%is_thermo) fl%mz_rhs = fl%mz_rhs +         fl%rre * accp
    if ( dm%is_thermo) fl%mz_rhs = fl%mz_rhs + m_ccp * fl%rre * accp
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Z-mom diffusion term (z-v1-2/1), \mu * L22(uz) at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    !call Get_y_2nd_derivative_C2C_3D( qz_ypencil, accp_ypencil, dm, dm%ibcy(:, 3), dm%fbcy(:, 3))
    call Get_y_1st_derivative_C2P_3D( qz_ypencil,   acpp_ypencil, dm, dm%ibcy(:, 3), dm%fbcy(:, 3))
    call Get_y_1st_derivative_P2C_3D( acpp_ypencil, accp_ypencil, dm, dm%ibcy(:, 3))
    if ( .not. dm%is_thermo) mz_rhs_ypencil = mz_rhs_ypencil +                 fl%rre * accp_ypencil
    if ( dm%is_thermo) mz_rhs_ypencil = mz_rhs_ypencil + m_ccp_ypencil * fl%rre * accp_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v1-3/7), \mu * L33(uz) at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    !call Get_z_2nd_derivative_P2P_3D(qz_zpencil, accp_zpencil, dm, dm%ibcz(:, 3))
    call Get_z_1st_derivative_P2C_3D(qz_zpencil,   accc_zpencil, dm, dm%ibcz(:, 3))
    call Get_z_1st_derivative_C2P_3D(accc_zpencil, accp_zpencil, dm, dm%ibcz(:, 3), dm%fbcz(:, 3))
    if ( .not. dm%is_thermo) mz_rhs_zpencil = mz_rhs_zpencil +                 fl%rre * accp_zpencil
    if ( dm%is_thermo) mz_rhs_zpencil = mz_rhs_zpencil + m_ccp_zpencil * fl%rre * accp_zpencil

    if ( dm%is_thermo) then
!----------------------------------------------------------------------------------------------------------
! z-pencil : Z-mom gravity force in z direction
!----------------------------------------------------------------------------------------------------------
      if( fl%igravity == 3 .or. fl%igravity == -3) then
        call Get_z_midp_C2P_3D(dDens_zpencil, accp_zpencil, dm, dm%ibcz(:, 5), dm%fbc_dend(:, 3) )
        mz_rhs_pfc_zpencil =  mz_rhs_pfc_zpencil + fl%fgravity(3) * accp_zpencil
      end if
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v2/7), \mu * 1/3 * d (div)/dz at (i, j, k')
!----------------------------------------------------------------------------------------------------------
      fbc = ZERO
      call Get_z_1st_derivative_C2P_3D(div_zpencil, accp_zpencil, dm, dm%ibcz(:, 3), fbc )
      mz_rhs_zpencil = mz_rhs_zpencil + one_third_rre * m_ccp_zpencil * accp_zpencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v3/7), 2d\mu/dz * (-1/3 * div(u)) +  2d\mu/dz * dw/dz
!----------------------------------------------------------------------------------------------------------
      fbc = ZERO
      call Get_z_midp_C2P_3D(div_zpencil, accp_zpencil, dm,  dm%ibcz(:, 3), fbc )
      mz_rhs_zpencil = mz_rhs_zpencil - two_third_rre * dmdz_ccp_zpencil * accp_zpencil
      call Get_z_1st_derivative_P2P_3D( qz_zpencil,  accp_zpencil, dm, dm%ibcz(:, 3), dm%fbcz(:, 3) )
      mz_rhs_zpencil = mz_rhs_zpencil + two_rre * dmdz_ccp_zpencil * accp_zpencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v4/7), d(mu^z)/dx * d(qx^x)/dz at (i, j, k')
!----------------------------------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2P_3D(qx_ccc_zpencil, accp_zpencil, dm, dm%ibcz(:, 1), dm%fbcz(:, 1) )
      mz_rhs_zpencil = mz_rhs_zpencil + fl%rre * dmdx_ccp_zpencil * accp_zpencil
!----------------------------------------------------------------------------------------------------------
! X-pencil : Z-mom diffusion term (z-v5/7), d(mu^z)/dx * d(qz)/dx at (i, j, k')
!----------------------------------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2C_3D(fl%qz, accp, dm, dm%ibcx(:, 3), dm%fbcx(:, 3) )
      fl%mz_rhs =  fl%mz_rhs + fl%rre * dmdx_ccp * accp
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v6/7), d(mu^z)/dy * d(qy^y)/dz at (i, j, k')
!----------------------------------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2P_3D(qy_ccc_zpencil, accp_zpencil, dm, dm%ibcz(:, 2), dm%fbcz(:, 2) )
      mz_rhs_zpencil = mz_rhs_zpencil + fl%rre * dmdy_ccp_zpencil * accp_zpencil
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Z-mom diffusion term (z-v7/7), d(mu^z)/dy * d(qz)/dy at (i, j, k')
!----------------------------------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2C_3D(qz_ypencil, accp_ypencil, dm, dm%ibcy(:, 3), dm%fbcy(:, 3) )
      mz_rhs_ypencil =  mz_rhs_ypencil + fl%rre * dmdy_ccp_ypencil * accp_ypencil
    end if
!----------------------------------------------------------------------------------------------------------
! z-mom: convert all terms to x-pencil
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_x (mz_rhs_ypencil, accp, dm%dccp)
    fl%mz_rhs =  fl%mz_rhs + accp
    call transpose_z_to_y (mz_rhs_zpencil, accp_ypencil, dm%dccp)
    call transpose_y_to_x (accp_ypencil, accp, dm%dccp)
    fl%mz_rhs =  fl%mz_rhs + accp

    call transpose_z_to_y (mz_rhs_pfc_zpencil, accp_ypencil, dm%dccp)
    call transpose_y_to_x (accp_ypencil, mz_rhs_pfc, dm%dccp)
    
!==========================================================================================================
! x-pencil : to build up rhs in total, in all directions
!==========================================================================================================
    !write(*,*) nrank,  'test-7'
!----------------------------------------------------------------------------------------------------------
! x-pencil : x-momentum
!----------------------------------------------------------------------------------------------------------
    call Calculate_momentum_fractional_step(fl%mx_rhs0, fl%mx_rhs, mx_rhs_pfc, dm%dpcc, dm, isub)  
!----------------------------------------------------------------------------------------------------------
! x-pencil : y-momentum
!----------------------------------------------------------------------------------------------------------
    call Calculate_momentum_fractional_step(fl%my_rhs0, fl%my_rhs, my_rhs_pfc, dm%dcpc, dm, isub)
!----------------------------------------------------------------------------------------------------------
! x-pencil : z-momentum
!----------------------------------------------------------------------------------------------------------
    call Calculate_momentum_fractional_step(fl%mz_rhs0, fl%mz_rhs, mz_rhs_pfc, dm%dccp, dm, isub)
!==========================================================================================================
! x-pencil : flow drive terms (source terms) in periodic Streamwise flow
!==========================================================================================================
    if (fl%idriven == IDRVF_X_MASSFLUX) then
      call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), dm%fbcy(:, 1), dm, dm%dpcc, fl%mx_rhs, rhsx_bulk, "mx_rhs")
      fl%mx_rhs(:, :, :) = fl%mx_rhs(:, :, :) - rhsx_bulk
    else if (fl%idriven == IDRVF_X_Cf) then
      rhsx_bulk = - HALF * fl%drvfc * dm%tAlpha(isub) * dm%dt
      fl%mx_rhs(:, :, :) = fl%mx_rhs(:, :, :) - rhsx_bulk
    else if (fl%idriven == IDRVF_Z_MASSFLUX) then
      call Get_volumetric_average_3d(.false., dm%ibcy(:, 3), dm%fbcy(:, 3), dm, dm%dccp, fl%mz_rhs, rhsz_bulk, "mz_rhs")
      fl%mz_rhs(:, :, :) = fl%mz_rhs(:, :, :) - rhsz_bulk
    else if (fl%idriven == IDRVF_Z_Cf) then
      rhsz_bulk = - HALF * fl%drvfc * dm%tAlpha(isub) * dm%dt
      fl%mz_rhs(:, :, :) = fl%mz_rhs(:, :, :) - rhsz_bulk
    else
    end if
 
    return
  end subroutine Compute_momentum_rhs

!==========================================================================================================
!==========================================================================================================
  subroutine Correct_massflux(ux, uy, uz, phi_ccc, dm, isub)
    use udf_type_mod
    use input_general_mod
    use operations
    use parameters_constant_mod
    implicit none

    type(t_domain), intent(in) :: dm
    integer,        intent(in) :: isub
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ), intent(inout) :: ux
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ), intent(inout) :: uy
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ), intent(inout) :: uz
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ), intent(in ) :: phi_ccc

    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: dphidx_pcc
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: dphidy_cpc
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: dphidz_ccp

    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: phi_ccc_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: dphidy_cpc_ypencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: dphidz_ccp_ypencil
    
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: pphi_ccc_zpencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: dphidz_ccp_zpencil
    
#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_mid_msg("Updating the velocity/mass flux ...")
#endif
!----------------------------------------------------------------------------------------------------------
!   x-pencil, ux = ux - dt * alpha * d(phi_ccc)/dx
!_______________________________________________________________________________
    dphidx_pcc = ZERO
    call Get_x_1st_derivative_C2P_3D(phi_ccc,  dphidx_pcc, dm, dm%ibcx(:, 4), dm%fbcx(:, 4) )
    ux = ux - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidx_pcc
!----------------------------------------------------------------------------------------------------------
!   y-pencil, uy = uy - dt * alpha * d(phi_ccc)/dy
!_______________________________________________________________________________
    phi_ccc_ypencil = ZERO
    dphidy_cpc_ypencil = ZERO
    dphidy_cpc = ZERO
    call transpose_x_to_y (phi_ccc, phi_ccc_ypencil, dm%dccc)
    call Get_y_1st_derivative_C2P_3D(phi_ccc_ypencil, dphidy_cpc_ypencil, dm, dm%ibcy(:, 4), dm%fbcy(:, 4) )
    call transpose_y_to_x (dphidy_cpc_ypencil, dphidy_cpc, dm%dcpc)
    uy = uy - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidy_cpc
!----------------------------------------------------------------------------------------------------------
!   z-pencil, uz = uz - dt * alpha * d(phi_ccc)/dz
!_______________________________________________________________________________
    pphi_ccc_zpencil =  ZERO
    dphidz_ccp_zpencil = ZERO
    dphidz_ccp_ypencil = ZERO
    dphidz_ccp = ZERO
    call transpose_y_to_z (phi_ccc_ypencil, pphi_ccc_zpencil, dm%dccc)
    call Get_z_1st_derivative_C2P_3D(pphi_ccc_zpencil, dphidz_ccp_zpencil, dm, dm%ibcz(:, 4), dm%fbcz(:, 4) )
    call transpose_z_to_y (dphidz_ccp_zpencil, dphidz_ccp_ypencil, dm%dccp)
    call transpose_y_to_x (dphidz_ccp_ypencil, dphidz_ccp,         dm%dccp)
    uz = uz - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidz_ccp

    return
  end subroutine Correct_massflux

!==========================================================================================================
  subroutine solve_poisson(fl, dm, isub)
    use udf_type_mod
    use parameters_constant_mod
    use decomp_2d_poisson
    use decomp_extended_mod
    use continuity_eq_mod

    implicit none
    type(t_domain), intent( in    ) :: dm
    type(t_flow),   intent( inout ) :: fl                  
    integer,        intent( in    ) :: isub

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: div
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: rhs_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: rhs_zpencil

  

    real(WP), dimension( dm%dccc%zst(1) : dm%dccc%zen(1), &
                         dm%dccc%zst(2) : dm%dccc%zen(2), &
                         dm%dccc%zst(3) : dm%dccc%zen(3) ) :: rhs_zpencil_ggg
    integer :: i, j, k, jj, ii


#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_start_msg("Calculating the RHS of Poisson Equation ...")
#endif

!==========================================================================================================
! RHS of Poisson Eq.
!==========================================================================================================
    fl%pcor = ZERO
!----------------------------------------------------------------------------------------------------------
! $d\rho / dt$ at cell centre
!----------------------------------------------------------------------------------------------------------
    if (dm%is_thermo) then
      call Calculate_drhodt(dm, fl%dDens, fl%dDensm1, fl%dDensm2, fl%pcor)
    end if
!----------------------------------------------------------------------------------------------------------
! $d(\rho u_i)) / dx_i $ at cell centre
!----------------------------------------------------------------------------------------------------------
    div  = ZERO
    if (dm%is_thermo) then
      call Get_divergence_vel(fl%gx, fl%gy, fl%gz, div, dm)
    else
      call Get_divergence_vel(fl%qx, fl%qy, fl%qz, div, dm)
    end if

    fl%pcor = fl%pcor + div
    fl%pcor = fl%pcor / (dm%tAlpha(isub) * dm%sigma2p * dm%dt)

!==========================================================================================================
!   convert RHS from xpencil gll to zpencil ggg
!==========================================================================================================
    call transpose_x_to_y (fl%pcor    , rhs_ypencil, dm%dccc)
    call transpose_y_to_z (rhs_ypencil, rhs_zpencil, dm%dccc)
    call zpencil_index_llg2ggg(rhs_zpencil, rhs_zpencil_ggg, dm%dccc)

!==========================================================================================================
!   solve Poisson
!==========================================================================================================
#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_mid_msg("Solving the Poisson Equation ...")
#endif
    call poisson(rhs_zpencil_ggg)
!==========================================================================================================
!   convert back RHS from zpencil ggg to xpencil gll
!==========================================================================================================
    call zpencil_index_ggg2llg(rhs_zpencil_ggg, rhs_zpencil, dm%dccc)
    call transpose_z_to_y (rhs_zpencil, rhs_ypencil, dm%dccc)
    call transpose_y_to_x (rhs_ypencil, fl%pcor,     dm%dccc)

    return
  end subroutine
!==========================================================================================================
  subroutine solve_poisson_x2z(fl, dm, isub)
    use udf_type_mod
    use parameters_constant_mod
    use decomp_2d_poisson
    use decomp_extended_mod
    use continuity_eq_mod

    implicit none
    type(t_domain), intent( in    ) :: dm
    type(t_flow),   intent( inout ) :: fl                  
    integer,        intent( in    ) :: isub

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: rhs
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: rhs_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: rhs_zpencil

  

    real(WP), dimension( dm%dccc%zst(1) : dm%dccc%zen(1), &
                         dm%dccc%zst(2) : dm%dccc%zen(2), &
                         dm%dccc%zst(3) : dm%dccc%zen(3) ) :: rhs_zpencil_ggg
    integer :: i, j, k, jj, ii

!==========================================================================================================
! RHS of Poisson Eq.
!==========================================================================================================
    fl%pcor_zpencil_ggg = ZERO
!----------------------------------------------------------------------------------------------------------
! $d\rho / dt$ at cell centre
!----------------------------------------------------------------------------------------------------------
    if (dm%is_thermo) then
      rhs = ZERO
      rhs_ypencil = ZERO
      rhs_zpencil = ZERO
      rhs_zpencil_ggg = ZERO
      call Calculate_drhodt(dm, fl%dDens, fl%dDensm1, fl%dDensm2, rhs)
      call transpose_x_to_y(rhs,         rhs_ypencil)
      call transpose_y_to_z(rhs_ypencil, rhs_zpencil)
      call zpencil_index_llg2ggg(rhs_zpencil, rhs_zpencil_ggg, dm%dccc)
      fl%pcor_zpencil_ggg = fl%pcor_zpencil_ggg + rhs_zpencil_ggg
    end if
!----------------------------------------------------------------------------------------------------------
! $d(\rho u_i)) / dx_i $ at cell centre
!----------------------------------------------------------------------------------------------------------
    rhs_zpencil_ggg  = ZERO
    if (dm%is_thermo) then
      call Get_divergence_vel_x2z(fl%gx, fl%gy, fl%gz, rhs_zpencil_ggg, dm)
    else
      call Get_divergence_vel_x2z(fl%qx, fl%qy, fl%qz, rhs_zpencil_ggg, dm)
    end if
    fl%pcor_zpencil_ggg = fl%pcor_zpencil_ggg + rhs_zpencil_ggg
    fl%pcor_zpencil_ggg = fl%pcor_zpencil_ggg / (dm%tAlpha(isub) * dm%sigma2p * dm%dt)
!==========================================================================================================
!   solve Poisson
!==========================================================================================================
    call poisson(fl%pcor_zpencil_ggg)
!==========================================================================================================
!   convert back RHS from zpencil ggg to xpencil gll
!==========================================================================================================
    call zpencil_index_ggg2llg(fl%pcor_zpencil_ggg, rhs_zpencil, dm%dccc)
    call transpose_z_to_y (rhs_zpencil, rhs_ypencil, dm%dccc)
    call transpose_y_to_x (rhs_ypencil, fl%pcor,     dm%dccc)

    return
  end subroutine

!==========================================================================================================
!> \brief To update the provisional u or rho u.
!>
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                         !
!______________________________________________________________________________!
!> \param[inout]  fl            flow field
!> \param[inout]  dm            domain
!> \param[in]     isub         RK sub-iteration
!==========================================================================================================
  subroutine Solve_momentum_eq(fl, dm, isub)
    use udf_type_mod,      only : t_flow, t_domain
    use typeconvert_mod
    use continuity_eq_mod
    use boundary_conditions_mod
    use parameters_constant_mod
    use mpi_mod
    use solver_tools_mod
#ifdef DEBUG_STEPS
    use io_visulisation_mod
    use typeconvert_mod
#endif
    implicit none

    type(t_flow), intent(inout) :: fl
    type(t_domain),  intent(in) :: dm
    integer,      intent(in) :: isub

!----------------------------------------------------------------------------------------------------------
! to calculate the rhs of the momenturn equation in stepping method
!_______________________________________________________________________________ 
    call Compute_momentum_rhs(fl, dm, isub)
!----------------------------------------------------------------------------------------------------------
! to update intermediate (\hat{q}) or (\hat{g})
!_______________________________________________________________________________
#ifdef DEBUG_STEPS
  !call write_snapshot_any3darray(fl%mx_rhs, 'mx_rhs_RK'//trim(int2str(isub)), 'debug', dm%dpcc, dm, fl%iteration)
  !call write_snapshot_any3darray(fl%my_rhs, 'my_rhs_RK'//trim(int2str(isub)), 'debug', dm%dcpc, dm, fl%iteration)
  !call write_snapshot_any3darray(fl%mz_rhs, 'mz_rhs_RK'//trim(int2str(isub)), 'debug', dm%dccp, dm, fl%iteration)
#endif
!----------------------------------------------------------------------------------------------------------
! to update velocity to get intermediate data
!_______________________________________________________________________________ 
    if ( .not. dm%is_thermo) then 
      fl%qx = fl%qx + fl%mx_rhs
      fl%qy = fl%qy + fl%my_rhs
      fl%qz = fl%qz + fl%mz_rhs
    else if ( dm%is_thermo) then 
      fl%gx = fl%gx + fl%mx_rhs
      fl%gy = fl%gy + fl%my_rhs
      fl%gz = fl%gz + fl%mz_rhs
    else
      call Print_error_msg("Error in velocity updating")
    end if

    !in order for a high order spacial accuracy
    ! to use Alternating direction implicit method
    ! ref: Cui2013: Convergence analysis of high-order compact 
    ! alternating direction implicit schemes for the two-dimensional 
    ! time fractional equation

!----------------------------------------------------------------------------------------------------------
! to update b.c. values
!----------------------------------------------------------------------------------------------------------
    call Apply_BC_velocity (dm, fl%qx, fl%qy, fl%qz)
    if(dm%is_thermo) call Apply_BC_velocity (dm, fl%gx, fl%gy, fl%gz)
!----------------------------------------------------------------------------------------------------------
! to solve Poisson equation
!----------------------------------------------------------------------------------------------------------
    !if(nrank == 0) call Print_debug_mid_msg("  Solving Poisson Equation ...") 
    !call solve_poisson_x2z(fl, dm, isub) ! mpi=4 error, to check
    call solve_poisson(fl, dm, isub) ! test show above two methods gave the same results. 
!----------------------------------------------------------------------------------------------------------
! to update velocity/massflux correction
!----------------------------------------------------------------------------------------------------------
    !if(nrank == 0) call Print_debug_mid_msg("  Updating velocity/mass flux ...")
    if ( .not. dm%is_thermo) then 
      call Correct_massflux(fl%qx, fl%qy, fl%qz, fl%pcor, dm, isub)
    else
      call Correct_massflux(fl%gx, fl%gy, fl%gz, fl%pcor, dm, isub)
    end if
!----------------------------------------------------------------------------------------------------------
! to update pressure
!----------------------------------------------------------------------------------------------------------
#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_mid_msg("Correcting the pressure term ...")
#endif
    fl%pres(:, :, :) = fl%pres(:, :, :) + fl%pcor(:, :, :)
!----------------------------------------------------------------------------------------------------------
! to update b.c. values
!----------------------------------------------------------------------------------------------------------
    call Apply_BC_velocity (dm, fl%qx, fl%qy, fl%qz)
    if ( dm%is_thermo) call Apply_BC_velocity (dm, fl%gx, fl%gy, fl%gz)
    
#ifdef DEBUG_STEPS
    !call write_snapshot_any3darray(fl%qx, 'qx_RK'//trim(int2str(isub)), 'debug', dm%dpcc, dm, fl%iteration)
    !call write_snapshot_any3darray(fl%qy, 'qy_RK'//trim(int2str(isub)), 'debug', dm%dcpc, dm, fl%iteration)
    !call write_snapshot_any3darray(fl%qz, 'qz_RK'//trim(int2str(isub)), 'debug', dm%dccp, dm, fl%iteration)
#endif

    return
  end subroutine Solve_momentum_eq

end module eq_momentum_mod
