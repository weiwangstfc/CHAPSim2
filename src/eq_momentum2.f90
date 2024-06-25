module eq_momentum_mod
  use operations
  use precision_mod
  use decomp_2d
  implicit none

  private :: Calculate_momentum_fractional_step
  private :: Compute_momentum_rhs
  private :: Correct_massflux
  private :: solve_poisson
  !private :: solve_poisson_x2z
  
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
    use typeconvert_mod
    use boundary_conditions_mod
    implicit none

    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(in )   :: dm
    integer,     intent(in ) :: isub

!----------------------------------------------------------------------------------------------------------
! result variables
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: mx_rhs_ypencil
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: mx_rhs_zpencil
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: mx_rhs_pfc_xpencil

    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: my_rhs_ypencil ! 
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: my_rhs_zpencil ! 
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: my_rhs_pfc_xpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: my_rhs_pfc_ypencil
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: my_rhs_pfc_zpencil

    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: mz_rhs_ypencil !
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: mz_rhs_zpencil ! 
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: mz_rhs_pfc_xpencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: mz_rhs_pfc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: mz_rhs_pfc_zpencil
!----------------------------------------------------------------------------------------------------------
! intermediate variables : common
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: pres_ypencil ! p
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: pres_zpencil ! p

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: div_ccc_xpencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: div_ccc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: div_ccc_zpencil

    real(WP), dimension(              4, dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: div_4cc_xpencil ! for dirichlet bc.
    real(WP), dimension( dm%dccc%ysz(1),              4, dm%dccc%ysz(3) ) :: div_c4c_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2),              4 ) :: div_cc4_zpencil
    
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: qxix_ccc_xpencil ! qx^x, x-c1, common
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: qxiy_ppc_xpencil ! qx^y, y-c1, no-thermal
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: qxiy_ppc_ypencil ! qx^y, x-c2, common
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: qxiz_pcp_xpencil ! qx^z, z-c1, no-thermal
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: qxiz_pcp_zpencil ! qx^z, x-c3, common
    
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: qxdx_ccc_xpencil ! d(qx)/dx, x-v1, common
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: qxdy_ppc_xpencil ! d(qx)/dy, y-v1, common
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: qxdy_ppc_ypencil ! d(qx)/dy, x-v2, common
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: qxdz_pcp_xpencil ! d(qx)/dz, z-v1, common
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: qxdz_pcp_zpencil ! d(qx)/dz, x-v3, common
    real(WP), dimension(              4, dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: qxdx_4cc_xpencil !
    real(WP), dimension( dm%dccc%ysz(1),              4, dm%dccc%ysz(3) ) :: qxdx_c4c_ypencil ! 
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2),              4 ) :: qxdx_cc4_zpencil ! 
    
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: qyix_ppc_xpencil ! qy^x, y-c1, common
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: qyix_ppc_ypencil ! qy^x, x-c2, no-thermal
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qyiy_ccc_ypencil ! qy^y, y-c2, no-thermal OR no-cyl
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qyiz_cpp_ypencil ! qy^z, z-c2, no-thermal
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qyiz_cpp_zpencil ! qy^z, y-c3, no-cyl
    
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: qydx_ppc_xpencil ! d(qy)/dx, y-v1, common
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: qydx_ppc_ypencil ! d(qy)/dx, x-v2, common
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qydy_ccc_ypencil ! d(qy)/dy, y-v2, no-cly
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qydz_cpp_ypencil ! d(qy)/dz, z-v2, z-v4, no-cly
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qydz_cpp_zpencil ! d(qy)/dz, y-v3, no-cly
    real(WP), dimension(              4, dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: qydy_4cc_xpencil !
    real(WP), dimension( dm%dccc%ysz(1),              4, dm%dccc%ysz(3) ) :: qydy_c4c_ypencil ! 
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2),              4 ) :: qydy_cc4_zpencil ! 

    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: qzix_pcp_xpencil ! qz^x, z-c1, common
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: qzix_pcp_zpencil ! qz^x, x-c3, no-thermal
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qziy_cpp_ypencil ! qz^y, z-c2, no-cly
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qziy_cpp_zpencil ! qz^y, y-c3, no-thermal and no-cly
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: qziz_ccc_zpencil ! qz^z, z-c3, common
    
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: qzdx_pcp_xpencil ! d(qz)/dx, z-v1, common
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: qzdx_pcp_zpencil ! d(qz)/dx, x-v3, common
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qzdy_cpp_ypencil ! d(qz)/dy, z-v2, no-cly
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qzdy_cpp_zpencil ! d(qz)/dy, y-v3, no-cly
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qzdz_ccc_ypencil ! d(qz)/dz, y-v4, cly-only
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: qzdz_ccc_zpencil ! d(qz)/dz, z-v3, common
    real(WP), dimension(              4, dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: qzdz_4cc_xpencil !
    real(WP), dimension( dm%dccc%ysz(1),              4, dm%dccc%ysz(3) ) :: qzdz_c4c_ypencil ! 
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2),              4 ) :: qzdz_cc4_zpencil ! 
!----------------------------------------------------------------------------------------------------------
!   intermediate variables : if dm%is_thermo = true
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: gxix_ccc_xpencil  ! gx^x, x-c1, thermal-only
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: gxiy_ppc_xpencil  ! gx^y, y-c1, thermal-only
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: gxiz_pcp_xpencil  ! gx^z, z-c1, thermal-only
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: gyix_ppc_ypencil  ! gy^x, x-c2, thermal-only
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: gyiy_ccc_ypencil  ! gy^y, y-c2, thermal-only
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: gyiz_cpp_ypencil  ! gy^z, z-c2, thermal-only
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: gzix_pcp_zpencil  ! gz^x, x-c3, thermal-only
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: gziy_cpp_zpencil  ! gz^y, y-c3, thermal-only, and no-cly
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: gziz_ccc_zpencil  ! gz^z, z-c3, thermal-only
    
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: mu_ccc_xpencil ! x-v1, thermal only
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: mu_ccc_ypencil ! y-v2, thermal only
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: mu_ccc_zpencil ! z-v3, thermal only
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: muixy_ppc_xpencil ! y-v1, thermal only
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: muixy_ppc_ypencil ! x-v2, thermal only
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: muixz_pcp_xpencil ! z-v1, thermal only
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: muixz_pcp_zpencil ! x-v3, thermal only
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: muiyz_cpp_ypencil ! z-v2, z-v4, thermal only
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: muiyz_cpp_zpencil ! y-v3, thermal only
    real(WP), dimension(              4, dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: mu_4cc_xpencil !
    real(WP), dimension( dm%dccc%ysz(1),              4, dm%dccc%ysz(3) ) :: mu_c4c_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2),              4 ) :: mu_cc4_zpencil
!----------------------------------------------------------------------------------------------------------
!   intermediate variables : if fl%igravity == i
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: dDens_ypencil  ! d 
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: dDens_zpencil  ! d 
!----------------------------------------------------------------------------------------------------------
!   intermediate variables : if dm%icoordinate == ICYLINDRICAL
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qyriy_ccc_ypencil ! (qy/r)^y, y-c2, y-v4, z-v3, cly-only
    real(WP), dimension( dm%dccc%ysz(1),              4, dm%dccc%ysz(3) ) :: qyriy_c4c_ypencil ! (qy/r)^y, y-c2, y-v4, z-v3, cly-only, BC
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: qyriy_ccc_zpencil !
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qyriz_cpp_ypencil ! (qy/r)^z, z-c4, no-thermal
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qyriz_cpp_zpencil ! (qy/r)^z, y-c3
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: qzriz_ccc_zpencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qyrdy_ccc_ypencil ! d(qy/r)/dy, y-v2
    real(WP), dimension( dm%dccc%ysz(1),              4, dm%dccc%ysz(3) ) :: qyrdy_c4c_ypencil ! 
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qyrdz_cpp_zpencil ! d(qy/r)/dz, y-v3
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qyrdz_cpp_ypencil ! d(qy/r)/dz, z-v2
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qzriz_ccc_ypencil ! (qz/r)^z, y-c4
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qzriy_cpp_ypencil ! (qz/r)^r, y-c2, z-c4, z-v2-f, z-v4
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qzriy_cpp_zpencil ! (qz/r)^y, y-c3, y-v3-f
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qzrdy_cpp_ypencil ! d(qz/r)/dy, z-v2, z-v4
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qzrdy_cpp_zpencil ! d(qz/r)/dy, y-v3
    
!----------------------------------------------------------------------------------------------------------
! intermediate variables : if dm%is_thermo = true && dm%icoordinate == ICYLINDRICAL
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: gzriy_cpp_zpencil ! (gz/r)^y, y-c3
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: gzriz_ccc_ypencil ! (gz/r)^z, y-c4
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: gyriz_cpp_ypencil ! (gy/r)^z, z-c4
!----------------------------------------------------------------------------------------------------------
! temporary variables
!----------------------------------------------------------------------------------------------------------
#ifdef DEBUG_STEPS  
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: apcc_test
#endif 
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: apcc_xpencil
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc_xpencil
    real(WP), dimension(              4, dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: a4cc_xpencil
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: appc_xpencil
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: appc_xpencil1
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: apcp_xpencil
    
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: appc_ypencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil1
    real(WP), dimension( dm%dccc%ysz(1),             4,  dm%dccc%ysz(3) ) :: ac4c_ypencil
    real(WP), dimension( dm%dccc%ysz(1),             4,  dm%dccc%ysz(3) ) :: ac4c_ypencil1
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: acpp_ypencil
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: acpp_ypencil1
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil1
    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: apcc_ypencil
    real(WP), dimension( dm%dpcp%ysz(1), dm%dpcp%ysz(2), dm%dpcp%ysz(3) ) :: apcp_ypencil
    
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: apcp_zpencil
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: acpp_zpencil
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: acpp_zpencil1
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2),              4 ) :: acc4_zpencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil1
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2),              4 ) :: acc4_zpencil1
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil1
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: acpc_zpencil
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: apcc_zpencil

    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: accp_xpencil
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: acpc_xpencil
!----------------------------------------------------------------------------------------------------------
! bc variables
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( 4, dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: fbcx_4cc
    real(WP), dimension( 4, dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: fbcx_4pc
    real(WP), dimension( 4, dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: fbcx_4cp
    
    real(WP), dimension( dm%dpcc%ysz(1), 4, dm%dpcc%ysz(3) ) :: fbcy_p4c
    real(WP), dimension( dm%dccc%ysz(1), 4, dm%dccc%ysz(3) ) :: fbcy_c4c
    real(WP), dimension( dm%dccp%ysz(1), 4, dm%dccp%ysz(3) ) :: fbcy_c4p
    
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), 4 ) :: fbcz_pc4
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), 4 ) :: fbcz_cp4
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), 4 ) :: fbcz_cc4
!----------------------------------------------------------------------------------------------------------
! others
!----------------------------------------------------------------------------------------------------------
    integer  :: i
    real(WP) :: rhsx_bulk, rhsz_bulk
    integer :: ibcy(2), ibcz(2)

!==========================================================================================================
! preparation of constants to be used.
!==========================================================================================================
#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_start_msg("Compute_momentum_rhs at isub = "//trim(int2str(isub)))
#endif
!==========================================================================================================
! preparation of intermediate variables to be used - common
! Note: 
!      Due to bc treatment, please do interp first, then derivative to minumize numerical 
!      errors for Dirichlet B.C.
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
!    p --> p_ypencil --> p_zpencil 
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y (fl%pres,      pres_ypencil, dm%dccc) 
    call transpose_y_to_z (pres_ypencil, pres_zpencil, dm%dccc)
!----------------------------------------------------------------------------------------------------------
!   d --> d_ypencil --> d_zpencil
!----------------------------------------------------------------------------------------------------------    
    if(dm%is_thermo .and. fl%igravity /= 0) then
      call transpose_x_to_y (fl%dDens,      dDens_ypencil, dm%dccc)
      call transpose_y_to_z (dDens_ypencil, dDens_zpencil, dm%dccc)
    end if
!----------------------------------------------------------------------------------------------------------
!    qx 
!    | -[ipx]-> qxix_ccc_xpencil(common)
!    | -[1dx]-> qxdx_ccc_xpencil(common)
!    | -[x2y]-> qx_ypencil(temp) 
!               | -[ipy]-> qxiy_ppc_ypencil(common) -[y2x]-> qxiy_ppc_xpencil(no thermal)
!               | -[1dy]-> qxdy_ppc_ypencil(common) -[y2x]-> qxdy_ppc_xpencil(common)
!               | -[y2z]-> qx_zpencil(temp) 
!                          | -[ipz]-> qxiz_pcp_zpencil(common) -[z2y]-> qxiz_pcp_ypencil(temp) -[y2x]-> qxiz_pcp_xpencil(no-thermal)
!                          | -[1dz]-> qxdz_pcp_zpencil(common) -[z2y]-> qxdz_pcp_ypencil(temp) -[y2x]-> qxdz_pcp_xpencil(common)
!----------------------------------------------------------------------------------------------------------
    call Get_x_midp_P2C_3D(fl%qx, qxix_ccc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_qx)
    call Get_x_1st_derivative_P2C_3D(fl%qx, qxdx_ccc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_qx)
    
    call transpose_x_to_y (fl%qx, apcc_ypencil, dm%dpcc) !apcc_ypencil = qx_ypencil
    call Get_y_1st_derivative_C2P_3D(apcc_ypencil, qxdy_ppc_ypencil, dm, dm%iAccuracy, dm%ibcy_qx, dm%fbcy_qx)
    call transpose_y_to_x(qxdy_ppc_ypencil, qxdy_ppc_xpencil, dm%dppc)
    call Get_y_midp_C2P_3D(apcc_ypencil, qxiy_ppc_ypencil, dm, dm%iAccuracy, dm%ibcy_qx, dm%fbcy_qx) 
    if(.not. dm%is_thermo) then
      call transpose_y_to_x (qxiy_ppc_ypencil, qxiy_ppc_xpencil, dm%dppc)
    end if

    call transpose_y_to_z (apcc_ypencil, apcc_zpencil, dm%dpcc)!apcc_zpencil = qx_zpencil
    call Get_z_midp_C2P_3D(apcc_zpencil, qxiz_pcp_zpencil, dm, dm%iAccuracy, dm%ibcz_qx, dm%fbcz_qx)
    if(.not. dm%is_thermo) then
      call transpose_z_to_y (qxiz_pcp_zpencil, apcp_ypencil, dm%dpcc)! apcp_ypencil = qxiz_pcp_ypencil
      call transpose_y_to_x (apcp_ypencil, qxiz_pcp_xpencil, dm%dpcp)
    end if

    call Get_z_1st_derivative_C2P_3D(apcc_zpencil, qxdz_pcp_zpencil, dm, dm%iAccuracy, dm%ibcz_qx, dm%fbcz_qx)
    call transpose_z_to_y(qxdz_pcp_zpencil, apcp_ypencil,     dm%dpcp) !qxdz_pcp_ypencil
    call transpose_y_to_x(apcp_ypencil,     qxdz_pcp_xpencil, dm%dpcp)

    if(is_fbcx_velo_required) then
      fbcx_qx_4cc(1:2, :, :) = dm%fbcx_qx(1:2, :, :)
      call Get_x_1st_derivative_C2P_3D(qxix_ccc_xpencil, qxdx_pcc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_qx)
      fbcx_qxdx_4cc(1, :, :) = qxdx_pcc_xpencil(1, :, :)
      fbcx_qxdx_4cc(2, :, :) = qxdx_pcc_xpencil(dm%dpcc%xsz(1), :, :)
    end if
    if(is_fbcy_velo_required) then
      fbcy_qx_p4c(:, 1, :) = qxiy_ppc_ypencil(:, 1,              :)
      fbcy_qx_p4c(:, 2, :) = qxiy_ppc_ypencil(:, dm%dppc%ysz(2), :)

      qxdx_cpc_ypencil
      call qxiy_ppc_ypencil, 
      call 



    end if
    if(is_fbcz_velo_required) then
      fbcz_qx_pc4(:, :, 1) = qxiz_pcp_zpencil(:, :,              1)
      fbcz_qx_pc4(:, :, 2) = qxiz_pcp_zpencil(:, :, dm%dpcp%zsz(3))
    end if
!----------------------------------------------------------------------------------------------------------
!    qy
!    | -[ipx]-> qyix_ppc_xpencil(common) -[x2y]-> qyix_ppc_ypencil(no thermal) 
!    | -[1dx]-> qydx_ppc_xpencil(common) -[x2y]-> qydx_ppc_ypencil(common)
!    | -[x2y]-> qy_ypencil(temp) 
!               | -[ipy]-> qyiy_ccc_ypencil(no thermal or no cyl)
!               | -[1dy]-> qydy_ccc_ypencil(no-cly) <need BCy> 
!               | -[y2z]-> qy_zpencil(temp) 
!                          | -[ipz]-> qyiz_cpp_zpencil(no-cly) -[z2y]-> qyiz_cpp_ypencil(no-thermal)    
!                          | -[1dz]-> qydz_cpp_zpencil(no-cly) -[z2y]-> qydz_cpp_ypencil(no-cly)
!----------------------------------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3D(fl%qy, qydx_ppc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_qy)
    call transpose_x_to_y(qydx_ppc_xpencil, qydx_ppc_ypencil, dm%dppc)
    call Get_x_midp_C2P_3D(fl%qy, qyix_ppc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_qy)
    if(.not. dm%is_thermo) then
      call transpose_x_to_y (qyix_ppc_xpencil, qyix_ppc_ypencil, dm%dppc)
    end if
    call transpose_x_to_y (fl%qy, acpc_ypencil, dm%dcpc) !acpc_ypencil = qy_ypencil
    if(dm%icoordinate == ICARTESIAN) then
      call Get_y_1st_derivative_P2C_3D(acpc_ypencil, qydy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)
    end if

    if((.not. dm%is_thermo) .or. (dm%icoordinate == ICARTESIAN)) then
      call Get_y_midp_P2C_3D(acpc_ypencil, qyiy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)
    end if

    call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc)!acpc_zpencil = qy_zpencil
    call Get_z_1st_derivative_C2P_3D(acpc_zpencil, qydz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qy)
    call transpose_z_to_y(qydz_cpp_zpencil, qydz_cpp_ypencil, dm%dcpp)

    if(dm%icoordinate == ICARTESIAN) then
      call Get_z_midp_C2P_3D(acpc_zpencil, qyiz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qy)
      if(.not. dm%is_thermo) then
        call transpose_z_to_y(qyiz_cpp_zpencil, qyiz_cpp_ypencil, dm%dcpp)
      end if
    end if


    if(is_fbcx_velo_required) then
      if(dm%is_thermo) then
        call transpose_x_to_y (qyix_ppc_xpencil, qyix_ppc_ypencil, dm%dppc)
      end if
      call Get_y_1st_derivative_P2C_3D(qyix_ppc_ypencil, qydy_pcc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy)
      call transport_y_to_x(qydy_pcc_ypencil, qydy_pcc_xpencil, dm%dpcc)
    end if

    if(is_fbcy_velo_required) then
      fbcy_qy_c4c(:, 1, :) = acpc_ypencil(:, 1,              :)
      fbcy_qy_c4c(:, 2, :) = acpc_ypencil(:, dm%dcpc%ysz(2), :)

      if(dm%icoordinate == ICARTESIAN) then
        call Get_y_1st_derivative_P2P_3D(acpc_ypencil, acpc_ypencil1, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)
        fbcy_qydy_c4c(:, 1, :) = acpc_ypencil1(:, 1,              :)
        fbcy_qydy_c4c(:, 2, :) = acpc_ypencil1(:, dm%dcpc%ysz(2), :)
      end if
    end if
!----------------------------------------------------------------------------------------------------------
!    qz 
!    | -[1dx]-> qzdx_pcp_xpencil(common) -[x2y]-> qzdx_pcp_ypencil(temp) -[y2z]-> qzdx_pcp_zpencil(common)
!    | -[ipx]-> qzix_pcp_xpencil(common) -[x2y]-> qzix_pcp_ypencil(temp) -[y2z]-> qzix_pcp_zpencil(no-thermal)
!    | -[x2y]-> qz_ypencil(temp) 
!               | -[1dy]-> qzdy_cpp_ypencil(no-cly)   -[y2z]-> qzdy_cpp_zpencil(no-cly)
!               | -[ipy]-> qziy_cpp_ypencil(no-cly)   -[y2z]-> qziy_cpp_zpencil(no cly, no thermal)
!               | -[y2z]-> qz_zpencil(temp) 
!                          | -[ipz]-> qziz_ccc_zpencil(common)
!                          | -[1dz]-> qzdz_ccc_zpencil(common) <need BCz> -[y2z]-> qzdz_ccc_ypencil(cly-only) <need BCy>
!----------------------------------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3D(fl%qz, qzdx_pcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_qz)
    call transpose_x_to_y(qzdx_pcp_xpencil, apcp_ypencil, dm%dpcp)
    call transpose_y_to_z(apcp_ypencil, qzdx_pcp_zpencil, dm%dpcp)

    call Get_x_midp_C2P_3D(fl%qz, qzix_pcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_qz)
    if(.not. dm%is_thermo) then
      call transpose_x_to_y(qzix_pcp_xpencil, apcp_ypencil, dm%dpcp) !qzix_pcp_ypencil
      call transpose_y_to_z(apcp_ypencil, qzix_pcp_zpencil, dm%dpcp)
    end if
    call transpose_x_to_y(fl%qz, accp_ypencil, dm%dccp) ! qz_ypencil
    if(dm%icoordinate == ICARTESIAN)then
      call Get_y_1st_derivative_C2P_3D(accp_ypencil, qzdy_cpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qz)
      call transpose_y_to_z(qzdy_cpp_ypencil, qzdy_cpp_zpencil, dm%dcpp)
    end if

    if((.not. dm%is_thermo) .or. (dm%icoordinate == ICARTESIAN)) then
      call Get_y_midp_C2P_3D(accp_ypencil, qziy_cpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qz)
      if(.not. dm%is_thermo) then
        call transpose_y_to_z(qziy_cpp_ypencil, qziy_cpp_zpencil, dm%dcpp)
      end if
    end if

    call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp) ! qz_zpencil
    call Get_z_midp_P2C_3D(accp_zpencil, qziz_ccc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qz)
    call Get_z_1st_derivative_P2C_3D(accp_zpencil, qzdz_ccc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)

    if(dm%icoordinate == ICYLINDRICAL) then
      call transpose_z_to_y(qzdz_ccc_zpencil, qzdz_ccc_ypencil, dm%dccc)
    end if

    if(is_fbcx_velo_required) then
      call Get_z_1st_derivative_P2C_3D(qzix_pcp_zpencil, qzdz_pcc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)
      call transpose_z_to_y(qzdz_pcc_zpencil, qzdz_pcc_ypencil, dm%dpcc )
      call transport_y_to_x(qzdz_pcc_ypencil, qzdz_pcc_xpencil, dm%dpcc)
    end if

    if(is_fbcz_velo_required) then
      fbcz_qz_pc4(:, :, 1) = qzix_pcp_zpencil(:, :,              1)
      fbcz_qz_pc4(:, :, 2) = qzix_pcp_zpencil(:, :, dm%dpcp%zsz(3))
    end if
!==========================================================================================================
! preparation of div
!==========================================================================================================
    div_ccc_xpencil  = ZERO 
    div_ccc_xpencil = div_ccc_xpencil + qxdx_ccc_xpencil ! div_tmp = du/dx
    accc_ypencil = qydy_ccc_ypencil
    if(dm%icoordinate == ICYLINDRICAL) then
      call multiple_cylindrical_rn(accc_ypencil, dm%dccc, dm%rci, 1, IPENCIL(2))
      ! accc_ypencil = qydy_ccc_ypencil * 1/r
    end if
    call transpose_y_to_x (accc_ypencil, accc_xpencil, dm%dccc)   
    div_ccc_xpencil = div_ccc_xpencil + accc_xpencil  ! div_tmp = dux/dx + 1/r * duy/dy

    accc_zpencil = qzdz_ccc_zpencil
    if(dm%icoordinate == ICYLINDRICAL) then
      call multiple_cylindrical_rn(accc_zpencil, dm%dccc, dm%rci, 2, IPENCIL(3))
      ! accc_zpencil = qzdz_ccc_zpencil * 1/r^2
    end if
    call transpose_z_to_y (accc_zpencil, accc_ypencil, dm%dccc)   
    call transpose_y_to_x (accc_ypencil, accc_xpencil, dm%dccc)   
    div_ccc_xpencil = div_ccc_xpencil + accc_xpencil   ! div_tmp = dux/dx + 1/r * duy/dy + 1/r^2 * duz/dz

    call transpose_x_to_y (div_ccc_xpencil, div_ccc_ypencil, dm%dccc)
    call transpose_y_to_z (div_ccc_ypencil, div_ccc_zpencil, dm%dccc)
    
    if(is_fbcx_velo_required) then
      fbcx_div_4cc(1:, :, :) = qxdx_pcc_xpencil(1:, :, :) + &
                               qydy_pcc_xpencil(1:, :, :) + &
                               qzdz_pcc_xpencil(1:, :, :)
      fbcx_div_4cc(2:, :, :) = qxdx_pcc_xpencil(dm%dpcc%xsz(1):, :, :) + &
                               qydy_pcc_xpencil(dm%dpcc%xsz(1):, :, :) + &
                               qzdz_pcc_xpencil(dm%dpcc%xsz(1):, :, :)
    end if

    if(is_fbcy_velo_required) then
      fbcy_div_c4c(:, 1, :) = qxdx_cpc_ypencil(:, 1, :) + &
                              qydy_cpc_ypencil(:, 1, :) + &
                              qzdz_cpc_ypencil(:, 1, :)
      fbcy_div_c4c(:, 2, :) = qxdx_cpc_ypencil(:, dm%dcpc%ysz(2), :) + &
                              qydy_cpc_ypencil(:, dm%dcpc%ysz(2), :) + &
                              qzdz_cpc_ypencil(:, dm%dcpc%ysz(2), :)

    end if

!==========================================================================================================
! preparation of intermediate variables to be used - thermal only
!==========================================================================================================
    !gxix_ccc_xpencil = qxix_ccc_xpencil ! default for flow only
    !gxiy_ppc_xpencil = qxiy_ppc_xpencil ! default for flow only
    !gxiz_pcp_xpencil = qxiz_pcp_xpencil ! default for flow only
    !gyix_ppc_ypencil = qyix_ppc_ypencil ! default for flow only
    !gyiy_ccc_ypencil = qyiy_ccc_ypencil ! default for flow only
    !gyiz_cpp_ypencil = qyiz_cpp_ypencil ! default for flow only
    !gzix_pcp_zpencil = qzix_pcp_zpencil ! default for flow only
    !gziy_cpp_zpencil = qziy_cpp_zpencil ! default for flow only
    !gziz_ccc_zpencil = qziz_ccc_zpencil ! default for flow only
    
    mu_ccc_xpencil    = ONE ! x-v1, default for flow only
    mu_ccc_ypencil    = ONE ! y-v2, default for flow only
    mu_ccc_zpencil    = ONE ! z-v3, default for flow only
    muixy_ppc_xpencil = ONE ! y-v1, default for flow only
    muixy_ppc_ypencil = ONE ! x-v2, default for flow only
    muixz_pcp_xpencil = ONE ! z-v1, default for flow only
    muixz_pcp_zpencil = ONE ! x-v3, default for flow only
    muiyz_cpp_ypencil = ONE ! z-v2, default for flow only
    muiyz_cpp_zpencil = ONE ! y-v3, default for flow only
    if(dm%is_thermo) then
!----------------------------------------------------------------------------------------------------------
!    gx 
!    | -[ipx]-> gxix_ccc_xpencil
!    | -[x2y]-> gx_ypencil(temp) 
!               | -[ipy]-> gxiy_ppc_ypencil(temp) -[y2x]-> gxiy_ppc_xpencil
!               | -[y2z]-> gx_zpencil(temp) 
!                          | -[ipz]-> gxiz_pcp_zpencil(temp) -[z2y]-> gxiz_pcp_ypencil(temp) -[y2x]-> gxiz_pcp_xpencil
!----------------------------------------------------------------------------------------------------------
      call Get_x_midp_P2C_3D(fl%gx, gxix_ccc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_gx)
      
      call transpose_x_to_y (fl%gx, apcc_ypencil, dm%dpcc) !gx_ypencil
      call Get_y_midp_C2P_3D(apcc_ypencil, appc_ypencil, dm, dm%iAccuracy, dm%ibcy_qx, dm%fbcy_gx)
      call transpose_y_to_x (appc_ypencil, gxiy_ppc_xpencil, dm%dppc)

      call transpose_y_to_z (apcc_ypencil, apcc_zpencil, dm%dpcc)!gx_zpencil
      call Get_z_midp_C2P_3D(apcc_zpencil, apcp_zpencil, dm, dm%iAccuracy, dm%ibcz_qx, dm%fbcz_gx)
      call transpose_z_to_y (apcp_zpencil, apcp_ypencil, dm%dpcc)!qxiz_pcp_ypencil
      call transpose_y_to_x (apcp_ypencil, gxiz_pcp_xpencil, dm%dpcp)

      if(is_fbcx_velo_required) then
        fbcx_gx_4cc(1:2, :, :) = dm%fbcx_gx(1:2, :, :)
      end if
!----------------------------------------------------------------------------------------------------------
!    gy
!    | -[ipx]-> gyix_ppc_xpencil(temp) -[x2y]-> gyix_ppc_ypencil
!    | -[x2y]-> gy_ypencil(temp) 
!               | -[ipy]-> gyiy_ccc_ypencil
!               | -[y2z]-> gy_zpencil(temp) 
!                          | -[ipz]-> gyiz_cpp_zpencil(temp) -[z2y]-> gyiz_cpp_ypencil          
!----------------------------------------------------------------------------------------------------------
      call Get_x_midp_C2P_3D(fl%gy, appc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_gy)
      call transpose_x_to_y (appc_xpencil, gyix_ppc_ypencil, dm%dppc)

      call transpose_x_to_y (fl%gy, acpc_ypencil, dm%dcpc) !acpc_ypencil = gy_ypencil
      call Get_y_midp_P2C_3D(acpc_ypencil, gyiy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_gy)

      call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc)!qy_zpencil
      call Get_z_midp_C2P_3D(acpc_zpencil, acpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_gy)
      call transpose_z_to_y(acpp_zpencil, gyiz_cpp_ypencil, dm%dcpp)

      if(is_fbcy_velo_required) then
        fbcy_gy_c4c(:, 1, :) = acpc_ypencil(:, 1,              :)
        fbcy_gy_c4c(:, 2, :) = acpc_ypencil(:, dm%dcpc%ysz(2), :)
      end if
!----------------------------------------------------------------------------------------------------------
!    gz 
!    | -[ipx]-> gzix_pcp_xpencil(temp) -[x2y]-> gzix_pcp_ypencil(temp) -[y2z]-> gzix_pcp_zpencil
!    | -[x2y]-> gz_ypencil(temp) 
!               | -[ipy]-> gziy_cpp_ypencil(temp) -[y2z]-> gziy_cpp_zpencil (no-cly)
!               | -[y2z]-> gz_zpencil(temp) 
!                          | -[ipz]-> gziz_ccc_zpencil
!----------------------------------------------------------------------------------------------------------
      call Get_x_midp_C2P_3D(fl%gz, apcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_gz)
      call transpose_x_to_y(apcp_xpencil, apcp_ypencil, dm%dpcp)
      call transpose_y_to_z(apcp_ypencil, gzix_pcp_zpencil, dm%dpcp)

      call transpose_x_to_y(fl%gz, accp_ypencil, dm%dccp) ! gz_ypencil
      
      if(dm%icoordinate == ICARTESIAN) then
        call Get_y_midp_C2P_3D(accp_ypencil, acpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_gz)
        call transpose_y_to_z(acpp_ypencil, gziy_cpp_zpencil, dm%dcpp)
      end if

      call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp) ! gz_zpencil
      call Get_z_midp_P2C_3D(accp_zpencil, gziz_ccc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_gz)
      
      if(is_fbcz_velo_required) then
        fbcz_gz_pc4(:, :, 1) = gzix_pcp_zpencil(:, :,              1)
        fbcz_gz_pc4(:, :, 2) = gzix_pcp_zpencil(:, :, dm%dpcp%zsz(3))
      end if
!----------------------------------------------------------------------------------------------------------
!    m = mu_ccc_xpencil <need BCx>
!             BCx:|-->mu_4cc_xpencil
!    | -[ipx]-> muix_pcc_xpencil(temp) -[x2y]-> muix_pcc_ypencil(temp)
!                                               | -[ipy]-> muixy_ppc_ypencil -[y2z]-> muixy_ppc_xpencil
!                                               | -[y2z]-> muix_pcc_zpencil(temp) -[ipz]-> muixz_pcp_zpencil -[z2y]-> muixz_pcp_ypencil (temp) -[y2z]-> muixz_pcp_xpencil
!    | -[x2y]-> mu_ccc_ypencil <need BCy> -[y2z]-> mu_ccc_zpencil <need BCz>
!                                             BCz:|-[ipz]->mu_ccp_zpencil -> mu_cc4_zpencil
!              | -[ipy]-> muiy_cpc_ypencil(temp) -[y2z]-> muiy_cpc_zpencil(temp) -[ipz]-> muiyz_cpp_zpencil -[z2y]-> muiyz_cpp_ypencil
!                     BCy:|->mu_c4c_ypencil                                           
!----------------------------------------------------------------------------------------------------------
      mu_ccc_xpencil = tm%mVisc
      call transpose_x_to_y(mu_ccc_xpencil, mu_ccc_ypencil, dm%dccc)
      call transpose_y_to_z(mu_ccc_ypencil, mu_ccc_zpencil, dm%dccc)

      call Get_x_midp_C2P_3D(mu_ccc_xpencil, apcc_xpencil, dm, dm%iAccuracy, dm%ibcx_Th, dm%fbcx_ftp(:, :, :)%m) ! muix_pcc_xpencil
      call transpose_x_to_y(apcc_xpencil, apcc_ypencil, dm%dpcc) ! muix_pcc_ypencil
      
      call Get_y_midp_C2P_3D(apcc_ypencil, muixy_ppc_ypencil, dm, dm%iAccuracy, dm%ibcy_Th, dm%fbcy_ftp(:, :, :)%m)
      call transpose_y_to_x(muixy_ppc_ypencil, muixy_ppc_xpencil, dm%dppc) 

      call transpose_y_to_z(apcc_ypencil, apcc_zpencil, dm%dpcc) ! muix_pcc_zpencil
      call Get_z_midp_C2P_3D(apcc_zpencil, muixz_pcp_zpencil, dm, dm%iAccuracy, dm%ibcz_Th, dm%fbcz_ftp(:, :, :)%m)
      call transpose_z_to_y(muixz_pcp_zpencil, apcp_zpencil, dm%dpcp)
      call transpose_y_to_x(apcp_zpencil, muixz_pcp_xpencil, dm%dpcp)

      call Get_y_midp_C2P_3D(mu_ccc_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcy_Th, dm%fbcy_ftp(:, :, :)%m) ! muiy_cpc_ypencil
      call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc) !muiy_cpc_zpencil
      call Get_z_midp_C2P_3D(acpc_zpencil, muiyz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_Th, dm%fbcz_ftp(:, :, :)%m)
      call transpose_z_to_y(muiyz_cpp_zpencil, muiyz_cpp_ypencil, dm%dcpp)

      if(is_fbcx_velo_required) then
        fbcx_mu_4cc(1, :, :) = apcc_xpencil(1, :, :)
        fbcx_mu_4cc(2, :, :) = apcc_xpencil(dm%dpcc%xsz(1), :, :)
      end if

    end if
    mu_ccc_xpencil    =  fl%rre * mu_ccc_xpencil    
    mu_ccc_ypencil    =  fl%rre * mu_ccc_ypencil    
    mu_ccc_zpencil    =  fl%rre * mu_ccc_zpencil    
    muixy_ppc_xpencil =  fl%rre * muixy_ppc_xpencil 
    muixy_ppc_ypencil =  fl%rre * muixy_ppc_ypencil 
    muixz_pcp_xpencil =  fl%rre * muixz_pcp_xpencil 
    muixz_pcp_zpencil =  fl%rre * muixz_pcp_zpencil 
    muiyz_cpp_ypencil =  fl%rre * muiyz_cpp_ypencil 
    muiyz_cpp_zpencil =  fl%rre * muiyz_cpp_zpencil 

!==========================================================================================================
! preparation of intermediate variables to be used - cylindrical only
!==========================================================================================================
    if(dm%icoordinate == ICYLINDRICAL) then
!----------------------------------------------------------------------------------------------------------
!    qy/r=qyr
!    | -[x2y]-> qyr_ypencil(temp) 
!               | -[ipy]-> qyriy_ccc_ypencil
!               | -[1dy]-> qyrdy_ccc_ypencil <need BCy>
!           BCy:| -[1dy] ->qyrdy_cpc_ypencil -> qyrdy_c4c_ypencil
!               | -[y2z]-> qyr_zpencil(temp) 
!                          | -[ipz]-> qyriz_cpp_zpencil -[z2y]-> qyriz_cpp_ypencil (no-thermal)        
!                          | -[1dz]-> qyrdz_cpp_zpencil -[z2y]-> qyrdz_cpp_ypencil
!----------------------------------------------------------------------------------------------------------
      acpc_xpencil = fl%qy
      call multiple_cylindrical_rn(acpc_xpencil, dm%dcpc, dm%rpi, 1, IPENCIL(1)) ! qr/r

      call transpose_x_to_y(acpc_xpencil, acpc_ypencil, dm%dcpc) ! acpc_ypencil = qyr_ypencil
      call Get_y_midp_P2C_3D(acpc_ypencil, qyriy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr)
      call Get_y_1st_derivative_P2C_3D(acpc_ypencil, qyrdy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr)
      if( any(dm%ibcy_qy == IBC_DIRICHLET) )then
        call Get_y_1st_derivative_P2P_3D(acpc_ypencil, acpc_ypencil1, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr
        call get_dirichlet_geo_bcy(dm%ibcy_qy, acpc_ypencil1, qyrdy_c4c_ypencil)
      end if

      call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc)
      call Get_z_midp_C2P_3D(acpc_zpencil, qyriz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qyr)
      call Get_z_1st_derivative_C2P_3D(acpc_zpencil, qyrdz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qyr)

      if(.not. dm%is_thermo) call transpose_z_to_y(qyriz_cpp_zpencil, qyriz_cpp_ypencil, dm%dcpp)
      call transpose_z_to_y(qyrdz_cpp_zpencil, qyrdz_cpp_ypencil, dm%dcpp)

      if(is_fbcy_velo_required) then
        fbcy_qyr_c4c(:, 1, :) = acpc_ypencil(:, 1,              :)
        fbcy_qyr_c4c(:, 2, :) = acpc_ypencil(:, dm%dcpc%ysz(2), :)

        call Get_y_1st_derivative_P2P_3D(acpc_ypencil, acpc_ypencil1, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr)
        fbcy_qyrdy_c4c(:, 1, :) = acpc_ypencil1(:, 1,              :)
        fbcy_qyrdy_c4c(:, 2, :) = acpc_ypencil1(:, dm%dcpc%ysz(2), :)
      end if
!----------------------------------------------------------------------------------------------------------
!    qz/r=qzr
!    | -[x2y]-> qzr_ypencil(temp) 
!               | -[ipy]-> qzriy_cpp_ypencil -[y2z]-> qzriy_cpp_zpencil
!               | -[1dy]-> qzrdy_cpp_ypencil -[y2z]-> qzrdy_cpp_zpencil
!               | -[y2z]-> qzr_zpencil(temp) -[ipz]-> qzriz_ccc_zpencil -[z2y]- qzriz_ccc_ypencil
!----------------------------------------------------------------------------------------------------------
      accp_xpencil = fl%qz
      call multiple_cylindrical_rn(accp_xpencil, dm%dccp, dm%rpi, 1, IPENCIL(1)) ! qz/r

      call transpose_x_to_y(accp_xpencil, accp_ypencil, dm%dccp)

      call Get_y_midp_C2P_3D(accp_ypencil, qzriy_cpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qzr(:, :, :))
      call transpose_y_to_z(qzriy_cpp_ypencil, qzriy_cpp_zpencil, dm%dcpp)

      call Get_y_1st_derivative_C2P_3D(accp_ypencil, qzrdy_cpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qzr(:, :, :))
      call transpose_y_to_z(qzrdy_cpp_ypencil, qzrdy_cpp_zpencil, dm%dcpp)

      call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
      call Get_z_midp_P2C_3D(accp_zpencil, qzriz_ccc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qzr(:, :, :))
      call transpose_z_to_y(qzriz_ccc_zpencil, qzriz_ccc_ypencil, dm%dccc)

      if(is_fbcy_velo_required) then
        call Get_y_midp_C2P_3D(qzriz_ccc_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qzr(:, :, :))
        fbcy_qzr_c4c(:, 1, :) = acpc_ypencil(:, 1,              :) 
        fbcy_qzr_c4c(:, 2, :) = acpc_ypencil(:, dm%dcpc%ysz(2), :) 
      end if

      if(dm%is_thermo) then
!----------------------------------------------------------------------------------------------------------
!    gy/r=gyr
!    | -[x2y]-> gyr_ypencil(temp) -[y2z]-> gyr_zpencil(temp) -[ipz]-> gyriz_cpp_zpencil(temp) -[z2y]-> gyriz_cpp_ypencil       
!----------------------------------------------------------------------------------------------------------
        acpc_xpencil = fl%gy
        call multiple_cylindrical_rn(acpc_xpencil, dm%dcpc, dm%rpi, 1, IPENCIL(1)) ! gr/r
        
        call transpose_x_to_y(acpc_xpencil, acpc_ypencil, dm%dcpc)
        call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc)
        call Get_z_midp_C2P_3D(acpc_zpencil, acpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_gyr(:, :, :))
        call transpose_z_to_y(acpp_zpencil, gyriz_cpp_ypencil, dm%dcpp)
!----------------------------------------------------------------------------------------------------------
!    gz/r=gzr
!    | -[x2y]-> gzr_ypencil(temp) 
!               | -[ipy]-> gzriy_cpp_ypencil (temp) -[y2z]-> gzriy_cpp_zpencil
!               | -[y2z]-> qzr_zpencil(temp) -[ipz]-> qzriz_ccc_zpencil (temp) -[y2z]-> gzriz_ccc_ypencil     
!----------------------------------------------------------------------------------------------------------
        accp_xpencil = fl%gz
        call multiple_cylindrical_rn(accp_xpencil, dm%dccp, dm%rpi, 1, IPENCIL(1)) ! gz/r

        call transpose_x_to_y(accp_xpencil, accp_ypencil, dm%dccp)

        call Get_y_midp_C2P_3D(accp_ypencil, acpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_gzr(:, :, :))
        call transpose_y_to_z(acpp_ypencil, gzriy_cpp_zpencil, dm%dccp)

        call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
        call Get_z_midp_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_gzr(:, :, :))
        call transpose_z_to_y(accc_zpencil, gzriz_ccc_ypencil, dm%dccc)

        if(is_fbcy_velo_required) then
          call Get_y_midp_C2P_3D(gzriz_ccc_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_gzr(:, :, :))
          fbcy_gzr_c4c(:, 1, :) = acpc_ypencil(:, 1,              :) 
          fbcy_gzr_c4c(:, 2, :) = acpc_ypencil(:, dm%dcpc%ysz(2), :) 
        end if
      end if

    end if

!==========================================================================================================
! the RHS of x-momentum equation
!==========================================================================================================
    i = 1
    fl%mx_rhs      = ZERO
    mx_rhs_ypencil = ZERO
    mx_rhs_zpencil = ZERO
    mx_rhs_pfc_xpencil = ZERO
!----------------------------------------------------------------------------------------------------------
! X-mom convection term (x-c1/3), X-pencil: -d(gx^x * qx^x)/dx at (i', j, k)
!----------------------------------------------------------------------------------------------------------  
    if ( .not. dm%is_thermo) then
      accc_xpencil = qxix_ccc_xpencil
      if(is_fbcx_velo_required) fbcx_4cc = -fbcx_qx_4cc * fbcx_qx_4cc
    else
      accc_xpencil = gxix_ccc_xpencil
      if(is_fbcx_velo_required) fbcx_4cc = -fbcx_gx_4cc * fbcx_qx_4cc
    end if
    accc_xpencil = -accc_xpencil * qxix_ccc_xpencil
    call Get_x_1st_derivative_C2P_3D(accc_xpencil, apcc_xpencil, dm, dm%iAccuracy, mbcx_cov1, fbcx_4cc)
    fl%mx_rhs = fl%mx_rhs + apcc_xpencil
!----------------------------------------------------------------------------------------------------------
! X-mom convection term (x-c2/3), Y-pencil: -d(<gy>^x * <qx>^y)/dy * (1/r) at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      appc_ypencil = qyix_ppc_ypencil
    else
      appc_ypencil = gyix_ppc_ypencil
    end if
    appc_ypencil = -appc_ypencil * qxiy_ppc_ypencil
    call Get_y_1st_derivative_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%iAccuracy, mbcy_cov1)
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_ypencil, dm%dpcc, dm%rci, 1, IPENCIL(2))
    mx_rhs_ypencil = mx_rhs_ypencil + apcc_ypencil
!----------------------------------------------------------------------------------------------------------
! X-mom convection term (x-c3/3), Z-pencil: -d(<gz>^x * <qx>^z)/dz * (1/r^2) at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      apcp_zpencil = qzix_pcp_zpencil
    else
      apcp_zpencil = gzix_pcp_zpencil
    end if
    apcp_zpencil = -apcp_zpencil * qxiz_pcp_zpencil
    call Get_z_1st_derivative_P2C_3D(apcp_zpencil, apcc_zpencil, dm, dm%iAccuracy, mbcz_cov1)
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_zpencil, dm%dpcc, dm%rci, 2, IPENCIL(3))
    mx_rhs_zpencil = mx_rhs_zpencil + apcc_zpencil

#ifdef DEBUG_STEPS
    apcc_test = fl%mx_rhs
    call transpose_y_to_x (mx_rhs_ypencil, apcc_xpencil, dm%dpcc)
    apcc_test =  apcc_test + apcc_xpencil
    call transpose_z_to_y (mx_rhs_zpencil, apcc_ypencil, dm%dpcc)
    call transpose_y_to_x (apcc_ypencil, apcc_xpencil, dm%dpcc)
    apcc_test =  apcc_test + apcc_xpencil 
    call wrt_3d_pt_debug(apcc_test,  dm%dpcc, fl%iteration, isub, '', 'ConX@bf st') ! debug_ww
#endif
!----------------------------------------------------------------------------------------------------------
! X-mom pressure gradient in x direction, X-pencil, d(sigma_1 p)
!----------------------------------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3D( -fl%pres, apcc_xpencil, dm, dm%iAccuracy, dm%ibcx_pr, -dm%fbcx_pr)
    mx_rhs_pfc_xpencil=  mx_rhs_pfc_xpencil + apcc_xpencil
!----------------------------------------------------------------------------------------------------------
! X-mom gravity in x direction, X-pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%is_thermo .and. (fl%igravity == i .or. fl%igravity == -i) )  then
      call Get_x_midp_C2P_3D(fl%dDens, apcc_xpencil, dm, dm%iAccuracy, dm%ibcx_Th, dm%ftpbcx_var(:, :, :)%d )
      mx_rhs_pfc_xpencil =  mx_rhs_pfc_xpencil + fl%fgravity(i) * apcc_xpencil
    end if
!----------------------------------------------------------------------------------------------------------
! X-mom diffusion term (x-v-1/3), X-pencil, d_x [2 * mu * (d_x qx - 1/3*div) ] at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    accc_xpencil = TWO * mu_ccc_xpencil * ( qxdx_ccc_xpencil - ONE_THIRD * div_ccc_xpencil)
    if(is_fbcx_velo_required) then
      fbcx_4cc = TWO * fbcx_mu_4cc * ( fbcx_qxdx_4cc - ONE_THIRD * fbcx_div_4cc)
    end if
    call Get_x_1st_derivative_C2P_3D(accc_xpencil,  apcc_xpencil, dm, dm%iAccuracy, mbcx_tau1, fbcx_4cc) 
    fl%mx_rhs = fl%mx_rhs + apcc_xpencil
!----------------------------------------------------------------------------------------------------------
! X-mom diffusion term (x-v-2/3), Y-pencil, (1/r) * d_y [mu^{x,y} * (d_x q_y + r * d_y q_x) ] at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    appc_ypencil = qxdy_ppc_ypencil
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(appc_ypencil, dm%dppc, ONE/dm%rpi, 1, IPENCIL(2))
    appc_ypencil = ( appc_ypencil + qydx_ppc_ypencil) * muixy_ppc_ypencil
    call Get_y_1st_derivative_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%iAccuracy, mbcy_tau1)
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_ypencil, dm%dpcc, dm%rci, 1, IPENCIL(2))
    mx_rhs_ypencil =  mx_rhs_ypencil + apcc_ypencil
!----------------------------------------------------------------------------------------------------------
! X-mom diffusion term (x-v-3/3), Z-pencil, (1/r^2) * d_z [mu^{x,z} * (d_x q_z + d_z q_x) ] at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    apcp_zpencil = ( qxdz_pcp_zpencil + qzdx_pcp_zpencil) * muixz_pcp_zpencil
    call Get_z_1st_derivative_P2C_3D(apcp_zpencil, apcc_zpencil, dm, dm%iAccuracy, mbcz_tau1)
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_zpencil, dm%dpcc, dm%rci, 2, IPENCIL(3))
    mx_rhs_zpencil =  mx_rhs_zpencil + apcc_zpencil
!----------------------------------------------------------------------------------------------------------
! x-mom: convert all terms to rhs
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_x (mx_rhs_ypencil, apcc_xpencil, dm%dpcc)
    fl%mx_rhs =  fl%mx_rhs + apcc_xpencil
    call transpose_z_to_y (mx_rhs_zpencil, apcc_ypencil, dm%dpcc)
    call transpose_y_to_x (apcc_ypencil, apcc_xpencil, dm%dpcc)
    fl%mx_rhs =  fl%mx_rhs + apcc_xpencil
    
!==========================================================================================================
! the RHS of y-momentum equation
!==========================================================================================================
    i = 2
    fl%my_rhs          = ZERO
    my_rhs_ypencil     = ZERO
    my_rhs_zpencil     = ZERO
    my_rhs_pfc_xpencil = ZERO
    my_rhs_pfc_ypencil = ZERO
!----------------------------------------------------------------------------------------------------------
! Y-mom convection term (y-c1/4), X-pencil: -d(gx^y * qy^x)/dx at (i, j', k)
!----------------------------------------------------------------------------------------------------------  
    if ( .not. dm%is_thermo) then
      appc_xpencil = qxiy_ppc_xpencil
    else
      appc_xpencil = gxiy_ppc_xpencil
    end if
    appc_xpencil = - appc_xpencil * qyix_ppc_xpencil
    call Get_x_1st_derivative_P2C_3D( appc_xpencil, acpc_xpencil, dm, dm%iAccuracy, mbcx_cov2)
    fl%my_rhs = fl%my_rhs + acpc_xpencil
!----------------------------------------------------------------------------------------------------------
! Y-mom convection term (y-c2/4), Y-pencil: -d(<gy>^y * <qy/r>^r)/dy at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      accc_ypencil1 = qyriy_ccc_ypencil
      if(is_fbcy_velo_required) fbcy_c4c = fbcy_qyr_c4c
    else
      accc_ypencil1 = qyiy_ccc_ypencil
      if(is_fbcy_velo_required) fbcy_c4c = fbcy_qy_c4c
    end if

    if ( .not. dm%is_thermo) then
      accc_ypencil = qyiy_ccc_ypencil
      if(is_fbcy_velo_required) fbcy_c4c = - fbcy_qy_c4c * fbcy_c4c
    else
      accc_ypencil = gyiy_ccc_ypencil
      if(is_fbcy_velo_required) fbcy_c4c = - fbcy_gy_c4c * fbcy_c4c
    end if

    accc_ypencil = - accc_ypencil1 * accc_ypencil
    call Get_y_1st_derivative_C2P_3D( accc_ypencil, acpc_ypencil, dm, dm%iAccuracy, mbcy_cov2, fbcy_c4c)
    my_rhs_ypencil = my_rhs_ypencil + acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Y-mom convection term (y-c3/4), Z-pencil: -d(<gz/r>^r * <qr/r>^z)/dz at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      acpp_zpencil1 = qyriz_cpp_zpencil
      if ( .not. dm%is_thermo) then
        acpp_zpencil = qzriy_cpp_zpencil
      else
        acpp_zpencil = gzriy_cpp_zpencil
      end if
    else
      acpp_zpencil1 = qyiz_cpp_zpencil
      if ( .not. dm%is_thermo) then
        acpp_zpencil = qziy_cpp_zpencil
      else
        acpp_zpencil = gziy_cpp_zpencil
      end if
    end if
    acpp_zpencil = - acpp_zpencil1 * acpp_zpencil
    call Get_z_1st_derivative_P2C_3D( acpp_zpencil, acpc_zpencil, dm, dm%iAccuracy, mbcz_cov2)
    my_rhs_zpencil = my_rhs_zpencil + acpc_zpencil
!----------------------------------------------------------------------------------------------------------
! Y-mom convection term (y-c4/4), Y-pencil: < <gz/r>^z * <qz/r>^z >^y
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      if ( .not. dm%is_thermo) then
        accc_ypencil = qzriz_ccc_ypencil
        if(is_fbcy_velo_required) fbcy_c4c = fbcy_qzr_c4c * fbcy_qzr_c4c
      else
        accc_ypencil = gzriz_ccc_ypencil
        if(is_fbcy_velo_required) fbcy_c4c = fbcy_gzr_c4c * fbcy_qzr_c4c
      end if
      accc_ypencil = accc_ypencil * qzriz_ccc_ypencil 
      
      call Get_y_midp_C2P_3D( accc_ypencil, acpc_ypencil, dm, dm%iAccuracy, mbcr_cov2, fbcy_c4c)
      my_rhs_ypencil = my_rhs_ypencil + acpc_ypencil
    end if
!----------------------------------------------------------------------------------------------------------
! Y-mom pressure gradient in y direction, Y-pencil, d(sigma_1 p)
!----------------------------------------------------------------------------------------------------------
    call Get_y_1st_derivative_C2P_3D( -pres_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcy_pr, -dm%fbcy_pr)
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(acpc_ypencil, dm%dcpc, ONE/dm%rpi, 1, IPENCIL(2))
    my_rhs_pfc_ypencil =  my_rhs_pfc_ypencil + acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Y-mom gravity in y direction, Y-pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%is_thermo .and. (fl%igravity == i .or. fl%igravity == -i) )  then
      call Get_y_midp_C2P_3D(dDens_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcy_Th, dm%fbcy_ftp(:, :, :)%d )
      my_rhs_pfc_ypencil =  my_rhs_pfc_ypencil + fl%fgravity(i) * acpc_ypencil
    end if
!----------------------------------------------------------------------------------------------------------
! Y-mom diffusion term (y-v-1/4), X-pencil, d_x [mu^{x,y} * (d_x q_y + r * d_y q_x) ] at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    appc_xpencil1 = qxdy_ppc_xpencil
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(appc_xpencil1, dm%dppc, ONE/dm%rpi, 1, IPENCIL(1))
    appc_xpencil = ( appc_xpencil1 + qydx_ppc_xpencil) * muixy_ppc_xpencil
    call Get_x_1st_derivative_P2C_3D(appc_xpencil, acpc_xpencil, dm, dm%iAccuracy, mbcx_tau2)
    fl%my_rhs = fl%my_rhs + acpc_xpencil
!----------------------------------------------------------------------------------------------------------
! Y-mom diffusion term (y-v-2/4), Y-pencil, d_y [2 * r * mu * (d_y (qy/r) - 1/3*div) ] at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      accc_ypencil1 = qyrdy_ccc_ypencil
      if(is_fbcy_velo_required) fbcy_c4c = fbcy_qyrdy_c4c
    else
      accc_ypencil1 = qydy_ccc_ypencil
      if(is_fbcy_velo_required) fbcy_c4c = fbcy_qydy_c4c
    end if
    
    accc_ypencil = ( accc_ypencil1 - ONE_THIRD * div_ccc_ypencil) * TWO * mu_ccc_ypencil
    if(is_fbcy_velo_required) then 
      fbcy_c4c = ( fbcy_c4c - ONE_THIRD * fbcy_div_c4c) * TWO * fbcy_mu_c4c
    end if
    
    if(dm%icoordinate == ICYLINDRICAL) then
      call multiple_cylindrical_rn(accc_ypencil, dm%dccc, ONE/dm%rci, 1, IPENCIL(2))
      if(is_fbcy_velo_required) &
      call multiple_cylindrical_rn_x4x(fbcy_c4c, dm%dcpc, ONE/dm%rci, 1, IPENCIL(2))
    end if
    call Get_y_1st_derivative_C2P_3D(accc_ypencil,  acpc_ypencil, dm, dm%iAccuracy, mbcy_tau2, fbcy_c4c) 
    my_rhs_ypencil = my_rhs_ypencil + acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Y-mom diffusion term (y-v-3/4), Z-pencil, d_z [mu^{y,z} * (dy(qz/r) + 1/r dz(qr/r) - 1/r * (qz/r)^r)] at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      acpp_zpencil  = qzrdy_cpp_zpencil ! acpp_zpencil  = d(qz/r)/dy
      acpp_zpencil1 = qyrdz_cpp_zpencil ! acpp_zpencil1 = d(qr/r)/dz
      call multiple_cylindrical_rn(acpp_zpencil1, dm%dcpp, dm%rpi, 1, IPENCIL(3)) ! acpp_zpencil1 = 1/r * d(qr/r)/dz
      acpp_zpencil = acpp_zpencil + acpp_zpencil1
      acpp_zpencil1 = qzriy_cpp_zpencil ! acpp_zpencil1 = (qz/r)^r
      call multiple_cylindrical_rn(acpp_zpencil1, dm%dcpp, dm%rpi, 1, IPENCIL(3)) ! acpp_zpencil1 = 1/r * (dz/r)^r
      acpp_zpencil = acpp_zpencil - acpp_zpencil1
    else
      acpp_zpencil  = qzdy_cpp_zpencil + qydz_cpp_zpencil - ZERO
    end if
    
    acpp_zpencil = acpp_zpencil * muiyz_cpp_zpencil
    call Get_z_1st_derivative_P2C_3D(acpp_zpencil, acpc_zpencil, dm, dm%iAccuracy, mbcz_tau2)
    my_rhs_zpencil =  my_rhs_zpencil + acpc_zpencil
!----------------------------------------------------------------------------------------------------------
! Y-mom diffusion term (y-v-4/4), Y-pencil, ( -2 * mu * [ 1/r^2 * d(qz)/dz - 1/3 div + 1/r * (qr/r)^r ] )^r at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      accc_ypencil = qzdz_ccc_ypencil
      call multiple_cylindrical_rn(accc_ypencil, dm%dccc, dm%rci, 2, IPENCIL(2)) ! accc_ypencil = 1/r^2 * d(qz)/dz
      if(is_fbcy_velo_required) then
        fbcy_c4c = fbcy_qzdz_c4c
        call multiple_cylindrical_rn_x4x(fbcy_c4c, dm%dcpc, dm%rci, 2, IPENCIL(2))
      end if
      
      accc_ypencil1 = qyriy_ccc_ypencil
      call multiple_cylindrical_rn(accc_ypencil1, dm%dccc, dm%rci, 1, IPENCIL(2)) !accc_ypencil1 = 1/r * (qr/r)^r
      if(is_fbcy_velo_required) then
        fbcy_c4c1 = fbcy_qyr_c4c
        call multiple_cylindrical_rn_x4x(fbcy_c4c1, dm%dcpc, dm%rci, 1, IPENCIL(2))
        fbcy_c4c = (fbcy_c4c - ONE_THIRD * fbcy_div_c4c + fbcy_c4c1) * TWO * fbcy_mu_c4c
      end if
      
      accc_ypencil = (accc_ypencil - ONE_THIRD * div_ccc_ypencil + accc_ypencil1) * TWO * mu_ccc_ypencil

      call Get_y_midp_C2P_3D(accc_ypencil, acpc_ypencil, dm, dm%iAccuracy, mbcr_tau2, fbcy_c4c)
      my_rhs_ypencil =  my_rhs_ypencil - acpc_ypencil
    end if
!----------------------------------------------------------------------------------------------------------
! Y-mom: convert all terms to rhs
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_x (my_rhs_ypencil, acpc_xpencil, dm%dcpc)
    fl%my_rhs =  fl%my_rhs + acpc_xpencil

    call transpose_z_to_y (my_rhs_zpencil, acpc_ypencil, dm%dcpc)
    call transpose_y_to_x (acpc_ypencil, acpc_xpencil, dm%dcpc)
    fl%my_rhs =  fl%my_rhs + acpc_xpencil

    call transpose_y_to_x (my_rhs_pfc_ypencil,  my_rhs_pfc_xpencil,  dm%dcpc)
!==========================================================================================================
! the RHS of z-momentum equation
!==========================================================================================================
    i = 3
    fl%mz_rhs          = ZERO
    mz_rhs_ypencil     = ZERO
    mz_rhs_zpencil     = ZERO
    mz_rhs_pfc_xpencil = ZERO
    mz_rhs_pfc_ypencil = ZERO
    mz_rhs_pfc_zpencil = ZERO
!----------------------------------------------------------------------------------------------------------
! Z-mom convection term (Z-c1/4), X-pencil: -d(gx^Z * qz^x)/dx at (i, j, k')
!----------------------------------------------------------------------------------------------------------  
    if ( .not. dm%is_thermo) then
      apcp_xpencil = qxiz_pcp_xpencil
    else
      apcp_xpencil = gxiz_pcp_xpencil
      if(is_fbcy_velo_required) 
    end if
    apcp_xpencil = - apcp_xpencil * qzix_pcp_xpencil
    call Get_x_1st_derivative_P2C_3D( appc_xpencil, accp_xpencil, dm, dm%iAccuracy, mbcx_cov3)
    fl%mz_rhs = fl%mz_rhs + accp_xpencil
!----------------------------------------------------------------------------------------------------------
! Z-mom convection term (y-c2/4), Y-pencil: -d(<gy>^z * <qz/r>^r)/dy at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      acpp_ypencil1 = qzriy_cpp_ypencil
    else
      acpp_ypencil1 = qziy_cpp_ypencil
    end if
    if ( .not. dm%is_thermo) then
      acpp_ypencil = qyiz_cpp_ypencil
    else
      acpp_ypencil = gyiz_cpp_ypencil
    end if
    acpp_ypencil = - acpp_ypencil * acpp_ypencil1
    call Get_y_1st_derivative_P2C_3D( acpp_ypencil, accp_ypencil, dm, dm%iAccuracy, mbcy_cov3)
    mz_rhs_ypencil = mz_rhs_ypencil + accp_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-mom convection term (z-c3/4), Z-pencil: -d(<gz>^z * <qz>^z)/dz * (1/r^2) at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      accc_zpencil = qziz_ccc_zpencil
      if(is_fbcz_velo_required) fbcz_cc4 = -fbcz_qz_cc4 * fbcz_qz_cc4
    else
      accc_zpencil = gziz_ccc_zpencil
      if(is_fbcz_velo_required) fbcz_cc4 = -fbcz_gz_cc4 * fbcz_qz_cc4
    end if
    accc_zpencil = -accc_zpencil * qziz_ccc_zpencil
    call Get_z_1st_derivative_C2P_3D( accc_zpencil, accp_zpencil, dm, dm%iAccuracy, mbcz_cov3, fbcz_cc4 )
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(accp_zpencil, dm%dccp, dm%rci, 2, IPENCIL(3))
    mz_rhs_zpencil = mz_rhs_zpencil + accp_zpencil
!----------------------------------------------------------------------------------------------------------
! Z-mom convection term (z-c4/4), Y-pencil: < <gy/r>^z * <qz/r>^r >^y
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      if ( .not. dm%is_thermo) then
        acpp_ypencil = qyriz_cpp_ypencil
      else
        acpp_ypencil = gyriz_cpp_ypencil
      end if
      acpp_ypencil = - acpp_ypencil * qzriy_cpp_ypencil
      call Get_y_midp_P2C_3D( acpp_ypencil, accp_ypencil, dm, dm%iAccuracy, mbcr_cov3)
      mz_rhs_ypencil = mz_rhs_ypencil + accp_ypencil
    end if
!----------------------------------------------------------------------------------------------------------
! Z-mom pressure gradient in z direction, Z-pencil, d(sigma_1 p)
!----------------------------------------------------------------------------------------------------------
    call Get_z_1st_derivative_C2P_3D( -pres_zpencil, accp_zpencil, dm, dm%iAccuracy, dm%ibcz_pr), -dm%fbcz_pr)
    mz_rhs_pfc_zpencil =  mz_rhs_pfc_zpencil + accp_zpencil
!----------------------------------------------------------------------------------------------------------
! Z-mom gravity in z direction, Z-pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%is_thermo .and. (fl%igravity == i .or. fl%igravity == -i) )  then
      call Get_z_midp_C2P_3D(dDens_zpencil, accp_zpencil, dm, dm%iAccuracy, dm%ibcz(:, 5), dm%fbcz_ftp(:, :, :)%d )
      mz_rhs_pfc_zpencil =  mz_rhs_pfc_zpencil + fl%fgravity(i) * accp_zpencil
    end if
!----------------------------------------------------------------------------------------------------------
! Z-mom diffusion term (z-v-1/4), X-pencil, d_x [mu^{x,z} * (d_x q_z + d_z q_x) ] at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    apcp_xpencil = ( qxdz_pcp_xpencil + qzdx_pcp_xpencil) * muixz_pcp_xpencil
    call Get_x_1st_derivative_P2C_3D(apcp_xpencil, accp_xpencil, dm, dm%iAccuracy, mbcx_tau3)
    fl%mz_rhs = fl%mz_rhs + accp_xpencil
!----------------------------------------------------------------------------------------------------------
! Z-mom diffusion term (z-v-2/4), Y-pencil, d_y [mu^{y,z} * (r * dy(qz/r) + dz(qr/r) - (qz/r)^r)] at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      acpp_ypencil = qzrdy_cpp_ypencil
      call multiple_cylindrical_rn(acpp_ypencil, dm%dcpp, ONE/dm%rpi, 1, IPENCIL(2))

      acpp_ypencil = acpp_ypencil + qyrdz_cpp_ypencil
      acpp_ypencil = acpp_ypencil - qzriy_cpp_ypencil
    else
      acpp_ypencil = qzdy_cpp_ypencil
      acpp_ypencil = acpp_ypencil + qydz_cpp_ypencil
      acpp_ypencil = acpp_ypencil - ZERO
    end if
    acpp_ypencil = acpp_ypencil * muiyz_cpp_ypencil
    call Get_y_1st_derivative_P2C_3D(acpp_ypencil, accp_ypencil, dm, dm%iAccuracy, mbcy_tau3)
    mz_rhs_ypencil =  mz_rhs_ypencil + accp_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-mom diffusion term (z-v-3/4), Z-pencil, d_z [2 * mu * (1/r^2 * d_z (qz) - 1/3*div + 1/r * (qy/r)^y) ] at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      accc_zpencil = qzdz_ccc_zpencil
      call multiple_cylindrical_rn(accc_zpencil, dm%dccc, dm%rci, 2, IPENCIL(3))
      accc_zpencil1 = qyriy_ccc_zpencil
      call multiple_cylindrical_rn(accc_zpencil1, dm%dccc, dm%rci, 1, IPENCIL(3))
      if(is_fbcz_velo_required) then
        fbcz_cc4 = fbcz_qzdz_cc4
        call multiple_cylindrical_rn_xx4(fbcz_cc4, dm%dccp, dm%rci, 2, IPENCIL(3))
        fbcz_cc41 = fbcz_qyr_cc4
        call multiple_cylindrical_rn_xx4(fbcz_cc41, dm%dccp, dm%rci, 1, IPENCIL(3))
      end if
    else
      accc_zpencil  = qzdz_ccc_zpencil
      accc_zpencil1 = ZERO
      if(is_fbcz_velo_required) then
        fbcz_cc4 = fbcz_qzdz_cc4
        fbcz_cc41 = ZERO
      end if
    end if

    accc_zpencil = ( accc_zpencil - ONE_THIRD * div_ccc_zpencil + accc_zpencil1) * TWO * mu_ccc_zpencil
    if(is_fbcz_velo_required) then
        fbcz_cc4 = ( fbcz_cc4 - ONE_THIRD * fbcz_div_cc4 + fbcz_cc41) * TWO * fbcz_mu_cc4
    end if
    call Get_z_1st_derivative_C2P_3D(accc_zpencil, accp_zpencil, dm, dm%iAccuracy, mbcz_tau3, fbcz_cc4) 
    mz_rhs_zpencil = mz_rhs_zpencil + accp_zpencil

!----------------------------------------------------------------------------------------------------------
! Z-mom diffusion term (z-v-4/4), Y-pencil, (mu^{y,z} * [ d(qz/r)/dy + 1/r^2 * d(qy)/dz - 1/r * (qz/r)^y ] )^r at (i, j, k')
!---------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      acpp_ypencil = qzrdy_cpp_ypencil
      acpp_ypencil1 = qydz_cpp_ypencil
      call multiple_cylindrical_rn(acpp_ypencil1, dm%dcpp, dm%rpi, 2, IPENCIL(2))
      acpp_ypencil = acpp_ypencil + acpp_ypencil1
      acpp_ypencil1 = qzriy_cpp_ypencil
      call multiple_cylindrical_rn(acpp_ypencil1, dm%dcpp, dm%rpi, 1, IPENCIL(2))
      acpp_ypencil = (acpp_ypencil - acpp_ypencil1) * muiyz_cpp_ypencil

      call Get_y_midp_P2C_3D(acpp_ypencil, accp_ypencil, dm, dm%iAccuracy, mbcr_tau3)
      mz_rhs_ypencil =  mz_rhs_ypencil + accp_ypencil
    end if
!----------------------------------------------------------------------------------------------------------
! Z-mom: convert all terms to rhs
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_x (mz_rhs_ypencil, accp_xpencil, dm%dccp)
    fl%mz_rhs =  fl%mz_rhs + accp_xpencil

    call transpose_z_to_y (mz_rhs_zpencil, accp_ypencil, dm%dccp)
    call transpose_y_to_x (accp_ypencil, accp_xpencil, dm%dccp)
    fl%mz_rhs =  fl%mz_rhs + accp_xpencil

    call transpose_z_to_y (mz_rhs_pfc_zpencil, accp_ypencil, dm%dccp)
    call transpose_y_to_x (accp_ypencil, mz_rhs_pfc_xpencil, dm%dccp)

!==========================================================================================================
! x-pencil : to build up rhs in total, in all directions
!==========================================================================================================
    !write(*,*) nrank,  'test-7'
!----------------------------------------------------------------------------------------------------------
! x-pencil : x-momentum
!----------------------------------------------------------------------------------------------------------
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%mx_rhs, dm%dpcc, fl%iteration, isub, 'ConVisX', '@bf stepping') ! debug_ww
    write(*,*) 'fl%mx_rhs', fl%mx_rhs(:, 1, 1), fl%mx_rhs(:, 8, 8)
#endif
    call Calculate_momentum_fractional_step(fl%mx_rhs0, fl%mx_rhs, mx_rhs_pfc_xpencil, dm%dpcc, dm, isub)  
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%mx_rhs, dm%dpcc, fl%iteration, isub, 'ConVisX', '@af stepping') ! debug_ww
#endif
!----------------------------------------------------------------------------------------------------------
! x-pencil : y-momentum
!----------------------------------------------------------------------------------------------------------
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%my_rhs, dm%dcpc, fl%iteration, isub, 'ConVisY', '@bf stepping') ! debug_ww
    write(*,*) 'fl%my_rhs', fl%my_rhs(:, 1, 1), fl%my_rhs(:, 8, 8)
#endif
    call Calculate_momentum_fractional_step(fl%my_rhs0, fl%my_rhs, my_rhs_pfc_xpencil, dm%dcpc, dm, isub)
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%my_rhs, dm%dcpc, fl%iteration, isub, 'ConVisY', '@af stepping') ! debug_ww
#endif
!----------------------------------------------------------------------------------------------------------
! x-pencil : z-momentum
!----------------------------------------------------------------------------------------------------------
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%mz_rhs, dm%dccp, fl%iteration, isub, 'ConVisZ', '@bf stepping') ! debug_ww
    write(*,*) 'fl%mz_rhs', fl%mz_rhs(:, 1, 1), fl%mz_rhs(:, 8, 8)
#endif
    call Calculate_momentum_fractional_step(fl%mz_rhs0, fl%mz_rhs, mz_rhs_pfc_xpencil, dm%dccp, dm, isub)
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%mz_rhs, dm%dccp, fl%iteration, isub, 'ConVisZ', '@af stepping') ! debug_ww
#endif
!==========================================================================================================
! x-pencil : flow drive terms (source terms) in periodic Streamwise flow
!==========================================================================================================
    if (fl%idriven == IDRVF_X_MASSFLUX) then
      call Get_volumetric_average_3d_for_var_xcx(dm, dm%dpcc, fl%mx_rhs, rhsx_bulk, "mx_rhs")
      !call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), dm%fbcy_var(:, :, :, 1), dm, dm%dpcc, fl%mx_rhs, rhsx_bulk, "mx_rhs")
      fl%mx_rhs(:, :, :) = fl%mx_rhs(:, :, :) - rhsx_bulk
    else if (fl%idriven == IDRVF_X_Cf) then
      rhsx_bulk = - HALF * fl%drvfc * dm%tAlpha(isub) * dm%dt
      fl%mx_rhs(:, :, :) = fl%mx_rhs(:, :, :) - rhsx_bulk
    else if (fl%idriven == IDRVF_Z_MASSFLUX) then
      call Get_volumetric_average_3d_for_var_xcx(dm, dm%dccp, fl%mz_rhs, rhsz_bulk, "mz_rhs")
      !call Get_volumetric_average_3d(.false., dm%ibcy(:, 3), dm%fbcy_var(:, :, :, 3), dm, dm%dccp, fl%mz_rhs, rhsz_bulk, "mz_rhs")
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

    real(WP), dimension( 4,              dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: fbcx_4cc
    real(WP), dimension( dm%dccc%ysz(1), 4,              dm%dccc%ysz(3) ) :: fbcy_c4c
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), 4              ) :: fbcz_cc4
    
#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_mid_msg("Updating the velocity/mass flux ...")
#endif
!----------------------------------------------------------------------------------------------------------
!   x-pencil, ux = ux - dt * alpha * d(phi_ccc)/dx
!----------------------------------------------------------------------------------------------------------
    dphidx_pcc = ZERO
    fbcx_4cc = ZERO ! check
    call Get_x_1st_derivative_C2P_3D(phi_ccc,  dphidx_pcc, dm, dm%iAccuracy, dm%ibcx(:, 4), fbcx_4cc )
    ux = ux - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidx_pcc
!----------------------------------------------------------------------------------------------------------
!   y-pencil, uy = uy - dt * alpha * d(phi_ccc)/dy
!----------------------------------------------------------------------------------------------------------
    phi_ccc_ypencil = ZERO
    dphidy_cpc_ypencil = ZERO
    dphidy_cpc = ZERO
    fbcy_c4c = ZERO
    call transpose_x_to_y (phi_ccc, phi_ccc_ypencil, dm%dccc)
    call Get_y_1st_derivative_C2P_3D(phi_ccc_ypencil, dphidy_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy(:, 4), fbcy_c4c )
    call transpose_y_to_x (dphidy_cpc_ypencil, dphidy_cpc, dm%dcpc)
    uy = uy - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidy_cpc
!----------------------------------------------------------------------------------------------------------
!   z-pencil, uz = uz - dt * alpha * d(phi_ccc)/dz
!----------------------------------------------------------------------------------------------------------
    pphi_ccc_zpencil =  ZERO
    dphidz_ccp_zpencil = ZERO
    dphidz_ccp_ypencil = ZERO
    dphidz_ccp = ZERO
    fbcz_cc4 = ZERO
    call transpose_y_to_z (phi_ccc_ypencil, pphi_ccc_zpencil, dm%dccc)
    call Get_z_1st_derivative_C2P_3D(pphi_ccc_zpencil, dphidz_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz(:, 4), fbcz_cc4 )
    call transpose_z_to_y (dphidz_ccp_zpencil, dphidz_ccp_ypencil, dm%dccp)
    call transpose_y_to_x (dphidz_ccp_ypencil, dphidz_ccp,         dm%dccp)
    uz = uz - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidz_ccp

    return
  end subroutine Correct_massflux

!==========================================================================================================
  subroutine solve_poisson(fl, dm, isub, tm)
    use udf_type_mod
    use parameters_constant_mod
    use decomp_2d_poisson
    use decomp_extended_mod
    use continuity_eq_mod

    implicit none
    type(t_domain), intent( in    ) :: dm
    type(t_flow),   intent( inout ) :: fl                  
    integer,        intent( in    ) :: isub
    type(t_thermo), optional, intent( in    ) :: tm

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: div
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: rhs_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: rhs_zpencil

  

    real(WP), dimension( dm%dccc%zst(1) : dm%dccc%zen(1), &
                         dm%dccc%zst(2) : dm%dccc%zen(2), &
                         dm%dccc%zst(3) : dm%dccc%zen(3) ) :: rhs_zpencil_ggg
    !integer :: i, j, k, jj, ii
    real(WP) :: coeff

#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_mid_msg("Calculating the RHS of Poisson Equation ...")
#endif

!==========================================================================================================
! RHS of Poisson Eq.
!==========================================================================================================
    fl%pcor = ZERO
!----------------------------------------------------------------------------------------------------------
! $d\rho / dt$ at cell centre
!----------------------------------------------------------------------------------------------------------
    if (dm%is_thermo) then
      call Calculate_drhodt(dm, tm%dDens, tm%dDensm1, tm%dDensm2, fl%pcor)
    end if
!----------------------------------------------------------------------------------------------------------
! $d(\rho u_i)) / dx_i $ at cell centre
!----------------------------------------------------------------------------------------------------------
    div  = ZERO
    call Get_divergence(fl, div, dm)
    coeff = ONE / (dm%tAlpha(isub) * dm%sigma2p * dm%dt)
    fl%pcor = fl%pcor + div
    fl%pcor = fl%pcor * coeff

#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug (fl%pcor, dm%dccc,   fl%iteration, isub, 'PhiRHS', '@RHS phi') ! debug_ww
    !call wrt_3d_all_debug(fl%pcor, dm%dccc,   fl%iteration, isub, 'PhiRHS', '@RHS phi') ! debug_ww
#endif
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
    !write(*, *) 'fft0', rhs_zpencil_ggg(:, 1, 1)
#endif
    call poisson(rhs_zpencil_ggg)
! #ifdef DEBUG_STEPS
!     write(*, *) 'fft1', rhs_zpencil_ggg(:, 1, 1)
! #endif
!==========================================================================================================
!   convert back RHS from zpencil ggg to xpencil gll
!==========================================================================================================
    call zpencil_index_ggg2llg(rhs_zpencil_ggg, rhs_zpencil, dm%dccc)
    call transpose_z_to_y (rhs_zpencil, rhs_ypencil, dm%dccc)
    call transpose_y_to_x (rhs_ypencil, fl%pcor,     dm%dccc)

#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug (fl%pcor, dm%dccc,   fl%iteration, isub, 'phi', '@sol phi') ! debug_ww
    !call wrt_3d_all_debug(fl%pcor, dm%dccc,   fl%iteration, isub, 'phi', '@sol phi') ! debug_ww
    !write(*,*) 'fft2', fl%pcor(:, 1, 1)
#endif
    return
  end subroutine
! !==========================================================================================================
!   subroutine solve_poisson_x2z(fl, dm, isub)
!     use udf_type_mod
!     use parameters_constant_mod
!     use decomp_2d_poisson
!     use decomp_extended_mod
!     use continuity_eq_mod

!     implicit none
!     type(t_domain), intent( in    ) :: dm
!     type(t_flow),   intent( inout ) :: fl                  
!     integer,        intent( in    ) :: isub

!     real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: rhs
!     real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: rhs_ypencil
!     real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: rhs_zpencil

  

!     real(WP), dimension( dm%dccc%zst(1) : dm%dccc%zen(1), &
!                          dm%dccc%zst(2) : dm%dccc%zen(2), &
!                          dm%dccc%zst(3) : dm%dccc%zen(3) ) :: rhs_zpencil_ggg
!     !integer :: i, j, k, jj, ii

! !==========================================================================================================
! ! RHS of Poisson Eq.
! !==========================================================================================================
!     fl%pcor_zpencil_ggg = ZERO
! !----------------------------------------------------------------------------------------------------------
! ! $d\rho / dt$ at cell centre
! !----------------------------------------------------------------------------------------------------------
!     if (dm%is_thermo) then
!       rhs = ZERO
!       rhs_ypencil = ZERO
!       rhs_zpencil = ZERO
!       rhs_zpencil_ggg = ZERO
!       call Calculate_drhodt(dm, tm%dDens, tm%dDensm1, tm%dDensm2, rhs)
!       call transpose_x_to_y(rhs,         rhs_ypencil)
!       call transpose_y_to_z(rhs_ypencil, rhs_zpencil)
!       call zpencil_index_llg2ggg(rhs_zpencil, rhs_zpencil_ggg, dm%dccc)
!       fl%pcor_zpencil_ggg = fl%pcor_zpencil_ggg + rhs_zpencil_ggg
!     end if
! !----------------------------------------------------------------------------------------------------------
! ! $d(\rho u_i)) / dx_i $ at cell centre
! !----------------------------------------------------------------------------------------------------------
!     rhs_zpencil_ggg  = ZERO
!     if (dm%is_thermo) then
!       call Get_divergence_x2z(fl%gx, fl%gy, fl%gz, rhs_zpencil_ggg, dm)
!     else
!       call Get_divergence_x2z(fl%qx, fl%qy, fl%qz, rhs_zpencil_ggg, dm)
!     end if
!     fl%pcor_zpencil_ggg = fl%pcor_zpencil_ggg + rhs_zpencil_ggg
!     fl%pcor_zpencil_ggg = fl%pcor_zpencil_ggg / (dm%tAlpha(isub) * dm%sigma2p * dm%dt)
! !==========================================================================================================
! !   solve Poisson
! !==========================================================================================================
!     call poisson(fl%pcor_zpencil_ggg)
! !==========================================================================================================
! !   convert back RHS from zpencil ggg to xpencil gll
! !==========================================================================================================
!     call zpencil_index_ggg2llg(fl%pcor_zpencil_ggg, rhs_zpencil, dm%dccc)
!     call transpose_z_to_y (rhs_zpencil, rhs_ypencil, dm%dccc)
!     call transpose_y_to_x (rhs_ypencil, fl%pcor,     dm%dccc)

!     return
!   end subroutine

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
  subroutine Solve_momentum_eq(dm, fl, isub, tm)
    use udf_type_mod
    use typeconvert_mod
    use continuity_eq_mod
    use boundary_conditions_mod
    use parameters_constant_mod
    use mpi_mod
    use solver_tools_mod
#ifdef DEBUG_STEPS
    use io_tools_mod
    use typeconvert_mod
    use wtformat_mod
#endif
    implicit none

    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(inout) :: dm
    type(t_thermo), intent(in), optional :: tm
    integer,        intent(in)    :: isub

!----------------------------------------------------------------------------------------------------------
! to calculate the rhs of the momenturn equation in stepping method
!----------------------------------------------------------------------------------------------------------
    call Compute_momentum_rhs(fl, dm, isub, tm)
!----------------------------------------------------------------------------------------------------------
! to update intermediate (\hat{q}) or (\hat{g})
!----------------------------------------------------------------------------------------------------------
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

#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%qx, dm%dpcc,   fl%iteration, isub, 'ux', '@bf divg') ! debug_ww
    call wrt_3d_pt_debug(fl%qy, dm%dcpc,   fl%iteration, isub, 'uy', '@bf divg') ! debug_ww
    call wrt_3d_pt_debug(fl%qz, dm%dccp,   fl%iteration, isub, 'uz', '@bf divg') ! debug_ww
    write(*,*) 'qx', fl%qx(:, 1, 1), fl%qx(:, 8, 8)
    write(*,*) 'qy', fl%qy(:, 1, 1), fl%qy(:, 8, 8)
    write(*,*) 'qz', fl%qz(:, 1, 1), fl%qz(:, 8, 8)
#endif

    !in order for a high order spacial accuracy
    ! to use Alternating direction implicit method
    ! ref: Cui2013: Convergence analysis of high-order compact 
    ! alternating direction implicit nnn/schemes for the two-dimensional 
    ! time fractional equation

! #ifdef DEBUG_STEPS
!   call write_snapshot_any3darray(fl%qx, 'qxs_RK'//trim(int2str(isub)), 'debug', dm%dpcc, dm, fl%iteration)
!   call write_snapshot_any3darray(fl%qy, 'qys_RK'//trim(int2str(isub)), 'debug', dm%dcpc, dm, fl%iteration)
!   call write_snapshot_any3darray(fl%qz, 'qzs_RK'//trim(int2str(isub)), 'debug', dm%dccp, dm, fl%iteration)
! #endif
!----------------------------------------------------------------------------------------------------------
! to solve Poisson equation
!----------------------------------------------------------------------------------------------------------
    !if(nrank == 0) call Print_debug_mid_msg("  Solving Poisson Equation ...") 
    !call solve_poisson_x2z(fl, dm, isub) !
    if ( .not. dm%is_thermo) then 
      call solve_poisson(fl, dm, isub) ! test show above two methods gave the same results. 
    else
      call solve_poisson(fl, dm, isub, tm) ! test show above two methods gave the same results. 
    end if
    
#ifdef DEBUG_STEPS
    call write_snapshot_any3darray(fl%pcor, 'pcor'//trim(int2str(isub)), 'debug', dm%dccc, dm, fl%iteration)
#endif
!----------------------------------------------------------------------------------------------------------
! to update pressure
!----------------------------------------------------------------------------------------------------------
#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_mid_msg("Correcting the pressure term ...")
#endif
    fl%pres(:, :, :) = fl%pres(:, :, :) + fl%pcor(:, :, :)
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%pres, dm%dccc,   fl%iteration, isub, 'pr', '@updated') ! debug_ww
#endif
!----------------------------------------------------------------------------------------------------------
! to update velocity/massflux correction
!----------------------------------------------------------------------------------------------------------
    !if(nrank == 0) call Print_debug_mid_msg("  Updating velocity/mass flux ...")
    if ( .not. dm%is_thermo) then 
      call Correct_massflux(fl%qx, fl%qy, fl%qz, fl%pcor, dm, isub)
    else
      call Correct_massflux(fl%gx, fl%gy, fl%gz, fl%pcor, dm, isub)
    end if
#ifdef DEBUG_STEPS
  call wrt_3d_pt_debug(fl%qx, dm%dpcc,   fl%iteration, isub, 'ux', '@updated') ! debug_ww
  call wrt_3d_pt_debug(fl%qy, dm%dcpc,   fl%iteration, isub, 'uy', '@updated') ! debug_ww
  call wrt_3d_pt_debug(fl%qz, dm%dccp,   fl%iteration, isub, 'uz', '@updated') ! debug_ww
#endif
!----------------------------------------------------------------------------------------------------------
! to update velocity from gx gy gz 
!----------------------------------------------------------------------------------------------------------
  if(dm%is_thermo) then
    call Calculate_velocity_from_massflux(dm, fl, tm)
    call apply_gxgygz_bc_geo(dm, fl, tm)
  end if

    
#ifdef DEBUG_STEPS
    if(dm%is_thermo) then
      call Find_maximum_velocity(dm, fl%gx, fl%gy, fl%gz)
      call Check_mass_conservation(dm, fl, "isub"//trim(int2str(isub)), tm) 
    else
      call Find_maximum_velocity(dm, fl%qx, fl%qy, fl%qz)
      call Check_mass_conservation(dm, fl, "isub"//trim(int2str(isub))) 
    end if
#endif

    return
  end subroutine Solve_momentum_eq

end module eq_momentum_mod