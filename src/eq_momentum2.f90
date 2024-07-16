module eq_momentum_mod
  use operations
  use precision_mod
  use decomp_2d
  use print_msg_mod
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
          rhs_explicit_current = rhs1(i, j, k) ! not (* dt)
          rhs_explicit_last    = rhs0(i, j, k) ! not (* dt)
          rhs_total = dm%tGamma(isub) * rhs_explicit_current + &
                      dm%tZeta (isub) * rhs_explicit_last
          rhs0(i, j, k) = rhs_explicit_current ! not (* dt)

      ! add pressure gradient
          rhs_total = rhs_total + &
                      dm%tAlpha(isub) * rhs1_pfc(i, j, k)

      ! times the time step 
          rhs1(i, j, k) = dm%dt * rhs_total ! * dt 

        end do
      end do
    end do

    return
  end subroutine
!==========================================================================================================
!> \brief To calcuate all rhs of momentum eq.
! The governing equation x is 
! d(gx)/dt = -        d(gxix * qxix)/dx                               ! conv-x-m1         
!            - 1/r  * d(gyix * qxiy)/dy                               ! conv-y-m1       
!            - 1/r2 * d(gzix * qxiz)/dz                               ! conv-z-m1        
!            - dpdx                                                   ! p-m1
!            + d[ 2 * mu * (qxdx - 1/3 * div)]/dx                     ! diff-x-m1
!            + 1/r  * d[muixy * (qydx + r * qxdy)]/dy                  ! diff-y-m1
!            + 1/r2 * d[muixz * (qzdx +     qxdz)]/dz                  ! diff-z-m1
! d(gy)/dt = -        d(gxiy * qyix)/dx                               ! conv-x-m2         
!            -        d(gyiy * qyriy)/dy                              ! conv-y-m2       
!            -        d(gzriy * qyriz)/dz                             ! conv-z-m2  
!            + (gzriz * qzriz)^y                                      ! conv-r-m2      
!            - dpdy                                                   ! p-m2
!            + d[muixy * (qydx + r * qxdy)]/dx                        ! diff-x-m2
!            + d[r * 2 * mu * (qyrdy - 1/3 * div)]/dy                 ! diff-y-m2
!            + d[muiyz * (qzrdy + 1/r * qyrdz - 1/r * qzriy]/dz       ! diff-z-m2
!            - [2 mu * (1/r2 * qzdz - 1/3 * div + 1/r * qyriy)]^y     ! diff-r-m2
! d(gz)/dt = -        d(gxiz * qzix)/dx                               ! conv-x-m3         
!            -        d(gyiz * qzriy)/dy                              ! conv-y-m3       
!            - 1/r2 * d(gziz * qziz)/dz                               ! conv-z-m3  
!            - (gyriz * qzriy)^y                                      ! conv-r-m3      
!            - dpdz                                                   ! p-m3
!            + d[muixz * (qzdx + qxdz)]/dx                            ! diff-x-m3
!            + d[muiyz * (r * qzrdy + qyrdz - qzriy]/dy               ! diff-y-m3
!            + d[2 * mu * (1/r2 * qzdz - 1/3 * div + 1/r * qyriy)]/dz ! diff-z-m3
!            + [muiyz * (qzrdy + 1/r2 * qydz - 1/r * qzriy)]^y        ! diff-r-m3
!______________________________________________________________________________
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

    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: my_rhs_ypencil
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: my_rhs_zpencil
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: my_rhs_pfc_xpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: my_rhs_pfc_ypencil

    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: mz_rhs_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: mz_rhs_zpencil
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: mz_rhs_pfc_xpencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: mz_rhs_pfc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: mz_rhs_pfc_zpencil
!----------------------------------------------------------------------------------------------------------
! intermediate variables
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: pres_ypencil     ! p-m2, common
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: pres_zpencil     ! p-m3, common

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: div_ccc_xpencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: div_ccc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: div_ccc_zpencil  
!----------------------------------------------------------------------------------------------------------
! qx, gx
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: qxix_ccc_xpencil ! conv-x-m1, common
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: qxiy_ppc_ypencil ! conv-y-m1, common
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: qxiy_ppc_xpencil ! conv-x-m2, no-thermo
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: qxiz_pcp_zpencil ! conv-z-m1, common
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: qxiz_pcp_xpencil ! conv-x-m3, no-thermo

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: gxix_ccc_xpencil ! conv-x-m1, thermo
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: gxiy_ppc_xpencil ! conv-x-m2, <=> qxiy_ppc_xpencil
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: gxiz_pcp_xpencil ! conv-x-m3, <=> qxiz_pcp_xpencil
    
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: qxdx_pcc_xpencil ! for div b.c., common
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: qxdx_cpc_ypencil ! for div b.c., common
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: qxdx_ccp_zpencil ! for div b.c., dommon
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: qxdx_ccc_xpencil ! diff-x-m1, common
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: qxdy_ppc_xpencil ! diff-x-m2, common
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: qxdy_ppc_ypencil ! diff-y-m1, common
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: qxdz_pcp_xpencil ! diff-x-m3, common
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: qxdz_pcp_zpencil ! diff-z-m1, common
!----------------------------------------------------------------------------------------------------------
! qy, gy
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: qyix_ppc_xpencil  ! conv-x-m2, common
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: qyix_ppc_ypencil  ! conv-y-m1, no-thermo
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qyiy_ccc_ypencil  ! conv-y-m2, no-thermo
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qyiz_cpp_ypencil  ! conv-y-m3, no-thermo
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qyiz_cpp_zpencil  ! conv-z-m2, no-cly1

    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: gyix_ppc_ypencil  ! conv-y-m1, <=> qyix_ppc_ypencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: gyiy_ccc_ypencil  ! conv-y-m2, <=> qyiy_ccc_ypencil
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: gyiz_cpp_ypencil  ! conv-y-m3, <=> qyiz_cpp_ypencil
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: gyriz_cpp_ypencil ! conv-r-m3, <=> qyriz_cpp_ypencil

    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: qyriy_ccc_zpencil ! diff-z-m3, cly
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qyriy_ccc_ypencil ! conv-y-m2, cly
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qyriz_cpp_ypencil ! conv-r-m3, no-thermal
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qyriz_cpp_zpencil ! conv-z-m2, cly1
    
    

    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: qydx_ppc_xpencil  ! diff-x-m2, common
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: qydx_ppc_ypencil  ! diff-y-m1, common
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qydz_cpp_ypencil  ! diff-r-m3, cly
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qydz_cpp_zpencil  ! diff-y-m3, no-cly2
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qydy_ccc_ypencil  ! diff-y-m2, no-cly1, div
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: qydy_pcc_xpencil  ! for div b.c., common
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: qydy_cpc_ypencil  ! for div b.c., common
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: qydy_ccp_zpencil  ! for div b.c., common

    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qyrdy_ccc_ypencil ! diff-y-m2, cly1
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qyrdz_cpp_zpencil ! diff-z-m2, cly2
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qyrdz_cpp_ypencil ! diff-y-m3, cly3
!----------------------------------------------------------------------------------------------------------
! qz, gz
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: qzix_pcp_xpencil  ! conv-x-m3, common
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: qzix_pcp_zpencil  ! conv-z-m1, no-thermo
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qziy_cpp_ypencil  ! conv-y-m3, no-cly
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qziy_cpp_zpencil  ! conv-z-m2, no-thermo, no-cly
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: qziz_ccc_zpencil  ! conv-z-m3, common
    
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: gzix_pcp_zpencil  ! conv-z-m1, <=> qzix_pcp_zpencil
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: gziy_cpp_zpencil  ! conv-z-m2
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: gziz_ccc_zpencil  ! conv-z-m3, <=> qziz_ccc_zpencil

    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qzriy_cpp_ypencil ! conv-y-m3, conv-r-m3, diff-y-m3, diff-r-m3, common
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qzriy_cpp_zpencil ! diff-z-m2, cly
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qzriz_ccc_ypencil ! conv-r-m2, cly
    
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: gzriy_cpp_zpencil ! conv-z-m2, 
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: gzriz_ccc_ypencil ! conv-r-m2, thermo+cly, 
    
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: qzdx_pcp_xpencil  ! diff-x-m3, common
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: qzdx_pcp_zpencil  ! diff-z-m1, common
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qzdy_cpp_ypencil  ! diff-y-m3, no-cly
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qzdy_cpp_zpencil  ! diff-z-m2, no-cly
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qzdz_ccc_ypencil  ! diff-r-m2, cly
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: qzdz_ccc_zpencil  ! diff-z-m3, common
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: qzdz_pcc_xpencil  ! for div b.c.
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: qzdz_cpc_ypencil  ! for div b.c.
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: qzdz_ccp_zpencil  ! for div b.c.

    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qzrdy_cpp_ypencil ! diff-y-m3, diff-r-m3
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qzrdy_cpp_zpencil ! diff-z-m2, cly
!----------------------------------------------------------------------------------------------------------
!   intermediate variables : if dm%is_thermo = true
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: mu_ccc_xpencil    ! x-v1, thermal only
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: mu_ccc_ypencil    ! y-v2, thermal only
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: mu_ccc_zpencil    ! z-v3, thermal only
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: muixy_ppc_xpencil ! diff-x-m2, thermal only
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: muixy_ppc_ypencil ! diff-y-m1, thermal only
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: muixz_pcp_xpencil ! diff-x-m3, thermal only
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: muixz_pcp_zpencil ! diff-z-m1, thermal only
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: muiyz_cpp_ypencil ! z-v2, z-v4, thermal only
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: muiyz_cpp_zpencil ! diff-z-m2, thermal only
!----------------------------------------------------------------------------------------------------------
!   intermediate variables : if fl%igravity == i
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: dDens_ypencil  ! d 
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: dDens_zpencil  ! d  
!----------------------------------------------------------------------------------------------------------
! temporary variables
!----------------------------------------------------------------------------------------------------------
#ifdef DEBUG_STEPS  
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: apcc_test
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: accp_test
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: acpc_test
#endif 
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: apcc_xpencil
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc_xpencil
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: appc_xpencil
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: appc_xpencil1
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: apcp_xpencil
    
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: appc_ypencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil1
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
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil1
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: acpc_zpencil
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: apcc_zpencil

    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: accp_xpencil
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: acpc_xpencil

!----------------------------------------------------------------------------------------------------------
! bc variables
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( 4, dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: fbcx_4cc 
    real(WP), dimension( 4, dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: fbcx_qx_4cc  
    real(WP), dimension( 4, dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: fbcx_qxdx_4cc
    real(WP), dimension( 4, dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: fbcx_gx_4cc   ! thermo
    real(WP), dimension( 4, dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: fbcx_div_4cc  ! common
    real(WP), dimension( 4, dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: fbcx_mu_4cc   ! thermo

    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_c4c
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_c4c1
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_qy_c4c
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_qyr_c4c
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_qydy_c4c
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_qyrdy_c4c
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_qzdz_c4c
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_qzr_c4c
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_gzr_c4c
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_gy_c4c
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_div_c4c
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_mu_c4c

    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), 4 ) :: fbcz_cc4
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), 4 ) :: fbcz_cc41
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), 4 ) :: fbcz_qyr_cc4
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), 4 ) :: fbcz_qz_cc4
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), 4 ) :: fbcz_qzdz_cc4
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), 4 ) :: fbcz_gz_cc4
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), 4 ) :: fbcz_div_cc4
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), 4 ) :: fbcz_mu_cc4

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
!               | -[ipy]-> qxiy_ppc_ypencil(common) -[y2x]-> qxiy_ppc_xpencil(no thermal)-[1dx]-> qxdx_cpc_xpencil(temp) -[x2y]-> qxdx_cpc_ypencil (for div bc)
!               | -[1dy]-> qxdy_ppc_ypencil(common) -[y2x]-> qxdy_ppc_xpencil(common)
!               | -[y2z]-> qx_zpencil(temp) 
!                          | -[ipz]-> qxiz_pcp_zpencil(common) -[z2x]-> qxiz_pcp_xpencil(no-thermal) 
!                                                            | -[1dx]-> qxdx_ccp_xpencil(temp) -[x2z]-> qxdx_ccp_zpencil (for div bc)
!                          | -[1dz]-> qxdz_pcp_zpencil(common) -[z2y]-> qxdz_pcp_ypencil(temp) -[y2x]-> qxdz_pcp_xpencil(common)
!----------------------------------------------------------------------------------------------------------
    call Get_x_midp_P2C_3D(fl%qx, qxix_ccc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_qx)
    call Get_x_1der_P2C_3D(fl%qx, qxdx_ccc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_qx)
    
    call transpose_x_to_y (fl%qx, apcc_ypencil, dm%dpcc) !apcc_ypencil = qx_ypencil
    call Get_y_1der_C2P_3D(apcc_ypencil, qxdy_ppc_ypencil, dm, dm%iAccuracy, dm%ibcy_qx, dm%fbcy_qx)
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

    call Get_z_1der_C2P_3D(apcc_zpencil, qxdz_pcp_zpencil, dm, dm%iAccuracy, dm%ibcz_qx, dm%fbcz_qx)
    call transpose_z_to_y(qxdz_pcp_zpencil, apcp_ypencil,     dm%dpcp) !qxdz_pcp_ypencil
    call transpose_y_to_x(apcp_ypencil,     qxdz_pcp_xpencil, dm%dpcp)
!----------------------------------------------------------------------------------------------------------
!   qx related b.c.
!----------------------------------------------------------------------------------------------------------
    if(is_fbcx_velo_required) then
      call Get_x_1der_C2P_3D(qxix_ccc_xpencil, qxdx_pcc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_qx)
      fbcx_qxdx_4cc(1, :, :) = qxdx_pcc_xpencil(1, :, :)
      fbcx_qxdx_4cc(2, :, :) = qxdx_pcc_xpencil(dm%dpcc%xsz(1), :, :)
      fbcx_qxdx_4cc(3:4, :, :) = fbcx_qxdx_4cc(1:2, :, :)

      fbcx_qx_4cc(1:2, :, :) = dm%fbcx_qx(1:2, :, :)
      fbcx_qx_4cc(3:4, :, :) = fbcx_qx_4cc(1:2, :, :)
    end if

    if(is_fbcy_velo_required) then
      ! qxiy_ppc_ypencil --> qxdx_cpc_ypencil (for div. b.c.)
      call transpose_y_to_x (qxiy_ppc_ypencil, appc_xpencil, dm%dppc)
      call Get_x_1der_P2C_3D(appc_xpencil, acpc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx)
      call transpose_x_to_y(acpc_xpencil, qxdx_cpc_ypencil, dm%dcpc)
    end if

    if(is_fbcz_velo_required) then
      ! qxiz_pcp_zpencil --> qxdx_ccp_zpencil (for div. b.c.)
      call transpose_z_to_y(qxiz_pcp_zpencil, apcp_ypencil, dm%dpcp)
      call transpose_y_to_x(apcp_ypencil,     apcp_xpencil, dm%dpcp)
      call Get_x_1der_P2C_3D(apcp_xpencil, accp_xpencil, dm, dm%iAccuracy, dm%ibcx_qx)
      call transpose_x_to_y(accp_xpencil, accp_ypencil,     dm%dccp)
      call transpose_y_to_z(accp_ypencil, qxdx_ccp_zpencil, dm%dccp) 
    end if

!----------------------------------------------------------------------------------------------------------
!    qy
!    | -[ipx]-> qyix_ppc_xpencil(common) -[x2y]-> qyix_ppc_ypencil(no thermal) 
!    | -[1dx]-> qydx_ppc_xpencil(common) -[x2y]-> qydx_ppc_ypencil(common)
!    | -[x2y]-> qy_ypencil(temp) 
!               | -[ipy]-> qyiy_ccc_ypencil(no thermal)
!               | -[1dy]-> qydy_ccc_ypencil(no-cly, div) 
!               | -[y2z]-> qy_zpencil(temp) 
!                          | -[ipz]-> qyiz_cpp_zpencil(no-cly) -[z2y]-> qyiz_cpp_ypencil(no-thermal, cly)    
!                          | -[1dz]-> qydz_cpp_zpencil(no-cly) -[z2y]-> qydz_cpp_ypencil(cly)
!----------------------------------------------------------------------------------------------------------
    call Get_x_1der_C2P_3D(fl%qy, qydx_ppc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_qy)
    call transpose_x_to_y(qydx_ppc_xpencil, qydx_ppc_ypencil, dm%dppc)
    call Get_x_midp_C2P_3D(fl%qy, qyix_ppc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_qy)
    if(.not. dm%is_thermo) then
      call transpose_x_to_y (qyix_ppc_xpencil, qyix_ppc_ypencil, dm%dppc)
    end if
    call transpose_x_to_y (fl%qy, acpc_ypencil, dm%dcpc) !acpc_ypencil = qy_ypencil
    call Get_y_1der_P2C_3D(acpc_ypencil, qydy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy)

    if(.not. dm%is_thermo) then
      call Get_y_midp_P2C_3D(acpc_ypencil, qyiy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)
    end if

    call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc)!acpc_zpencil = qy_zpencil
    call Get_z_1der_C2P_3D(acpc_zpencil, acpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qy) ! acpp_zpencil = qydz
    if(dm%icoordinate == ICARTESIAN) then
      qydz_cpp_zpencil = acpp_zpencil
    end if
    if(dm%icoordinate == ICYLINDRICAL) then
      call transpose_z_to_y(acpp_zpencil, qydz_cpp_ypencil, dm%dcpp)
    end if
    call Get_z_midp_C2P_3D(acpc_zpencil, acpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qy) ! acpp_zpencil = qyiz
    if(dm%icoordinate == ICARTESIAN) then
      qyiz_cpp_zpencil = acpp_zpencil
    end if
    if(.not. dm%is_thermo) then
      call transpose_z_to_y(acpp_zpencil, qyiz_cpp_ypencil, dm%dcpp)
    end if
!----------------------------------------------------------------------------------------------------------
!   qy related b.c.
!----------------------------------------------------------------------------------------------------------
    if(is_fbcx_velo_required) then
      !qyix_ppc_xpencil --> qydy_pcc_xpencil  for div b.c.
      call transpose_x_to_y (qyix_ppc_xpencil, appc_ypencil, dm%dppc)
      call Get_y_1der_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy)
      call transpose_y_to_x(apcc_ypencil, qydy_pcc_xpencil, dm%dpcc) 
    end if

    if(is_fbcy_velo_required) then
      !qydy_cpc_ypencil for div b.c.
      !acpc_ypencil = qy_ypencil
      call Get_y_midp_P2C_3D(acpc_ypencil, accc_ypencil,     dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)
      call Get_y_1der_C2P_3D(accc_ypencil, qydy_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)
      if(dm%icoordinate == ICARTESIAN) then
        fbcy_qy_c4c(:, :, :) = dm%fbcy_qy(:, :, :)
        fbcy_qydy_c4c(:, 1, :) = qydy_cpc_ypencil(:, 1,              :)
        fbcy_qydy_c4c(:, 2, :) = qydy_cpc_ypencil(:, dm%dcpc%ysz(2), :)
        fbcy_qydy_c4c(:, 3:4, :) = fbcy_qydy_c4c(:, 1:2, :)
      else if (dm%icoordinate == ICYLINDRICAL) then
        fbcy_qyr_c4c(:, :, :) = dm%fbcy_qyr(:, :, :)
      else
      end if
    end if

    if(is_fbcz_velo_required) then
      !qyiz_cpp_zpencil --> qydy_ccp_zpencil  for div b.c.
      call transpose_z_to_y(acpp_zpencil, acpp_ypencil, dm%dcpp)
      call Get_y_1der_P2C_3D(acpp_ypencil, accp_ypencil, dm, dm%iAccuracy, dm%ibcy_qy)
      call transpose_y_to_z(accp_ypencil, qydy_ccp_zpencil, dm%dpcc) 
    end if
!----------------------------------------------------------------------------------------------------------
!    qz 
!    | -[1dx]-> qzdx_pcp_xpencil(common) -[x2y]-> qzdx_pcp_ypencil(temp) -[y2z]-> qzdx_pcp_zpencil(common)
!    | -[ipx]-> qzix_pcp_xpencil(common) -[x2y]-> qzix_pcp_ypencil(temp) -[y2z]-> qzix_pcp_zpencil(no-thermal)
!    | -[x2y]-> qz_ypencil(temp) 
!               | -[1dy]-> qzdy_cpp_ypencil(no-cly) -[y2z]-> qzdy_cpp_zpencil(no-cly)
!               | -[ipy]-> qziy_cpp_ypencil(no-cly) -[y2z]-> qziy_cpp_zpencil(no cly, no thermal)
!               | -[y2z]-> qz_zpencil(temp) 
!                          | -[ipz]-> qziz_ccc_zpencil(common)
!                          | -[1dz]-> qzdz_ccc_zpencil(common) -[y2z]-> qzdz_ccc_ypencil(cly)
!----------------------------------------------------------------------------------------------------------
    call Get_x_1der_C2P_3D(fl%qz, qzdx_pcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_qz)
    call transpose_x_to_y(qzdx_pcp_xpencil, apcp_ypencil, dm%dpcp)
    call transpose_y_to_z(apcp_ypencil, qzdx_pcp_zpencil, dm%dpcp)

    call Get_x_midp_C2P_3D(fl%qz, qzix_pcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_qz)
    if(.not. dm%is_thermo) then
      call transpose_x_to_y(qzix_pcp_xpencil, apcp_ypencil, dm%dpcp) !qzix_pcp_ypencil
      call transpose_y_to_z(apcp_ypencil, qzix_pcp_zpencil, dm%dpcp)
    end if
    call transpose_x_to_y(fl%qz, accp_ypencil, dm%dccp) ! qz_ypencil
    if(dm%icoordinate == ICARTESIAN)then
      call Get_y_midp_C2P_3D(accp_ypencil, qziy_cpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qz)
      call Get_y_1der_C2P_3D(accp_ypencil, qzdy_cpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qz)
      call transpose_y_to_z(qzdy_cpp_ypencil, qzdy_cpp_zpencil, dm%dcpp)
    end if

    if((.not. dm%is_thermo) .and. (dm%icoordinate == ICARTESIAN)) then
        call transpose_y_to_z(qziy_cpp_ypencil, qziy_cpp_zpencil, dm%dcpp)
    end if

    call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp) ! qz_zpencil
    call Get_z_midp_P2C_3D(accp_zpencil, qziz_ccc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qz)
    call Get_z_1der_P2C_3D(accp_zpencil, qzdz_ccc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)

    if(dm%icoordinate == ICYLINDRICAL) then
      call transpose_z_to_y(qzdz_ccc_zpencil, qzdz_ccc_ypencil, dm%dccc)
    end if
!----------------------------------------------------------------------------------------------------------
!   qz related b.c.
!----------------------------------------------------------------------------------------------------------
    if(is_fbcx_velo_required) then
      call Get_z_1der_P2C_3D(qzix_pcp_zpencil, apcc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)
      call transpose_z_to_y(apcc_zpencil, apcc_ypencil, dm%dpcc )
      call transpose_y_to_x(apcc_ypencil, qzdz_pcc_xpencil, dm%dpcc)
    end if

    if(is_fbcy_velo_required) then
      call Get_y_midp_C2P_3D(accp_ypencil, acpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qz) !acpp_ypencil = qziy
      call transpose_y_to_z(acpp_ypencil, acpp_zpencil, dm%dcpp)
      call Get_z_1der_P2C_3D(acpp_zpencil, acpc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz) !acpc_zpencil = qzdz
      call transpose_z_to_y(acpc_zpencil, qzdz_cpc_ypencil, dm%dcpc)
      
      if(dm%icoordinate == ICYLINDRICAL) then
        call transpose_z_to_y(qziz_ccc_zpencil, accc_ypencil, dm%dccc)
        call Get_y_midp_C2P_3D(accc_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcy_qz)
        fbcy_qzdz_c4c(:, 1, :) = acpc_ypencil(:, 1, :)
        fbcy_qzdz_c4c(:, 2, :) = acpc_ypencil(:, dm%dcpc%ysz(2), :)
        fbcy_qzdz_c4c(:, 3:4, :) = fbcy_qzdz_c4c(:, 1:2, :)
      end if
    end if

    if(is_fbcz_velo_required) then
      fbcz_qz_cc4(:, :, :) = dm%fbcz_qz(:, :, :)
      call Get_z_1der_C2P_3D(qziz_ccc_zpencil, qzdz_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qz)
    end if
!==========================================================================================================
! preparation of div
! div = qxdx + 1/r * qydy + 1/r^2 qzdz
!==========================================================================================================
    ! div_tmp = qxdx
    div_ccc_xpencil  = ZERO 
    div_ccc_xpencil = div_ccc_xpencil + qxdx_ccc_xpencil !
    ! div_tmp = qxdx + 1/r * qydy
    accc_ypencil = qydy_ccc_ypencil
    if(dm%icoordinate == ICYLINDRICAL) then
      call multiple_cylindrical_rn(accc_ypencil, dm%dccc, dm%rci, 1, IPENCIL(2))
    end if
    call transpose_y_to_x (accc_ypencil, accc_xpencil, dm%dccc)   
    div_ccc_xpencil = div_ccc_xpencil + accc_xpencil
    ! div_tmp = qxdx + 1/r * qydy + 1/r^2 qzdz
    accc_zpencil = qzdz_ccc_zpencil
    if(dm%icoordinate == ICYLINDRICAL) then
      call multiple_cylindrical_rn(accc_zpencil, dm%dccc, dm%rci, 2, IPENCIL(3))
    end if
    call transpose_z_to_y (accc_zpencil, accc_ypencil, dm%dccc)   
    call transpose_y_to_x (accc_ypencil, accc_xpencil, dm%dccc)   
    div_ccc_xpencil = div_ccc_xpencil + accc_xpencil

    call transpose_x_to_y (div_ccc_xpencil, div_ccc_ypencil, dm%dccc)
    call transpose_y_to_z (div_ccc_ypencil, div_ccc_zpencil, dm%dccc)
    
    if(is_fbcx_velo_required) then
      fbcx_div_4cc(1, :, :) = qxdx_pcc_xpencil(1, :, :) + &
                              qydy_pcc_xpencil(1, :, :) + &
                              qzdz_pcc_xpencil(1, :, :)
      fbcx_div_4cc(2, :, :) = qxdx_pcc_xpencil(dm%dpcc%xsz(1), :, :) + &
                              qydy_pcc_xpencil(dm%dpcc%xsz(1), :, :) + &
                              qzdz_pcc_xpencil(dm%dpcc%xsz(1), :, :)
      fbcx_div_4cc(3:4, :, :) = fbcx_div_4cc(1:2, :, :)
    end if
    if(is_fbcy_velo_required) then
      fbcy_div_c4c(:, 1, :) = qxdx_cpc_ypencil(:, 1, :) + &
                              qydy_cpc_ypencil(:, 1, :) + &
                              qzdz_cpc_ypencil(:, 1, :)
      fbcy_div_c4c(:, 2, :) = qxdx_cpc_ypencil(:, dm%dcpc%ysz(2), :) + &
                              qydy_cpc_ypencil(:, dm%dcpc%ysz(2), :) + &
                              qzdz_cpc_ypencil(:, dm%dcpc%ysz(2), :)
      fbcy_div_c4c(:, 3:4, :) = fbcy_div_c4c(:, 1:2, :)
    end if
    if(is_fbcz_velo_required) then
      fbcz_div_cc4(:, :, 1) = qxdx_ccp_zpencil(:, :, 1) + &
                              qydy_ccp_zpencil(:, :, 1) + &
                              qzdz_ccp_zpencil(:, :, 1)
      fbcz_div_cc4(:, :, 2) = qxdx_ccp_zpencil(:, :, dm%dccp%zsz(3)) + &
                              qydy_ccp_zpencil(:, :, dm%dccp%zsz(3)) + &
                              qzdz_ccp_zpencil(:, :, dm%dccp%zsz(3))
      fbcz_div_cc4(:, :, 3:4) = fbcz_div_cc4(:, :, 1:2)
    end if
!==========================================================================================================
! preparation of intermediate variables to be used - thermal only
!==========================================================================================================
        mu_ccc_xpencil = fl%rre ! x-v1, default for flow only
        mu_ccc_ypencil = fl%rre ! y-v2, default for flow only
        mu_ccc_zpencil = fl%rre ! z-v3, default for flow only
    muixy_ppc_xpencil = fl%rre ! y-v1, default for flow only
    muixy_ppc_ypencil = fl%rre ! x-v2, default for flow only
    muixz_pcp_xpencil = fl%rre ! z-v1, default for flow only
    muixz_pcp_zpencil = fl%rre ! x-v3, default for flow only
    muiyz_cpp_ypencil = fl%rre ! z-v2, default for flow only
    muiyz_cpp_zpencil = fl%rre ! y-v3, default for flow only
          fbcx_mu_4cc = fl%rre
          fbcy_mu_c4c = fl%rre
          fbcz_mu_cc4 = fl%rre
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
      mu_ccc_xpencil = fl%mVisc
      call transpose_x_to_y(mu_ccc_xpencil, mu_ccc_ypencil, dm%dccc)
      call transpose_y_to_z(mu_ccc_ypencil, mu_ccc_zpencil, dm%dccc)

      call Get_x_midp_C2P_3D(mu_ccc_xpencil, apcc_xpencil, dm, dm%iAccuracy, dm%ibcx_Th, dm%fbcx_ftp%m) ! muix_pcc_xpencil
      call transpose_x_to_y(apcc_xpencil, apcc_ypencil, dm%dpcc) ! muix_pcc_ypencil
      
      call Get_y_midp_C2P_3D(apcc_ypencil, muixy_ppc_ypencil, dm, dm%iAccuracy, dm%ibcy_Th, dm%fbcy_ftp%m)
      call transpose_y_to_x(muixy_ppc_ypencil, muixy_ppc_xpencil, dm%dppc) 

      call transpose_y_to_z(apcc_ypencil, apcc_zpencil, dm%dpcc) ! muix_pcc_zpencil
      call Get_z_midp_C2P_3D(apcc_zpencil, muixz_pcp_zpencil, dm, dm%iAccuracy, dm%ibcz_Th, dm%fbcz_ftp%m)
      call transpose_z_to_y(muixz_pcp_zpencil, apcp_zpencil, dm%dpcp)
      call transpose_y_to_x(apcp_zpencil, muixz_pcp_xpencil, dm%dpcp)

      call Get_y_midp_C2P_3D(mu_ccc_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcy_Th, dm%fbcy_ftp%m) ! muiy_cpc_ypencil
      call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc) !muiy_cpc_zpencil
      call Get_z_midp_C2P_3D(acpc_zpencil, muiyz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_Th, dm%fbcz_ftp%m)
      call transpose_z_to_y(muiyz_cpp_zpencil, muiyz_cpp_ypencil, dm%dcpp)

      call Get_z_midp_C2P_3D(mu_ccc_zpencil, accp_zpencil, dm, dm%iAccuracy, dm%ibcz_Th, dm%fbcz_ftp%m)

      if(is_fbcx_velo_required) then
        fbcx_mu_4cc(1, :, :) = apcc_xpencil(1, :, :)
        fbcx_mu_4cc(2, :, :) = apcc_xpencil(dm%dpcc%xsz(1), :, :)
      end if

      if(is_fbcy_velo_required) then
        fbcy_mu_c4c(:, 1, :) = acpc_ypencil(:, 1, :)
        fbcy_mu_c4c(:, 2, :) = acpc_ypencil(:, dm%dcpc%ysz(2), :)
      end if

      if(is_fbcz_velo_required) then
        fbcz_mu_cc4(:, :, 1) = accp_zpencil(:, :, 1)
        fbcz_mu_cc4(:, :, 2) = accp_zpencil(:, :, dm%dccp%zsz(3))
      end if

    end if

!==========================================================================================================
! preparation of intermediate variables to be used - cylindrical only
!==========================================================================================================
    if(dm%icoordinate == ICYLINDRICAL) then
!----------------------------------------------------------------------------------------------------------
!    qy/r=qyr
!    | -[x2y]-> qyr_ypencil(temp) 
!               | -[ipy]-> qyriy_ccc_ypencil-[y2z]->qyriy_ccc_zpencil
!               | -[1dy]-> qyrdy_ccc_ypencil
!               | -[y2z]-> qyr_zpencil(temp) 
!                          | -[ipz]-> qyriz_cpp_zpencil -[z2y]-> qyriz_cpp_ypencil (no-thermal)        
!                          | -[1dz]-> qyrdz_cpp_zpencil -[z2y]-> qyrdz_cpp_ypencil
!----------------------------------------------------------------------------------------------------------
      acpc_xpencil = fl%qy
      call multiple_cylindrical_rn(acpc_xpencil, dm%dcpc, dm%rpi, 1, IPENCIL(1)) ! qr/r

      call transpose_x_to_y(acpc_xpencil, acpc_ypencil, dm%dcpc) ! acpc_ypencil = qyr_ypencil
      call Get_y_midp_P2C_3D(acpc_ypencil, qyriy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr)
      call transpose_y_to_z(qyriy_ccc_ypencil, qyriy_ccc_zpencil, dm%dccc)
      call Get_y_1der_P2C_3D(acpc_ypencil, qyrdy_ccc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr) ! check r treatment

      call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc)
      call Get_z_midp_C2P_3D(acpc_zpencil, qyriz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qyr)
      call Get_z_1der_C2P_3D(acpc_zpencil, qyrdz_cpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qyr)

      if(.not. dm%is_thermo) call transpose_z_to_y(qyriz_cpp_zpencil, qyriz_cpp_ypencil, dm%dcpp)
      call transpose_z_to_y(qyrdz_cpp_zpencil, qyrdz_cpp_ypencil, dm%dcpp)

      if(is_fbcy_velo_required) then
        fbcy_qyr_c4c(:, 1, :) = acpc_ypencil(:, 1,              :)
        fbcy_qyr_c4c(:, 2, :) = acpc_ypencil(:, dm%dcpc%ysz(2), :)

        call Get_y_1der_P2P_3D(acpc_ypencil, acpc_ypencil1, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qyr)
        fbcy_qyrdy_c4c(:, 1, :) = acpc_ypencil1(:, 1,              :)
        fbcy_qyrdy_c4c(:, 2, :) = acpc_ypencil1(:, dm%dcpc%ysz(2), :)
      end if
!----------------------------------------------------------------------------------------------------------
!    qz/r=qzr
!    | -[x2y]-> qzr_ypencil(temp) 
!               | -[ipy]-> qzriy_cpp_ypencil -[y2z]-> qzriy_cpp_zpencil
!               | -[1dy]-> qzrdy_cpp_ypencil -[y2z]-> qzrdy_cpp_zpencil
!               | -[y2z]-> qzr_zpencil(temp) -[ipz]-> qzriz_ccc_zpencil(temp) -[z2y]- qzriz_ccc_ypencil
!----------------------------------------------------------------------------------------------------------
      accp_xpencil = fl%qz
      call multiple_cylindrical_rn(accp_xpencil, dm%dccp, dm%rpi, 1, IPENCIL(1)) ! qz/r

      call transpose_x_to_y(accp_xpencil, accp_ypencil, dm%dccp)

      call Get_y_midp_C2P_3D(accp_ypencil, qzriy_cpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qzr)
      call transpose_y_to_z(qzriy_cpp_ypencil, qzriy_cpp_zpencil, dm%dcpp)

      call Get_y_1der_C2P_3D(accp_ypencil, qzrdy_cpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qzr)
      call transpose_y_to_z(qzrdy_cpp_ypencil, qzrdy_cpp_zpencil, dm%dcpp)

      call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
      call Get_z_midp_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qzr)
      call transpose_z_to_y(accc_zpencil, qzriz_ccc_ypencil, dm%dccc)

      if(is_fbcy_velo_required) then
        call Get_y_midp_C2P_3D(qzriz_ccc_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qzr)
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
        call Get_z_midp_C2P_3D(acpc_zpencil, acpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_gyr)
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

        call Get_y_midp_C2P_3D(accp_ypencil, acpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_gzr)
        call transpose_y_to_z(acpp_ypencil, gzriy_cpp_zpencil, dm%dccp)

        call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
        call Get_z_midp_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_gzr)
        call transpose_z_to_y(accc_zpencil, gzriz_ccc_ypencil, dm%dccc)

        if(is_fbcy_velo_required) then
          call Get_y_midp_C2P_3D(gzriz_ccc_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_gzr)
          fbcy_gzr_c4c(:, 1, :) = acpc_ypencil(:, 1,              :) 
          fbcy_gzr_c4c(:, 2, :) = acpc_ypencil(:, dm%dcpc%ysz(2), :) 
        end if
      end if

    end if
!==========================================================================================================
! The governing equation x is 
! d(gx)/dt = -        d(gxix * qxix)/dx                               ! conv-x-m1         
!            - 1/r  * d(gyix * qxiy)/dy                               ! conv-y-m1       
!            - 1/r2 * d(gzix * qxiz)/dz                               ! conv-z-m1        
!            - dpdx                                                   ! p-m1
!            + d[ 2 * mu * (qxdx - 1/3 * div)]/dx                     ! diff-x-m1
!            + 1/r  * d[muixy * (qydx + r * qxdy)]/dy                  ! diff-y-m1
!            + 1/r2 * d[muixz * (qzdx +     qxdz)]/dz                  ! diff-z-m1
!==========================================================================================================
    i = 1
    fl%mx_rhs      = ZERO
    mx_rhs_ypencil = ZERO
    mx_rhs_zpencil = ZERO
    mx_rhs_pfc_xpencil = ZERO
!----------------------------------------------------------------------------------------------------------
! X-mom convection term 1/3 at (i', j, k): 
! conv-x-m1 = - d(gxix * qxix)/dx
!----------------------------------------------------------------------------------------------------------  
    if ( .not. dm%is_thermo) then
      accc_xpencil = qxix_ccc_xpencil
      if(is_fbcx_velo_required) fbcx_4cc = -fbcx_qx_4cc * fbcx_qx_4cc
    else
      accc_xpencil = gxix_ccc_xpencil
      if(is_fbcx_velo_required) fbcx_4cc = -fbcx_gx_4cc * fbcx_qx_4cc
    end if
    accc_xpencil = - accc_xpencil * qxix_ccc_xpencil
    call Get_x_1der_C2P_3D(accc_xpencil, apcc_xpencil, dm, dm%iAccuracy, mbcx_cov1, fbcx_4cc)
    fl%mx_rhs = fl%mx_rhs + apcc_xpencil
!----------------------------------------------------------------------------------------------------------
! X-mom convection term 2/3 at (i', j, k): 
! conv-y-m1 = - 1/r  * d(gyix * qxiy)/dy
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      appc_ypencil = qyix_ppc_ypencil
    else
      appc_ypencil = gyix_ppc_ypencil
    end if
    appc_ypencil = -appc_ypencil * qxiy_ppc_ypencil
    call Get_y_1der_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%iAccuracy, mbcy_cov1)
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_ypencil, dm%dpcc, dm%rci, 1, IPENCIL(2))
    mx_rhs_ypencil = mx_rhs_ypencil + apcc_ypencil
!----------------------------------------------------------------------------------------------------------
! X-mom convection term 3/3 at (i', j, k) : 
! conv-z-m1 = - 1/r2 * d(gzix * qxiz)/dz
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      apcp_zpencil = qzix_pcp_zpencil
    else
      apcp_zpencil = gzix_pcp_zpencil
    end if
    apcp_zpencil = -apcp_zpencil * qxiz_pcp_zpencil
    call Get_z_1der_P2C_3D(apcp_zpencil, apcc_zpencil, dm, dm%iAccuracy, mbcz_cov1)
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_zpencil, dm%dpcc, dm%rci, 2, IPENCIL(3))
    mx_rhs_zpencil = mx_rhs_zpencil + apcc_zpencil

#ifdef DEBUG_STEPS
    apcc_test = fl%mx_rhs
    call transpose_y_to_x (mx_rhs_ypencil, apcc_xpencil, dm%dpcc)
    apcc_test =  apcc_test + apcc_xpencil
    call transpose_z_to_y (mx_rhs_zpencil, apcc_ypencil, dm%dpcc)
    call transpose_y_to_x (apcc_ypencil, apcc_xpencil, dm%dpcc)
    apcc_test =  apcc_test + apcc_xpencil 
    call wrt_3d_pt_debug(apcc_test,  dm%dpcc, fl%iteration, isub, 'ConX@bf st') ! debug_ww
#endif
!----------------------------------------------------------------------------------------------------------
! X-mom pressure gradient: 
! p-m1 = - dpdx
!----------------------------------------------------------------------------------------------------------
    call Get_x_1der_C2P_3D( -fl%pres, apcc_xpencil, dm, dm%iAccuracy, dm%ibcx_pr, -dm%fbcx_pr)
    mx_rhs_pfc_xpencil=  mx_rhs_pfc_xpencil + apcc_xpencil
!----------------------------------------------------------------------------------------------------------
! X-mom gravity in x direction, X-pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%is_thermo .and. (fl%igravity == i .or. fl%igravity == -i) )  then
      call Get_x_midp_C2P_3D(fl%dDens, apcc_xpencil, dm, dm%iAccuracy, dm%ibcx_Th, dm%fbcx_ftp%d )
      mx_rhs_pfc_xpencil =  mx_rhs_pfc_xpencil + fl%fgravity(i) * apcc_xpencil
    end if
!----------------------------------------------------------------------------------------------------------
! X-mom diffusion term 1/3 at (i', j, k)
! diff-x-m1 = d[ 2 * mu * (qxdx - 1/3 * div)]/dx
!----------------------------------------------------------------------------------------------------------
    accc_xpencil = TWO * mu_ccc_xpencil * ( qxdx_ccc_xpencil - ONE_THIRD * div_ccc_xpencil)
    if(is_fbcx_velo_required) then
      fbcx_4cc = TWO * fbcx_mu_4cc * ( fbcx_qxdx_4cc - ONE_THIRD * fbcx_div_4cc)
    end if
    call Get_x_1der_C2P_3D(accc_xpencil,  apcc_xpencil, dm, dm%iAccuracy, mbcx_tau1, fbcx_4cc) 
    fl%mx_rhs = fl%mx_rhs + apcc_xpencil
#ifdef DEBUG_STEPS
    write(*,*) 'visx-11', apcc_xpencil(1, 1:4, 1)
#endif
!----------------------------------------------------------------------------------------------------------
! X-mom diffusion term 2/3 at (i', j, k)
! diff-y-m1 = 1/r  * d[muixy * (qydx + r * qxdy)]/dy
!----------------------------------------------------------------------------------------------------------
    appc_ypencil = qxdy_ppc_ypencil
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(appc_ypencil, dm%dppc, ONE/dm%rpi, 1, IPENCIL(2))
    appc_ypencil = ( appc_ypencil + qydx_ppc_ypencil) * muixy_ppc_ypencil
    call Get_y_1der_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%iAccuracy, mbcy_tau1)
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_ypencil, dm%dpcc, dm%rci, 1, IPENCIL(2))
    mx_rhs_ypencil =  mx_rhs_ypencil + apcc_ypencil
#ifdef DEBUG_STEPS
    write(*,*) 'visx-12', apcc_ypencil(1, 1:4, 1)
#endif
!----------------------------------------------------------------------------------------------------------
! X-mom diffusion term 3/3 at (i', j, k)
! diff-z-m1 = 1/r2 * d[muixz * (qzdx +     qxdz)]/dz
!----------------------------------------------------------------------------------------------------------
    apcp_zpencil = ( qxdz_pcp_zpencil + qzdx_pcp_zpencil) * muixz_pcp_zpencil
    call Get_z_1der_P2C_3D(apcp_zpencil, apcc_zpencil, dm, dm%iAccuracy, mbcz_tau1)
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_zpencil, dm%dpcc, dm%rci, 2, IPENCIL(3))
    mx_rhs_zpencil =  mx_rhs_zpencil + apcc_zpencil
#ifdef DEBUG_STEPS
    write(*,*) 'visx-13', apcc_zpencil(1, 1:4, 1)
#endif
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
! d(gy)/dt = -        d(gxiy * qyix)/dx                               ! conv-x-m2         
!            -        d(gyiy * qyriy)/dy                              ! conv-y-m2       
!            -        d(gzriy * qyriz)/dz                             ! conv-z-m2  
!            + (gzriz * qzriz)^y                                      ! conv-r-m2      
!            - dpdy                                                   ! p-m2
!            + d[muixy * (qydx + r * qxdy)]/dx                        ! diff-x-m2
!            + d[r * 2 * mu * (qyrdy - 1/3 * div)]/dy                 ! diff-y-m2
!            + d[muiyz * (qzrdy + 1/r * qyrdz - 1/r * qzriy]/dz       ! diff-z-m2
!            - [2 mu * (1/r2 * qzdz - 1/3 * div + 1/r * qyriy)]^y     ! diff-r-m2
!==========================================================================================================
    i = 2
    fl%my_rhs          = ZERO
    my_rhs_ypencil     = ZERO
    my_rhs_zpencil     = ZERO
    my_rhs_pfc_xpencil = ZERO
    my_rhs_pfc_ypencil = ZERO
!----------------------------------------------------------------------------------------------------------
! Y-mom convection term 1/4 at (i, j', k)
! conv-x-m2 = - d(gxiy * qyix)/dx
!----------------------------------------------------------------------------------------------------------  
    if ( .not. dm%is_thermo) then
      appc_xpencil = qxiy_ppc_xpencil
    else
      appc_xpencil = gxiy_ppc_xpencil
    end if
    appc_xpencil = - appc_xpencil * qyix_ppc_xpencil
    call Get_x_1der_P2C_3D( appc_xpencil, acpc_xpencil, dm, dm%iAccuracy, mbcx_cov2)
    fl%my_rhs = fl%my_rhs + acpc_xpencil
!----------------------------------------------------------------------------------------------------------
! Y-mom convection term 2/4 at (i, j', k)
! conv-y-m2 = - d(gyiy * qyriy)/dy
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
    call Get_y_1der_C2P_3D( accc_ypencil, acpc_ypencil, dm, dm%iAccuracy, mbcy_cov2, fbcy_c4c)
    my_rhs_ypencil = my_rhs_ypencil + acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Y-mom convection term 3/4 at (i, j', k)
! conv-z-m2 = - d(gzriy * qyriz)/dz
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
    call Get_z_1der_P2C_3D( acpp_zpencil, acpc_zpencil, dm, dm%iAccuracy, mbcz_cov2)
    my_rhs_zpencil = my_rhs_zpencil + acpc_zpencil
!----------------------------------------------------------------------------------------------------------
! Y-mom convection term 4/4 at (i, j', k)
! conv-r-m2 = (gzriz * qzriz)^y
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
#ifdef DEBUG_STEPS
    acpc_test = fl%my_rhs
    call transpose_y_to_x (my_rhs_ypencil, acpc_xpencil, dm%dcpc)
    acpc_test =  acpc_test + acpc_xpencil
    call transpose_z_to_y (my_rhs_zpencil, acpc_ypencil, dm%dcpc)
    call transpose_y_to_x (acpc_ypencil, acpc_xpencil, dm%dcpc)
    acpc_test =  acpc_test + acpc_xpencil 
    call wrt_3d_pt_debug(acpc_test,  dm%dcpc, fl%iteration, isub, 'ConY@bf st') ! debug_ww
#endif
!----------------------------------------------------------------------------------------------------------
! Y-mom pressure gradient in y direction, Y-pencil, d(sigma_1 p)
! p-m2 = - dpdy
!----------------------------------------------------------------------------------------------------------
    call Get_y_1der_C2P_3D( -pres_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcy_pr, -dm%fbcy_pr)
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(acpc_ypencil, dm%dcpc, ONE/dm%rpi, 1, IPENCIL(2))
    my_rhs_pfc_ypencil =  my_rhs_pfc_ypencil + acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Y-mom gravity in y direction, Y-pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%is_thermo .and. (fl%igravity == i .or. fl%igravity == -i) )  then
      call Get_y_midp_C2P_3D(dDens_ypencil, acpc_ypencil, dm, dm%iAccuracy, dm%ibcy_Th, dm%fbcy_ftp%d )
      my_rhs_pfc_ypencil =  my_rhs_pfc_ypencil + fl%fgravity(i) * acpc_ypencil
    end if
!----------------------------------------------------------------------------------------------------------
! Y-mom diffusion term 1/4 at (i, j', k)
! diff-x-m2 = d[muixy * (qydx + r * qxdy)]/dx                        
!----------------------------------------------------------------------------------------------------------
    appc_xpencil1 = qxdy_ppc_xpencil
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(appc_xpencil1, dm%dppc, ONE/dm%rpi, 1, IPENCIL(1))
    appc_xpencil = ( appc_xpencil1 + qydx_ppc_xpencil) * muixy_ppc_xpencil
    call Get_x_1der_P2C_3D(appc_xpencil, acpc_xpencil, dm, dm%iAccuracy, mbcx_tau2)
    fl%my_rhs = fl%my_rhs + acpc_xpencil
#ifdef DEBUG_STEPS
    write(*,*) 'visy-21', acpc_xpencil(1, 1:4, 1)
#endif
!----------------------------------------------------------------------------------------------------------
! Y-mom diffusion term 2/4 at (i, j', k)
! diff-y-m2 = d[r * 2 * mu * (qyrdy - 1/3 * div)]/dy                 
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      accc_ypencil1 = qyrdy_ccc_ypencil
      if(is_fbcy_velo_required) then
        fbcy_c4c = ( fbcy_qyrdy_c4c - ONE_THIRD * fbcy_div_c4c) * TWO * fbcy_mu_c4c
        call multiple_cylindrical_rn_x4x(fbcy_c4c, dm%dcpc, ONE/dm%rci, 1, IPENCIL(2))
      end if
    else if(dm%icoordinate == ICARTESIAN) then
      accc_ypencil1 = qydy_ccc_ypencil
      if(is_fbcy_velo_required) then
        fbcy_c4c = ( fbcy_qydy_c4c - ONE_THIRD * fbcy_div_c4c) * TWO * fbcy_mu_c4c
      end if
    else
    end if
    
    accc_ypencil = ( accc_ypencil1 - ONE_THIRD * div_ccc_ypencil) * TWO * mu_ccc_ypencil    
    if(dm%icoordinate == ICYLINDRICAL) then
      call multiple_cylindrical_rn(accc_ypencil, dm%dccc, ONE/dm%rci, 1, IPENCIL(2))
    end if
    call Get_y_1der_C2P_3D(accc_ypencil,  acpc_ypencil, dm, dm%iAccuracy, mbcy_tau2, fbcy_c4c) 
    my_rhs_ypencil = my_rhs_ypencil + acpc_ypencil
#ifdef DEBUG_STEPS
    write(*,*) 'visy-22', acpc_ypencil(1, 1:4, 1)
#endif
!----------------------------------------------------------------------------------------------------------
! Y-mom diffusion term 3/4 at (i, j', k)
! diff-z-m2 = d[muiyz * (qzrdy + 1/r * qyrdz - 1/r * qzriy]/dz       
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
    call Get_z_1der_P2C_3D(acpp_zpencil, acpc_zpencil, dm, dm%iAccuracy, mbcz_tau2)
    my_rhs_zpencil =  my_rhs_zpencil + acpc_zpencil
#ifdef DEBUG_STEPS
    write(*,*) 'visy-23', acpc_zpencil(1, 1:4, 1)
#endif
!----------------------------------------------------------------------------------------------------------
! Y-mom diffusion term 4/4 at (i, j', k)
! diff-r-m2 = - [2 mu * (1/r2 * qzdz - 1/3 * div + 1/r * qyriy)]^y     
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
#ifdef DEBUG_STEPS
      write(*,*) 'visy-24', acpc_ypencil(1, 1:4, 1)
#endif
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
! d(gz)/dt = -        d(gxiz * qzix)/dx                               ! conv-x-m3         
!            -        d(gyiz * qzriy)/dy                              ! conv-y-m3       
!            - 1/r2 * d(gziz * qziz)/dz                               ! conv-z-m3  
!            - (gyriz * qzriy)^y                                      ! conv-r-m3      
!            - dpdz                                                   ! p-m3
!            + d[muixz * (qzdx + qxdz)]/dx                            ! diff-x-m3
!            + d[muiyz * (r * qzrdy + qyrdz - qzriy]/dy               ! diff-y-m3
!            + d[2 * mu * (1/r2 * qzdz - 1/3 * div + 1/r * qyriy)]/dz ! diff-z-m3
!            + [muiyz * (qzrdy + 1/r2 * qydz - 1/r * qzriy)]^y        ! diff-r-m3
!==========================================================================================================
    i = 3
    fl%mz_rhs          = ZERO
    mz_rhs_ypencil     = ZERO
    mz_rhs_zpencil     = ZERO
    mz_rhs_pfc_xpencil = ZERO
    mz_rhs_pfc_ypencil = ZERO
    mz_rhs_pfc_zpencil = ZERO
!----------------------------------------------------------------------------------------------------------
! Z-mom convection term 1/4 at (i, j, k')
! conv-x-m3 = - d(gxiz * qzix)/dx                               
!----------------------------------------------------------------------------------------------------------  
    if ( .not. dm%is_thermo) then
      apcp_xpencil = qxiz_pcp_xpencil
    else
      apcp_xpencil = gxiz_pcp_xpencil
    end if
    apcp_xpencil = - apcp_xpencil * qzix_pcp_xpencil
    call Get_x_1der_P2C_3D(apcp_xpencil, accp_xpencil, dm, dm%iAccuracy, mbcx_cov3)
    fl%mz_rhs = fl%mz_rhs + accp_xpencil
!----------------------------------------------------------------------------------------------------------
! Z-mom convection term 2/4 at (i, j, k')
! conv-y-m3 = - d(gyiz * qzriy)/dy                              
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
    call Get_y_1der_P2C_3D( acpp_ypencil, accp_ypencil, dm, dm%iAccuracy, mbcy_cov3)
    mz_rhs_ypencil = mz_rhs_ypencil + accp_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-mom convection term 3/4 at (i, j, k')
! conv-z-m3 = - 1/r2 * d(gziz * qziz)/dz                               
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      accc_zpencil = qziz_ccc_zpencil
      if(is_fbcz_velo_required) fbcz_cc4 = -fbcz_qz_cc4 * fbcz_qz_cc4
    else
      accc_zpencil = gziz_ccc_zpencil
      if(is_fbcz_velo_required) fbcz_cc4 = -fbcz_gz_cc4 * fbcz_qz_cc4
    end if
    accc_zpencil = -accc_zpencil * qziz_ccc_zpencil
    call Get_z_1der_C2P_3D( accc_zpencil, accp_zpencil, dm, dm%iAccuracy, mbcz_cov3, fbcz_cc4 )
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(accp_zpencil, dm%dccp, dm%rci, 2, IPENCIL(3))
    mz_rhs_zpencil = mz_rhs_zpencil + accp_zpencil
!----------------------------------------------------------------------------------------------------------
! Z-mom convection term 4/4 at (i, j, k')
! conv-r-m3 = - (gyriz * qzriy)^y                                      
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
#ifdef DEBUG_STEPS
    accp_test = fl%mz_rhs
    call transpose_y_to_x (mz_rhs_ypencil, accp_xpencil, dm%dccp)
    accp_test =  accp_test + accp_xpencil
    call transpose_z_to_y (mz_rhs_zpencil, accp_ypencil, dm%dccp)
    call transpose_y_to_x (accp_ypencil, accp_xpencil, dm%dccp)
    accp_test =  accp_test + accp_xpencil 
    call wrt_3d_pt_debug(accp_test,  dm%dccp, fl%iteration, isub, 'ConZ@bf st') ! debug_ww
#endif
!----------------------------------------------------------------------------------------------------------
! Z-mom pressure gradient in z direction, Z-pencil, d(sigma_1 p)
! p-m3 = -dpdz
!----------------------------------------------------------------------------------------------------------
    call Get_z_1der_C2P_3D( -pres_zpencil, accp_zpencil, dm, dm%iAccuracy, dm%ibcz_pr, -dm%fbcz_pr)
    mz_rhs_pfc_zpencil =  mz_rhs_pfc_zpencil + accp_zpencil
!----------------------------------------------------------------------------------------------------------
! Z-mom gravity in z direction, Z-pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%is_thermo .and. (fl%igravity == i .or. fl%igravity == -i) )  then
      call Get_z_midp_C2P_3D(dDens_zpencil, accp_zpencil, dm, dm%iAccuracy, dm%ibcz_Th, dm%fbcz_ftp%d )
      mz_rhs_pfc_zpencil =  mz_rhs_pfc_zpencil + fl%fgravity(i) * accp_zpencil
    end if
!----------------------------------------------------------------------------------------------------------
! Z-mom diffusion term 1/4  at (i, j, k')
! diff-x-m3 = d[muixz * (qzdx + qxdz)]/dx                            
!----------------------------------------------------------------------------------------------------------
    apcp_xpencil = ( qxdz_pcp_xpencil + qzdx_pcp_xpencil) * muixz_pcp_xpencil
    call Get_x_1der_P2C_3D(apcp_xpencil, accp_xpencil, dm, dm%iAccuracy, mbcx_tau3)
    fl%mz_rhs = fl%mz_rhs + accp_xpencil
#ifdef DEBUG_STEPS
      write(*,*) 'visz-31', accp_xpencil(1, 1:4, 1)
#endif
!----------------------------------------------------------------------------------------------------------
! Z-mom diffusion term 2/4 at (i, j, k')
! diff-y-m3 = d[muiyz * (r * qzrdy + qyrdz - qzriy]/dy               
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      acpp_ypencil = qzrdy_cpp_ypencil
      call multiple_cylindrical_rn(acpp_ypencil, dm%dcpp, ONE/dm%rpi, 1, IPENCIL(2))

      acpp_ypencil = acpp_ypencil + qyrdz_cpp_ypencil
      acpp_ypencil = acpp_ypencil - qzriy_cpp_ypencil
    else
      acpp_ypencil = qzdy_cpp_ypencil
      acpp_ypencil = acpp_ypencil + qydz_cpp_ypencil
    end if
    acpp_ypencil = acpp_ypencil * muiyz_cpp_ypencil
    call Get_y_1der_P2C_3D(acpp_ypencil, accp_ypencil, dm, dm%iAccuracy, mbcy_tau3)
    mz_rhs_ypencil =  mz_rhs_ypencil + accp_ypencil
#ifdef DEBUG_STEPS
    write(*,*) 'visz-32', accp_ypencil(1, 1:4, 1)
#endif
!----------------------------------------------------------------------------------------------------------
! Z-mom diffusion term 3/4 at (i, j, k')
! diff-z-m3 = d[2 * mu * (1/r2 * qzdz - 1/3 * div + 1/r * qyriy)]/dz 
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
    call Get_z_1der_C2P_3D(accc_zpencil, accp_zpencil, dm, dm%iAccuracy, mbcz_tau3, fbcz_cc4) 
    mz_rhs_zpencil = mz_rhs_zpencil + accp_zpencil
#ifdef DEBUG_STEPS
    write(*,*) 'visz-33', accp_zpencil(1, 1:4, 1)
#endif
!----------------------------------------------------------------------------------------------------------
! Z-mom diffusion term 4/4 at (i, j, k')
! diff-r-m3 = [muiyz * (qzrdy + 1/r2 * qydz - 1/r * qzriy)]^y        
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
#ifdef DEBUG_STEPS
      write(*,*) 'visz-34', accp_ypencil(1, 1:4, 1)
#endif
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
    call wrt_3d_pt_debug(fl%mx_rhs, dm%dpcc, fl%iteration, isub, 'ConVisX@bf stepping') ! debug_ww
    !write(*,*) 'fl%mx_rhs', fl%mx_rhs(:, 1, 1), fl%mx_rhs(:, 8, 8)
#endif
    call Calculate_momentum_fractional_step(fl%mx_rhs0, fl%mx_rhs, mx_rhs_pfc_xpencil, dm%dpcc, dm, isub)  
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%mx_rhs, dm%dpcc, fl%iteration, isub, 'ConVisPX@af stepping') ! debug_ww
#endif
!----------------------------------------------------------------------------------------------------------
! x-pencil : y-momentum
!----------------------------------------------------------------------------------------------------------
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%my_rhs, dm%dcpc, fl%iteration, isub, 'ConVisY@bf stepping') ! debug_ww
    !write(*,*) 'fl%my_rhs', fl%my_rhs(:, 1, 1), fl%my_rhs(:, 8, 8)
#endif
    call Calculate_momentum_fractional_step(fl%my_rhs0, fl%my_rhs, my_rhs_pfc_xpencil, dm%dcpc, dm, isub)
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%my_rhs, dm%dcpc, fl%iteration, isub, 'ConVisY@af stepping') ! debug_ww
#endif
!----------------------------------------------------------------------------------------------------------
! x-pencil : z-momentum
!----------------------------------------------------------------------------------------------------------
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%mz_rhs, dm%dccp, fl%iteration, isub, 'ConVisZ@bf stepping') ! debug_ww
    !write(*,*) 'fl%mz_rhs', fl%mz_rhs(:, 1, 1), fl%mz_rhs(:, 8, 8)
#endif
    call Calculate_momentum_fractional_step(fl%mz_rhs0, fl%mz_rhs, mz_rhs_pfc_xpencil, dm%dccp, dm, isub)
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%mz_rhs, dm%dccp, fl%iteration, isub, 'ConVisZ@af stepping') ! debug_ww
#endif
!==========================================================================================================
! x-pencil : flow drive terms (source terms) in periodic Streamwise flow
!==========================================================================================================
    if (fl%idriven == IDRVF_X_MASSFLUX) then
      call Get_volumetric_average_3d_for_var_xcx(dm, dm%dpcc, fl%mx_rhs, rhsx_bulk, "mx_rhs")
      !call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), dm%fbcy_var(:, :, :, 1), dm, dm%dpcc, fl%mx_rhs, rhsx_bulk, "mx_rhs")
      fl%mx_rhs = fl%mx_rhs - rhsx_bulk
    else if (fl%idriven == IDRVF_X_Cf) then
      rhsx_bulk = - HALF * fl%drvfc * dm%tAlpha(isub) * dm%dt
      fl%mx_rhs = fl%mx_rhs - rhsx_bulk
    else if (fl%idriven == IDRVF_Z_MASSFLUX) then
      call Get_volumetric_average_3d_for_var_xcx(dm, dm%dccp, fl%mz_rhs, rhsz_bulk, "mz_rhs")
      !call Get_volumetric_average_3d(.false., dm%ibcy(:, 3), dm%fbcy_var(:, :, :, 3), dm, dm%dccp, fl%mz_rhs, rhsz_bulk, "mz_rhs")
      fl%mz_rhs = fl%mz_rhs - rhsz_bulk
    else if (fl%idriven == IDRVF_Z_Cf) then
      rhsz_bulk = - HALF * fl%drvfc * dm%tAlpha(isub) * dm%dt
      fl%mz_rhs = fl%mz_rhs - rhsz_bulk
    else
    end if
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%mx_rhs, dm%dpcc, fl%iteration, isub, 'RHSX@total') ! debug_ww
    call wrt_3d_pt_debug(fl%mz_rhs, dm%dccp, fl%iteration, isub, 'RHSZ@total') ! debug_ww
#endif
 
    return
  end subroutine Compute_momentum_rhs

!==========================================================================================================
!==========================================================================================================
  subroutine Correct_massflux(fl, phi_ccc, dm, isub)
    use udf_type_mod
    use input_general_mod
    use operations
    use parameters_constant_mod
    implicit none
    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(in) :: dm
    integer,        intent(in) :: isub
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ), intent(in ) :: phi_ccc

    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: ux
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: uy
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: uz

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

    if(dm%is_thermo) then
      ux = fl%gx
      uy = fl%gy
      uz = fl%gz
    else 
      ux = fl%qx
      uy = fl%qy
      uz = fl%qz
    end if
!----------------------------------------------------------------------------------------------------------
!   x-pencil, ux = ux - dt * alpha * d(phi_ccc)/dx
!----------------------------------------------------------------------------------------------------------
    dphidx_pcc = ZERO
    fbcx_4cc = ZERO ! check
    call Get_x_1der_C2P_3D(phi_ccc,  dphidx_pcc, dm, dm%iAccuracy, dm%ibcx_pr, fbcx_4cc )
    ux = ux - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidx_pcc
!----------------------------------------------------------------------------------------------------------
!   y-pencil, uy = uy - dt * alpha * d(phi_ccc)/dy
!----------------------------------------------------------------------------------------------------------
    phi_ccc_ypencil = ZERO
    dphidy_cpc_ypencil = ZERO
    dphidy_cpc = ZERO
    fbcy_c4c = ZERO
    call transpose_x_to_y (phi_ccc, phi_ccc_ypencil, dm%dccc)
    call Get_y_1der_C2P_3D(phi_ccc_ypencil, dphidy_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_pr, fbcy_c4c )
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
    call Get_z_1der_C2P_3D(pphi_ccc_zpencil, dphidz_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_pr, fbcz_cc4 )
    call transpose_z_to_y (dphidz_ccp_zpencil, dphidz_ccp_ypencil, dm%dccp)
    call transpose_y_to_x (dphidz_ccp_ypencil, dphidz_ccp,         dm%dccp)
    uz = uz - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidz_ccp

    if(dm%is_thermo) then
      fl%gx = ux
      fl%gy = uy
      fl%gz = uz
    else 
      fl%qx = ux
      fl%qy = uy
      fl%qz = uz
    end if

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
      call Calculate_drhodt(fl, dm, isub)
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
    call wrt_3d_pt_debug (fl%pcor, dm%dccc, fl%iteration, isub, 'PhiRHS@RHS phi') ! debug_ww
    !call wrt_3d_all_debug(fl%pcor, dm%dccc,   fl%iteration, 'PhiRHS@RHS phi') ! debug_ww
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
    call wrt_3d_pt_debug (fl%pcor, dm%dccc,   fl%iteration, isub, 'phi@sol phi') ! debug_ww
    !call wrt_3d_all_debug(fl%pcor, dm%dccc,   fl%iteration, 'phi@sol phi') ! debug_ww
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
  subroutine Solve_momentum_eq(fl, dm, isub)
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
    integer,        intent(in)    :: isub

    call update_fbcy_cc_flow_halo(fl, dm)
!----------------------------------------------------------------------------------------------------------
! to calculate the rhs of the momenturn equation in stepping method
!----------------------------------------------------------------------------------------------------------
    call Compute_momentum_rhs(fl, dm, isub)
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
    call wrt_3d_pt_debug(fl%qx, dm%dpcc,   fl%iteration, isub, 'ux@bf divg') ! debug_ww
    call wrt_3d_pt_debug(fl%qy, dm%dcpc,   fl%iteration, isub, 'uy@bf divg') ! debug_ww
    call wrt_3d_pt_debug(fl%qz, dm%dccp,   fl%iteration, isub, 'uz@bf divg') ! debug_ww
    !write(*,*) 'qx', fl%qx(:, 1, 1), fl%qx(:, 8, 8)
    !write(*,*) 'qy', fl%qy(:, 1, 1), fl%qy(:, 8, 8)
    !write(*,*) 'qz', fl%qz(:, 1, 1), fl%qz(:, 8, 8)
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
    call solve_poisson(fl, dm, isub) ! test show above two methods gave the same results. 
!----------------------------------------------------------------------------------------------------------
! to update pressure
!----------------------------------------------------------------------------------------------------------
#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_mid_msg("Correcting the pressure term ...")
#endif
    fl%pres = fl%pres + fl%pcor
#ifdef DEBUG_STEPS
    call wrt_3d_pt_debug(fl%pres, dm%dccc,   fl%iteration, isub, 'pr@updated') ! debug_ww
#endif
!----------------------------------------------------------------------------------------------------------
! to update velocity/massflux correction
!----------------------------------------------------------------------------------------------------------
    !if(nrank == 0) call Print_debug_mid_msg("  Updating velocity/mass flux ...")
    call Correct_massflux(fl, fl%pcor, dm, isub)
#ifdef DEBUG_STEPS
  call wrt_3d_pt_debug(fl%qx, dm%dpcc,   fl%iteration, isub, 'ux@updated') ! debug_ww
  call wrt_3d_pt_debug(fl%qy, dm%dcpc,   fl%iteration, isub, 'uy@updated') ! debug_ww
  call wrt_3d_pt_debug(fl%qz, dm%dccp,   fl%iteration, isub, 'uz@updated') ! debug_ww
#endif
!----------------------------------------------------------------------------------------------------------
! to update velocity from gx gy gz 
!----------------------------------------------------------------------------------------------------------
  if(dm%is_thermo) then
    call Calculate_velocity_from_massflux(fl, dm)
  end if

    return
  end subroutine Solve_momentum_eq

end module eq_momentum_mod