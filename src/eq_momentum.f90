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
    use typeconvert_mod
    implicit none

    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(in ) :: dm
    integer,     intent(in ) :: isub

!----------------------------------------------------------------------------------------------------------
! result variables
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: mx_rhs_ypencil
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: mx_rhs_zpencil
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: mx_rhs_pfc

    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: my_rhs_ypencil ! 
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: my_rhs_zpencil ! 
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: my_rhs_pfc
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: my_rhs_pfc_ypencil

    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: mz_rhs_ypencil !
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: mz_rhs_zpencil ! 
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: mz_rhs_pfc
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: mz_rhs_pfc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: mz_rhs_pfc_zpencil
!----------------------------------------------------------------------------------------------------------
! temporary variables
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: apcc_xpencil
    apcc_ypencil
!----------------------------------------------------------------------------------------------------------
! intermediate variables : common
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: pres_ypencil ! p
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: pres_zpencil ! p

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: div_xpencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: div_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: div_zpencil
    
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

    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: qzix_pcp_xpencil ! qz^x, z-c1, common
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: qzix_pcp_zpencil ! qz^x, x-c3, no-thermal
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qziy_cpp_ypencil ! qz^y, z-c2, no-cly
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qziy_cpp_zpencil ! qz^y, y-c3, no-thermal, no-cly
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: qziz_ccc_zpencil ! qz^z, z-c3, common
    
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: qzdx_pcp_xpencil ! d(qz)/dx, z-v1, common
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: qzdx_pcp_zpencil ! d(qz)/dx, x-v3, common
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qzdy_cpp_ypencil ! d(qz)/dy, z-v2, common
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qzdy_cpp_zpencil ! d(qz)/dy, y-v3, no-cly
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qzdz_ccc_ypencil ! d(qz)/dz, y-v4, cly-only
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: qzdz_ccc_zpencil ! d(qz)/dz, z-v3, cly-only
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
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: gziy_cpp_zpencil  ! gz^y, y-c3, thermal-only
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: gziz_ccc_zpencil  ! gz^z, z-c3, thermal-only
    
    
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: mui00_Re_ccc_xpencil ! x-v1, thermal only
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: mui00_Re_ccc_ypencil ! y-v2, thermal only
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: mui00_Re_ccc_zpencil ! z-v3, thermal only
    real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: muixy_Re_ppc_xpencil ! y-v1, thermal only
    real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: muixy_Re_ppc_ypencil ! x-v2, thermal only
    real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: muixz_Re_pcp_xpencil ! z-v1, thermal only
    real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: muixz_Re_pcp_zpencil ! x-v3, thermal only
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: muiyz_Re_cpp_ypencil ! z-v2, z-v4, thermal only
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: muiyz_Re_cpp_zpencil ! y-v3, thermal only
!----------------------------------------------------------------------------------------------------------
!   intermediate variables : if fl%igravity == i
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: dDens_ypencil  ! d 
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: dDens_zpencil  ! d 
!----------------------------------------------------------------------------------------------------------
!   intermediate variables : if dm%icoordinate == ICYLINDRICAL
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qyriy_ccc_ypencil ! (qy/r)^y, y-c2, y-v4, z-v3, cly-only
    real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: qyriz_cpp_ypencil ! (qy/r)^z, z-c4, no-thermal
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: qyriz_cpp_zpencil ! (qy/r)^z, y-c3

    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: qyrdy_ccc_ypencil ! d(qy/r)/dy, y-v2
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
! bc variables
!----------------------------------------------------------------------------------------------------------
    real(WP), dimension( 4, dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: fbcx_pcc
    real(WP), dimension( 4, dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: fbcx_cpc
    real(WP), dimension( 4, dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: fbcx_ccp
    real(WP), dimension( 4, dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: fbcx_ccc

    real(WP), dimension( 4, dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: fbcy_pcc
    real(WP), dimension( 4, dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: fbcy_cpc
    real(WP), dimension( 4, dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: fbcy_ccp
    real(WP), dimension( 4, dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: fbcy_ccc

    real(WP), dimension( 4, dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: fbcz_pcc
    real(WP), dimension( 4, dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: fbcz_cpc
    real(WP), dimension( 4, dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: fbcz_ccp
    real(WP), dimension( 4, dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: fbcz_ccc
!----------------------------------------------------------------------------------------------------------
! others
!----------------------------------------------------------------------------------------------------------
    real(WP) :: one_third_rre, two_third_rre, two_rre
    integer  :: i
    real(WP) :: rhsx_bulk, rhsz_bulk

!==========================================================================================================
! preparation of constants to be used.
!==========================================================================================================
#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_start_msg("Compute_momentum_rhs at isub = "//trim(int2str(isub)))
#endif

    one_third_rre = ONE_THIRD * fl%rre
    two_third_rre = TWO_THIRD * fl%rre
          two_rre = TWO * fl%rre
!==========================================================================================================
! preparation of intermediate variables to be used - common
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
!    p --> p_ypencil --> p_zpencil 
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y (fl%pres,      pres_ypencil, dm%dccc) 
    call transpose_y_to_z (pres_ypencil, pres_zpencil, dm%dccc)
!----------------------------------------------------------------------------------------------------------
!   d --> d_ypencil --> d_zpencil
!----------------------------------------------------------------------------------------------------------    
    if(fl%igravity /= 0) then
      call transpose_x_to_y (tm%dDens,      dDens_ypencil, dm%dccc)
      call transpose_y_to_z (dDens_ypencil, dDens_zpencil, dm%dccc)
    end if
!----------------------------------------------------------------------------------------------------------
!    qx 
!    | -[ipx]-> qxix_ccc_xpencil (common)
!    | -[1dx]-> qxdx_ccc_xpencil(common)
!    | -[x2y]-> qx_ypencil(temp) 
!               | -[ipy]-> qxiy_ppc_ypencil(common) -[y2x]-> qxiy_ppc_xpencil(no thermal)
!               | -[1dy]-> qxdy_ppc_ypencil(common) -[y2x]-> qxdy_ppc_xpencil(common)
!               | -[y2z]-> qx_zpencil(temp) 
!                          | -[ipz]-> qxiz_pcp_zpencil(common) -[z2y]-> qxiz_pcp_ypencil(temp) -[y2x]-> qxiz_pcp_xpencil(no-thermal)
!                          | -[1dz]-> qxdz_pcp_zpencil(common) -[z2y]-> qxdz_pcp_ypencil(temp) -[y2x]-> qxdz_pcp_xpencil(common)
!----------------------------------------------------------------------------------------------------------
    call Get_x_midp_P2C_3D(fl%qx, qxix_ccc_xpencil, dm, dm%ibcx(:, 1), dm%fbcx_qx(:, :, :))
    call Get_x_1st_derivative_P2C_3D(fl%qx, qxdx_ccc_xpencil, dm, dm%ibcx(:, 1), dm%fbcx_qx(:, :, :))

    call transpose_x_to_y (fl%qx, apcc_ypencil, dm%dpcc) !qx_ypencil
    call Get_y_1st_derivative_C2P_3D(apcc_ypencil, qxdy_ppc_ypencil, dm, dm%ibcy(:, 1), dm%fbcy_qx(:, :, :))
    call transpose_y_to_z(qxdy_ppc_ypencil, qxdy_ppc_xpencil, dm%dppc)
    call Get_y_midp_C2P_3D(apcc_ypencil, qxiy_ppc_ypencil, dm, dm%ibcy(:, 1), dm%fbcy_qx(:, :, :))
    if(.not. dm%is_thermo) then
      call transpose_y_to_x (qxiy_ppc_ypencil, qxiy_ppc_xpencil, dm%dppc)
    end if

    call transpose_y_to_z (apcc_ypencil, apcc_zpencil, dm%dpcc)!qx_zpencil
    call Get_z_midp_C2P_3D(apcc_zpencil, qxiz_pcp_zpencil, dm, dm%ibcz(:, 1), dm%fbcz_qx(:, :, :))
    call transpose_z_to_y (qxiz_pcp_zpencil, apcp_ypencil, dm%dpcc)!qxiz_pcp_ypencil
    if(.not. dm%is_thermo) then
      call transpose_y_to_x (apcp_ypencil, qxiz_pcp_xpencil, dm%dpcp)
    end if

    call Get_z_1st_derivative_C2P_3D(apcc_zpencil, qxdz_pcp_zpencil, dm, dm%ibcz(:, 1), dm%fbcz_qx(:, :, :))
    call transpose_z_to_y(qxdz_pcp_zpencil, apcp_ypencil,     dm%dpc) !qxdz_pcp_ypencil
    call transpose_y_to_x(apcp_ypencil,     qxdz_pcp_xpencil, dm%dpc)
!----------------------------------------------------------------------------------------------------------
!    qy
!    | -[ipx]-> qyix_ppc_xpencil(common) -[x2y]-> qyix_ppc_ypencil(no thermal) 
!    | -[1dx]-> qydx_ppc_xpencil(common) -[x2y]-> qydx_ppc_ypencil(common)
!    | -[x2y]-> qy_ypencil(temp) 
!               | -[ipy]-> qyiy_ccc_ypencil(no thermal or no cyl)
!               | -[1dy]-> qydy_ccc_ypencil(no-cly)
!               | -[y2z]-> qy_zpencil(temp) 
!                          | -[ipz]-> qyiz_cpp_zpencil(no-cly) -[z2y]-> qyiz_cpp_ypencil(no-thermal)           
!                          | -[1dz]-> qydz_cpp_zpencil(no-cly) -[z2y]-> qydz_cpp_ypencil(no-cly)
!----------------------------------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3D(fl%qy, qydx_ppc_xpencil, dm, dm%ibcx(:, 2), dm%fbcx_qy(:, :, :))
    call transpose_x_to_y(qydx_ppc_xpencil, qydx_ppc_ypencil, dm%dppc)
    call Get_x_midp_C2P_3D(fl%qy, qyix_ppc_xpencil, dm, dm%ibcx(:, 2), dm%fbcx_qy(:, :, :))
    if(.not. dm%is_thermo) then
      call transpose_x_to_y (qyix_ppc_xpencil, qyix_ppc_ypencil, dm%dppc)
    end if

    call transpose_x_to_y (fl%qy, acpc_ypencil, dm%dcpc) !qy_ypencil
    if(dm%icoordinate == ICARTESIAN) then
      call Get_y_1st_derivative_P2C_3D(acpc_ypencil, qydy_ccc_ypencil, dm, dm%ibcy(:, 2), dm%fbcy_qy(:, :, :))
    end if
    if((.not. dm%is_thermo) .or. (dm%icoordinate == ICARTESIAN)) then
      call Get_y_midp_P2C_3D(acpc_ypencil, qyiy_ccc_ypencil, dm, dm%ibcy(:, 2), dm%fbcy_qy(:, :, :) )
    end if

    call transpose_z_to_y(acpc_ypencil, acpc_zpencil, dm%dcpc)!qy_zpencil
    if(dm%icoordinate == ICARTESIAN)then
      call Get_z_1st_derivative_C2P_3D(acpc_zpencil, qydz_cpp_zpencil, dm, dm%ibcz(:, 2), dm%fbcz_qy(:, :, :))
      call transpose_z_to_y(qydz_cpp_zpencil, qydz_cpp_ypencil, dm%dcpp)
    end if
    if((.not. dm%is_thermo) .or. (dm%icoordinate == ICARTESIAN)) then
      call Get_z_midp_C2P_3D(acpc_zpencil, qyiz_cpp_zpencil, dm, dm%ibcz(:, 2), dm%fbcz_qy(:, :, :) )
      if(.not. dm%is_thermo) then
        call transpose_z_to_y(qyiz_cpp_zpencil, qyiz_cpp_ypencil, dm%dcpp)
      end if
    end if
!----------------------------------------------------------------------------------------------------------
!    qz 
!    | -[1dx]-> qzdx_pcp_xpencil(common) -[x2y]-> qzdx_pcp_ypencil(temp) -[y2z]-> qzdx_pcp_zpencil(common)
!    | -[ipx]-> qzix_pcp_xpencil(common) -[x2y]-> qzix_pcp_ypencil(temp) -[y2z]-> qzix_pcp_zpencil(no-thermal)
!    | -[x2y]-> qz_ypencil(temp) 
!               | -[1dy]-> qzdy_cpp_ypencil(common) -[y2z]-> qzdy_cpp_zpencil(no-cly)
!               | -[ipy]-> qziy_cpp_ypencil(no-cly) -[y2z]-> qziy_cpp_zpencil(no cly, no thermal)
!               | -[y2z]-> qz_zpencil(temp) 
!                          | -[ipz]-> qziz_ccc_zpencil(common)
!                          | -[1dz]-> qzdz_ccc_zpencil(cly-only) -[y2z]-> qzdz_ccc_zpencil(cly-only)
!----------------------------------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3D(fl%qz, qzdx_pcp_xpencil, dm, dm%ibcx(:, 3), dm%fbcx_qz(:, :, :))
    call transpose_x_to_y(qzdx_pcp_xpencil, apcp_ypencil, dm%dpcp)
    call transpose_y_to_z(apcp_ypencil, qzdx_pcp_zpencil, dm%dpcp)

    call Get_x_midp_C2P_3D(fl%qz, qzix_pcp_xpencil, dm, dm%ibcx(:, 3), dm%fbcx_qz(:, :, :))
    if(.not. dm%is_thermo) then
      call transpose_x_to_y(qzix_pcp_xpencil, apcp_ypencil, dm%dpcp) !qzix_pcp_ypencil
      call transpose_y_to_z(apcp_ypencil, qzix_pcp_zpencil, dm%dpcp)
    end if

    call transpose_x_to_y(fl%qz,        accp_ypencil, dm%dccp) ! qz_ypencil
    call Get_y_1st_derivative_C2P_3D(accp_ypencil, qzdy_cpp_ypencil, dm, dm%ibcy(:, 3), dm%fbcy_qz(:, :, :))
    if(dm%icoordinate == ICARTESIAN)then
      call transpose_y_to_z(qzdy_cpp_ypencil, qzdy_cpp_zpencil, dm%dcpp)
    end if

    if((.not. dm%is_thermo) .or. (dm%icoordinate == ICARTESIAN)) then
      call Get_y_midp_C2P_3D(accp_ypencil, qziy_cpp_ypencil, dm, dm%ibcy(:, 3), dm%fbcy_qz(:, :, :) )
      if(.not. dm%is_thermo) then
        call transpose_y_to_z(qziy_cpp_ypencil, qziy_cpp_zpencil, dm%dcpp)
      end if
    end if

    call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp) ! qz_zpencil
    call Get_z_midp_P2C_3D(accp_zpencil, qziz_ccc_zpencil, dm, dm%ibcz(:, 3), dm%fbcz_qz(:, :, :))
    if(dm%icoordinate == ICYLINDRICAL) then
      call Get_z_1st_derivative_P2C_3D(accp_zpencil, qzdz_ccc_zpencil, dm, dm%ibcz(:, 3), dm%fbcz_qz(:, :, :))
      call transpose_z_to_y(qzdz_ccc_zpencil, qzdz_ccc_ypencil, dm%dccc)
    end if
!==========================================================================================================
! preparation of intermediate variables to be used - thermal only
!==========================================================================================================
    if((dm%is_thermo) then
!----------------------------------------------------------------------------------------------------------
!    gx 
!    | -[ipx]-> gxix_ccc_xpencil
!    | -[x2y]-> gx_ypencil(temp) 
!               | -[ipy]-> gxiy_ppc_ypencil(temp) -[y2x]-> gxiy_ppc_xpencil
!               | -[y2z]-> gx_zpencil(temp) 
!                          | -[ipz]-> gxiz_pcp_zpencil(temp) -[z2y]-> gxiz_pcp_ypencil(temp) -[y2x]-> gxiz_pcp_xpencil
!----------------------------------------------------------------------------------------------------------
      call Get_x_midp_P2C_3D(fl%gx, gxix_ccc_xpencil, dm, dm%ibcx(:, 1), dm%fbcx_gx(:, :, :))
      
      call transpose_x_to_y (fl%gx, apcc_ypencil, dm%dpcc) !gx_ypencil
      call Get_y_midp_C2P_3D(apcc_ypencil, appc_ypencil, dm, dm%ibcy(:, 1), dm%fbcy_gx(:, :, :))
      call transpose_y_to_x (appc_ypencil, gxiy_ppc_xpencil, dm%dppc)

      call transpose_y_to_z (apcc_ypencil, apcc_zpencil, dm%dpcc)!gx_zpencil
      call Get_z_midp_C2P_3D(apcc_zpencil, apcp_zpencil, dm, dm%ibcz(:, 1), dm%fbcz_gx(:, :, :))
      call transpose_z_to_y (apcp_zpencil, apcp_ypencil, dm%dpcc)!qxiz_pcp_ypencil
      call transpose_y_to_x (apcp_ypencil, gxiz_pcp_xpencil, dm%dpcp)
!----------------------------------------------------------------------------------------------------------
!    gy
!    | -[ipx]-> gyix_ppc_xpencil(temp) -[x2y]-> gyix_ppc_ypencil
!    | -[x2y]-> gy_ypencil(temp) 
!               | -[ipy]-> gyiy_ccc_ypencil
!               | -[y2z]-> gy_zpencil(temp) 
!                          | -[ipz]-> gyiz_cpp_zpencil(temp) -[z2y]-> gyiz_cpp_ypencil          
!----------------------------------------------------------------------------------------------------------
      call Get_x_midp_C2P_3D(fl%gy, appc_xpencil, dm, dm%ibcx(:, 2), dm%fbcx_gy(:, :, :))
      call transpose_x_to_y (appc_xpencil, gyix_ppc_ypencil, dm%dppc)

      call transpose_x_to_y (fl%gy, acpc_ypencil, dm%dcpc) !gy_ypencil
      call Get_y_midp_P2C_3D(acpc_ypencil, gyiy_ccc_ypencil, dm, dm%ibcy(:, 2), dm%fbcy_gy(:, :, :) )

      call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc)!qy_zpencil
      call Get_z_midp_C2P_3D(acpc_zpencil, acpp_zpencil, dm, dm%ibcz(:, 2), dm%fbcz_gy(:, :, :) )
      call transpose_z_to_y(acpp_zpencil, gyiz_cpp_ypencil, dm%dcpp)
!----------------------------------------------------------------------------------------------------------
!    gz 
!    | -[ipx]-> gzix_pcp_xpencil(temp) -[x2y]-> gzix_pcp_ypencil(temp) -[y2z]-> gzix_pcp_zpencil
!    | -[x2y]-> gz_ypencil(temp) 
!               | -[ipy]-> gziy_cpp_ypencil(temp) -[y2z]-> gziy_cpp_zpencil
!               | -[y2z]-> gz_zpencil(temp) 
!                          | -[ipz]-> gziz_ccc_zpencil
!----------------------------------------------------------------------------------------------------------
      call Get_x_midp_C2P_3D(fl%gz, apcp_xpencil, dm, dm%ibcx(:, 3), dm%fbcx_gz(:, :, :))
      call transpose_x_to_y(apcp_xpencil, apcp_ypencil, dm%dpcp)
      call transpose_y_to_z(apcp_ypencil, gzix_pcp_zpencil, dm%dpcp)

      call transpose_x_to_y(fl%gz, accp_ypencil, dm%dccp) ! gz_ypencil
      call Get_y_midp_C2P_3D(accp_ypencil, acpp_ypencil, dm, dm%ibcy(:, 3), dm%fbcy_gz(:, :, :) )
      call transpose_y_to_z(acpp_ypencil, gziy_cpp_zpencil, dm%dcpp)

      call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp) ! gz_zpencil
      call Get_z_midp_P2C_3D(accp_zpencil, gziz_ccc_zpencil, dm, dm%ibcz(:, 3), dm%fbcz_gz(:, :, :))
!----------------------------------------------------------------------------------------------------------
!    m = mui00_ccc_xpencil
!    | -[ipx]-> muix_pcc_xpencil(temp) -[x2y]-> muix_pcc_ypencil(temp)
!                                               | -[ipy]-> muixy_ppc_ypencil -[y2z]-> muixy_ppc_xpencil
!                                               | -[y2z]-> muix_pcc_zpencil(temp) -[ipz]-> muixz_pcp_zpencil -[z2y]-> muixz_pcp_ypencil (temp) -[y2z]-> muixz_pcp_xpencil
!    | -[x2y]-> mui00_ccc_ypencil -[y2z]-> mui00_ccc_zpencil
!              | -[ipy]-> muiy_cpc_ypencil(temp) -[y2z]-> muiy_cpc_zpencil(temp) -[ipz]-> muiyz_cpp_zpencil -[z2y]-> muiyz_cpp_ypencil
!----------------------------------------------------------------------------------------------------------
      mui00_ccc_xpencil = tm%mVisc
      call transpose_x_to_y(mui00_ccc_xpencil, mui00_ccc_ypencil, dm%dccc)
      call transpose_y_to_z(mui00_ccc_ypencil, mui00_ccc_zpencil, dm%dccc)

      call Get_x_midp_C2P_3D(mui00_ccc_xpencil, apcc_xpencil, dm, dm%ibcx(:, 5), tm%fbcx_ftp(:, :, :)%m) ! muix_pcc_xpencil
      call transpose_x_to_y(apcc_xpencil, apcc_ypencil, dm%dpcc) ! muix_pcc_ypencil
      
      call Get_y_midp_C2P_3D(apcc_ypencil, muixy_ppc_ypencil, dm, dm%ibcy(:, 5), tm%fbcy_ftp(:, :, :)%m)
      call transpose_y_to_x(muixy_ppc_ypencil, muixy_ppc_xpencil, dm%dppc) 

      call transpose_y_to_z(apcc_ypencil, apcc_zpencil, dm%dpcc) ! muix_pcc_zpencil
      call Get_z_midp_C2P_3D(apcc_zpencil, muixz_pcp_zpencil, dm, dm%ibcz(:, 5), tm%fbcz_ftp(:, :, :)%m)
      call transpose_z_to_y(muixz_pcp_zpencil, apcp_zpencil, dm%dpcp)
      call transpose_y_to_x(apcp_zpencil, muixz_pcp_xpencil, dm%dpcp)

      call Get_y_midp_C2P_3D(mui00_ccc_ypencil, acpc_ypencil, dm, dm%ibcy(:, 5), tm%fbcy_ftp(:, :, :)%m) ! muiy_cpc_ypencil
      call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc) !muiy_cpc_zpencil
      call Get_z_midp_C2P_3D(acpc_zpencil, muiyz_cpp_zpencil, dm, dm%ibcz(:, 5), tm%fbcz_ftp(:, :, :)%m)
      call transpose_z_to_y(muiyz_cpp_zpencil, muiyz_cpp_ypencil, dm%dcpp)
    end if

!==========================================================================================================
! preparation of intermediate variables to be used - cylindrical only
!==========================================================================================================
    if((dm%icoordinate == ICYLINDRICAL) then
!----------------------------------------------------------------------------------------------------------
!    qy/r=qyr
!    | -[x2y]-> qyr_ypencil(temp) 
!               | -[ipy]-> qyriy_ccc_ypencil
!               | -[1dy]-> qyrdy_ccc_ypencil
!               | -[y2z]-> qyr_zpencil(temp) 
!                          | -[ipz]-> qyriz_cpp_zpencil -[z2y]-> qyriz_cpp_ypencil        
!                          | -[1dz]-> qyrdz_cpp_zpencil -[z2y]-> qyrdz_cpp_ypencil
!----------------------------------------------------------------------------------------------------------
      acpc_xpencil = fl%qy * dm%
      call multiple_cylindrical_rn(acpc_xpencil, dm%dcpc, dm%rpi, 1, IPENCIL(1)) ! qr/r
      call transpose_x_to_y(acpc_xpencil, acpc_ypencil, dm%dcpc)
      
      call Get_y_midp_P2C_3D(acpc_ypencil, qyriy_ccc_ypencil, dm, dm%ibcy(:, 2), dm%fbcz_qyr(:, :, :))
       

    end if

!==========================================================================================================
! the RHS of x-momentum equation
!==========================================================================================================
    i = 1
    fl%mx_rhs      = ZERO
    mx_rhs_ypencil = ZERO
    mx_rhs_zpencil = ZERO
    mx_rhs_pfc     = ZERO
!----------------------------------------------------------------------------------------------------------
! X-mom convection term (x-c1/3), X-pencil: -d(gx^x * qx^x)/dx at (i', j, k)
!----------------------------------------------------------------------------------------------------------  
    if ( .not. dm%is_thermo) then
      call Get_x_1st_derivative_C2P_3D(-qxix_ccc_xpencil * qxix_ccc_xpencil, apcc_xpencil, dm, dm%ibcx(:, i), dm%fbcx_var(:, :, :, i) * dm%fbcx_var(:, :, :, 1) )
    else
      call Get_x_1st_derivative_C2P_3D(-gxix_ccc_xpencil * qxix_ccc_xpencil, apcc_xpencil, dm, dm%ibcx(:, i), dm%fbcx_var(:, :, :, i) * dm%fbcx_var(:, :, :, 1 + NBC) )
    end if
    fl%mx_rhs = fl%mx_rhs + apcc
!----------------------------------------------------------------------------------------------------------
! X-mom convection term (x-c2/3), Y-pencil: -d(<gy>^x * <qx>^y)/dy * (1/r) at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_y_1st_derivative_P2C_3D(-qy_ppc_ypencil * qx_ppc_ypencil, apcc_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i) * dm%fbcy_var(:, :, :, 2) )
    else
      call Get_y_1st_derivative_P2C_3D(-gy_ppc_ypencil * qx_ppc_ypencil, apcc_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i) * dm%fbcy_var(:, :, :, 2 + NBC) )
    end if
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_ypencil, dm%dpcc, dm%rci, 1, IPENCIL(2))
    mx_rhs_ypencil = mx_rhs_ypencil + apcc_ypencil
!----------------------------------------------------------------------------------------------------------
! X-mom convection term (x-c3/3), Z-pencil: -d(<gz>^x * <qx>^z)/dz at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_z_1st_derivative_P2C_3D(-qz_pcp_zpencil * qx_pcp_zpencil, apcc_zpencil, dm, dm%ibcz(:, i), dm%fbcz_var(:, :, :, i) * dm%fbcz_var(:, :, :, 3) )
    else
      call Get_z_1st_derivative_P2C_3D(-gz_pcp_zpencil * qx_pcp_zpencil, apcc_zpencil, dm, dm%ibcz(:, i), dm%fbcz_var(:, :, :, i) * dm%fbcz_var(:, :, :, 3 + NBC) )
    end if
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_zpencil, dm%dpcc, dm%rci, 2, IPENCIL(3))
    mx_rhs_zpencil = mx_rhs_zpencil + apcc_zpencil
!----------------------------------------------------------------------------------------------------------
! X-mom pressure gradient in x direction, X-pencil, d(sigma_1 p)
!----------------------------------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3D( -fl%pres, apcc, dm, dm%ibcx(:, 4), dm%fbcx_var(:, :, :, 4) )
    mx_rhs_pfc=  mx_rhs_pfc + apcc
!----------------------------------------------------------------------------------------------------------
! X-mom diffusion term (x-v-1/3), X-pencil, d_x [2 * mu * (d_x qx - 1/3*div) ] at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    if ( dm%is_thermo) then
      accc = (dqxdx_ccc_xpencil - ONE_THIRD * div) * TWO * mu_Re_ccc_xpencil
    else
      accc = (dqxdx_ccc_xpencil - ONE_THIRD * div) * TWO * fl%rre ! accc = 2 * mu * (d_x qx - 1/3*div)  at (i, j, k)
    end if
    call Get_x_1st_derivative_C2P_3D(accc,  apcc, dm, dm%ibcx(:, 1), dm%fbcx_var(:, :, :, 1)) ! check 
    fl%mx_rhs = fl%mx_rhs + apcc
!----------------------------------------------------------------------------------------------------------
! X-mom diffusion term (x-v-2/3), Y-pencil, (1/r) * d_y [mu^{x,y} * (d_x q_y + r * d_y q_x) ] at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    appc_ypencil2 = dqxdy_ppc_ypencil
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(appc_ypencil2, dm%dppc, ONE/dm%rci, 1, IPENCIL(2))
    if ( dm%is_thermo) then
      appc_ypencil = ( appc_ypencil2 + dqydx_ppc_ypencil) * mu_Re_ppc_ypencil
    else
      appc_ypencil = ( appc_ypencil2 + dqydx_ppc_ypencil) * fl%rre
    end if
    call Get_y_1st_derivative_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%ibcy(:, 1), dm%fbcy_var(:, :, :, 1))
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_ypencil, dm%dpcc, dm%rci, 1, IPENCIL(2))
    mx_rhs_ypencil =  mx_rhs_ypencil + apcc_ypencil
!----------------------------------------------------------------------------------------------------------
! X-mom diffusion term (x-v-3/3), Z-pencil, (1/r^2) * d_z [mu^{x,z} * (d_x q_z + d_z q_x) ] at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    if ( dm%is_thermo) then
      apcp_zpencil = ( dqxdz_pcp_zpencil + dqzdx_pcp_zpencil) * mu_Re_pcp_zpencil
    else
      apcp_zpencil = ( dqxdz_pcp_zpencil + dqzdx_pcp_zpencil) * fl%rre
    end if
    call Get_z_1st_derivative_P2C_3D(apcp_zpencil, apcc_zpencil, dm, dm%ibcz(:, 1), dm%fbcz_var(:, :, :, 1))
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_zpencil, dm%dpcc, dm%rci, 2, IPENCIL(3))
    mx_rhs_zpencil =  mx_rhs_zpencil + apcc_zpencil
!----------------------------------------------------------------------------------------------------------
! x-mom: convert all terms to rhs
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_x (mx_rhs_ypencil, apcc, dm%dpcc)
    fl%mx_rhs =  fl%mx_rhs + apcc
    call transpose_z_to_y (mx_rhs_zpencil, apcc_ypencil, dm%dpcc)
    call transpose_y_to_x (apcc_ypencil, apcc, dm%dpcc)
    fl%mx_rhs =  fl%mx_rhs + apcc 
    




































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
    i = 1

!==========================================================================================================
! the RHS of y-momentum equation
!==========================================================================================================
    fl%my_rhs = ZERO
    my_rhs_ypencil = ZERO
    my_rhs_zpencil = ZERO
    
    my_rhs_pfc  = ZERO
    my_rhs_pfc_ypencil = ZERO

    acpc = ZERO
    acpc_ypencil = ZERO
    acpc_zpencil = ZERO
    i = 2
!----------------------------------------------------------------------------------------------------------
! Y-mom convection term (y-c1/3), X-pencil: -d(gx^y * qy^x)/dx at (i, j', k)
!----------------------------------------------------------------------------------------------------------  
    if ( .not. dm%is_thermo) then
      call Get_x_1st_derivative_P2C_3D( -qx_ppc * qy_ppc, acpc, dm, dm%ibcx(:, i), dm%fbcx_var(:, :, :, i) * dm%fbcx_var(:, :, :, 1) )
    else
      call Get_x_1st_derivative_P2C_3D( -gx_ppc * qy_ppc, acpc, dm, dm%ibcx(:, i), dm%fbcx_var(:, :, :, i) * dm%fbcx_var(:, :, :, 1 + NBC) )
    end if
    fl%my_rhs = fl%my_rhs + acpc
!----------------------------------------------------------------------------------------------------------
! Y-mom convection term (x-c2/3), Y-pencil: -d(<gy>^y * <qy/r>^x)/dy * (1/r) at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_y_1st_derivative_C2P_3D(-qy_ccc_ypencil * qy_ccc_ypencil, acpc_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i) * dm%fbcy_var(:, :, :, 2) )
    else
      call Get_y_1st_derivative_C2P_3D(-gy_ccc_ypencil * qy_ccc_ypencil, acpc_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i) * dm%fbcy_var(:, :, :, 2 + NBC) )
    end if
    my_rhs_ypencil = my_rhs_ypencil + acpc_ypencil


    
    if ( .not. dm%is_thermo) then
      call Get_y_1st_derivative_P2C_3D(-qy_ppc_ypencil * qx_ppc_ypencil, apcc_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i) * dm%fbcy_var(:, :, :, 2) )
    else
      call Get_y_1st_derivative_P2C_3D(-gy_ppc_ypencil * qx_ppc_ypencil, apcc_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i) * dm%fbcy_var(:, :, :, 2 + NBC) )
    end if
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_ypencil, dm%dpcc, dm%rci, 1, IPENCIL(2))
    mx_rhs_ypencil = mx_rhs_ypencil + apcc_ypencil
!----------------------------------------------------------------------------------------------------------
! Y-mom convection term (x-c3/3), Z-pencil: -d(<gz>^x * <qx>^z)/dz at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_z_1st_derivative_P2C_3D(-qz_pcp_zpencil * qx_pcp_zpencil, apcc_zpencil, dm, dm%ibcz(:, i), dm%fbcz_var(:, :, :, i) * dm%fbcz_var(:, :, :, 3) )
    else
      call Get_z_1st_derivative_P2C_3D(-gz_pcp_zpencil * qx_pcp_zpencil, apcc_zpencil, dm, dm%ibcz(:, i), dm%fbcz_var(:, :, :, i) * dm%fbcz_var(:, :, :, 3 + NBC) )
    end if
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_zpencil, dm%dpcc, dm%rci, 2, IPENCIL(3))
    mx_rhs_zpencil = mx_rhs_zpencil + apcc_zpencil
!----------------------------------------------------------------------------------------------------------
! Y-mom pressure gradient in x direction, X-pencil, d(sigma_1 p)
!----------------------------------------------------------------------------------------------------------
    call Get_x_1st_derivative_C2P_3D( -fl%pres, apcc, dm, dm%ibcx(:, 4), dm%fbcx_var(:, :, :, 4) )
    mx_rhs_pfc=  mx_rhs_pfc + apcc
!----------------------------------------------------------------------------------------------------------
! Y-mom diffusion term (x-v-1/3), X-pencil, d_x [2 * mu * (d_x qx - 1/3*div) ] at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    if ( dm%is_thermo) then
      accc = (dqxdx_ccc_xpencil - ONE_THIRD * div) * TWO * mu_Re_ccc_xpencil
    else
      accc = (dqxdx_ccc_xpencil - ONE_THIRD * div) * TWO * fl%rre ! accc = 2 * mu * (d_x qx - 1/3*div)  at (i, j, k)
    end if
    call Get_x_1st_derivative_C2P_3D(accc,  apcc, dm, dm%ibcx(:, 1), dm%fbcx_var(:, :, :, 1)) ! check 
    fl%mx_rhs = fl%mx_rhs + apcc
!----------------------------------------------------------------------------------------------------------
! Y-mom diffusion term (x-v-2/3), Y-pencil, (1/r) * d_y [mu^{x,y} * (d_x q_y + r * d_y q_x) ] at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    appc_ypencil2 = dqxdy_ppc_ypencil
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(appc_ypencil2, dm%dppc, ONE/dm%rci, 1, IPENCIL(2))
    if ( dm%is_thermo) then
      appc_ypencil = ( appc_ypencil2 + dqydx_ppc_ypencil) * mu_Re_ppc_ypencil
    else
      appc_ypencil = ( appc_ypencil2 + dqydx_ppc_ypencil) * fl%rre
    end if
    call Get_y_1st_derivative_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%ibcy(:, 1), dm%fbcy_var(:, :, :, 1))
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_ypencil, dm%dpcc, dm%rci, 1, IPENCIL(2))
    mx_rhs_ypencil =  mx_rhs_ypencil + apcc_ypencil
!----------------------------------------------------------------------------------------------------------
! Y-mom diffusion term (x-v-3/3), Z-pencil, (1/r^2) * d_z [mu^{x,z} * (d_x q_z + d_z q_x) ] at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    if ( dm%is_thermo) then
      apcp_zpencil = ( dqxdz_pcp_zpencil + dqzdx_pcp_zpencil) * mu_Re_pcp_zpencil
    else
      apcp_zpencil = ( dqxdz_pcp_zpencil + dqzdx_pcp_zpencil) * fl%rre
    end if
    call Get_z_1st_derivative_P2C_3D(apcp_zpencil, apcc_zpencil, dm, dm%ibcz(:, 1), dm%fbcz_var(:, :, :, 1))
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(apcc_zpencil, dm%dpcc, dm%rci, 2, IPENCIL(3))
    mx_rhs_zpencil =  mx_rhs_zpencil + apcc_zpencil
!----------------------------------------------------------------------------------------------------------
! Y-mom: convert all terms to rhs
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_x (mx_rhs_ypencil, apcc, dm%dpcc)
    fl%mx_rhs =  fl%mx_rhs + apcc
    call transpose_z_to_y (mx_rhs_zpencil, apcc_ypencil, dm%dpcc)
    call transpose_y_to_x (apcc_ypencil, apcc, dm%dpcc)
    fl%mx_rhs =  fl%mx_rhs + apcc 













!----------------------------------------------------------------------------------------------------------
! X-pencil : X-mom diffusion term (x-v1-1/7), \mu^x * LL1(ux) at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    !call Get_x_2nd_derivative_P2P_3D(fl%qx, apcc, dm, dm%ibcx(:, 1)) ! check
    accc = ZERO
    apcc = ZERO
    call Get_x_1st_derivative_P2C_3D(fl%qx, accc, dm, dm%ibcx(:, i), dm%fbcx_var(:, :, :, i))
    call Get_x_1st_derivative_C2P_3D(accc,  apcc, dm, dm%ibcx(:, i), dm%fbcx_var(:, :, :, i))

    if ( .not. dm%is_thermo) then
      fl%mx_rhs = fl%mx_rhs +         fl%rre * apcc
    else
      fl%mx_rhs = fl%mx_rhs + m_pcc * fl%rre * apcc
    end if
!----------------------------------------------------------------------------------------------------------
! Y-pencil : X-mom diffusion term (x-v1-2/7), \mu^x * LL2(ux) at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    !call Get_y_2nd_derivative_C2C_3D(qx_ypencil, apcc_ypencil, dm, dm%ibcy(:, 1), dm%fbcy_var(:, :, :, 1))
    appc_ypencil = ZERO
    apcc_ypencil = ZERO
    call Get_y_1st_derivative_C2P_3D(qx_ypencil,   appc_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i))
    call Get_y_1st_derivative_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i))

    if ( .not. dm%is_thermo ) then
      mx_rhs_ypencil = mx_rhs_ypencil +                 fl%rre * apcc_ypencil
    else
      mx_rhs_ypencil = mx_rhs_ypencil + m_pcc_ypencil * fl%rre * apcc_ypencil
    end if
!----------------------------------------------------------------------------------------------------------
! Z-pencil : X-mom diffusion term (x-v1-3/7), \mu^x * LL3(ux) at (i', j, k)
!----------------------------------------------------------------------------------------------------------
    !call Get_z_2nd_derivative_C2C_3D(qx_zpencil, apcc_zpencil, dm, dm%ibcz(:, 1), dm%fbcz_var(:, 1))
    apcp_zpencil = ZERO 
    apcc_zpencil = ZERO
    call Get_z_1st_derivative_C2P_3D(qx_zpencil,   apcp_zpencil, dm, dm%ibcz(:, i), dm%fbcz_var(:, :, :, i))
    call Get_z_1st_derivative_P2C_3D(apcp_zpencil, apcc_zpencil, dm, dm%ibcz(:, i), dm%fbcz_var(:, :, :, i))

    if ( .not. dm%is_thermo) then
      mx_rhs_zpencil = mx_rhs_zpencil +                 fl%rre * apcc_zpencil
    else
      mx_rhs_zpencil = mx_rhs_zpencil + m_pcc_zpencil * fl%rre * apcc_zpencil
    end if

    if ( dm%is_thermo) then
!----------------------------------------------------------------------------------------------------------
! x-pencil : X-mom, gravity force in x direction
!----------------------------------------------------------------------------------------------------------
      if(fl%igravity == i .or. fl%igravity == -i)  then
        call Get_x_midp_C2P_3D(tm%dDens, apcc, dm, dm%ibcx(:, 5), dm%ftpbcx_var(:, :, :)%d )
        mx_rhs_pfc =  mx_rhs_pfc + fl%fgravity(i) * apcc
      end if
!----------------------------------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v2/7), \mu^x * 1/3 * d (div)/dx at (i', j, k)
!----------------------------------------------------------------------------------------------------------
      fbcx = ZERO
      call Get_x_1st_derivative_C2P_3D(div, apcc, dm, dm%ibcx(:, i), fbcx ) ! apcc = d(div)/dx_xpencil
      fl%mx_rhs = fl%mx_rhs + one_third_rre * m_pcc * apcc
!----------------------------------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v3/7), -2/3 * d\mu/dx * (div(u)^x) +  
!                                                2 * d\mu/dx * du/dx
!----------------------------------------------------------------------------------------------------------
      fbcx = ZERO ! check
      call Get_x_midp_C2P_3D (div, apcc, dm, dm%ibcx(:, i), fbcx ) ! apcc = div_pcc
      fl%mx_rhs = fl%mx_rhs - two_third_rre * dmdx_pcc * apcc
      call Get_x_1st_derivative_P2P_3D(fl%qx, apcc, dm, dm%ibcx(:, 1), dm%fbcx_var(:, :, :, 1) ) ! apcc = d(qx)/dx_pcc
      fl%mx_rhs = fl%mx_rhs + two_rre       * dmdx_pcc * apcc
!----------------------------------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v4/7), d(mu^x)/dy * d(qy^y)/dx at (i', j, k)
!----------------------------------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3D(qy_ccc, apcc, dm, dm%ibcx(:, 2), dm%fbcx_var(:, :, :, 2) ) !apcc = d(qy)/dx
      fl%mx_rhs =  fl%mx_rhs + fl%rre * dmdy_pcc * apcc
!----------------------------------------------------------------------------------------------------------
!   Y-pencil : X-mom diffusion term (x-v5/7), d(mu^x)/dy * d(qx)/dy at (i', j, k)
!----------------------------------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2C_3D(qx_ypencil, apcc_ypencil, dm, dm%ibcy(:, 1), dm%fbcy_var(:, :, :, 1) ) !apcc_ypencil = d(qx)/dy_ypencil
      mx_rhs_ypencil =  mx_rhs_ypencil + fl%rre * dmdy_pcc_ypencil * apcc_ypencil
!----------------------------------------------------------------------------------------------------------
!   X-pencil : X-mom diffusion term (x-v6/7), d(mu^x)/dz * d(qz^z)/dx at (i', j, k)
!----------------------------------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2P_3D(qz_ccc, apcc, dm, dm%ibcx(:, 3), dm%fbcx_var(:, :, :, 3) ) ! apcc = d(qz)/dx
      fl%mx_rhs =  fl%mx_rhs + fl%rre * dmdz_pcc * apcc
!----------------------------------------------------------------------------------------------------------
!   Z-pencil : X-mom diffusion term (x-v7/7), d(mu^x)/dz * d(qx)/dz at (i', j, k)
!----------------------------------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2C_3D(qx_zpencil, apcc_zpencil, dm, dm%ibcz(:, 1), dm%fbcz_var(:, :, :, 1) ) ! apcc_zpencil = d(qx)/dz
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
    i = 2
!----------------------------------------------------------------------------------------------------------
! X-pencil : Y-mom convection term (y-c1/3), d(gx^y * qy^x)/dx at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_x_1st_derivative_P2C_3D( -qx_ppc * qy_ppc, acpc, dm, dm%ibcx(:, i), dm%fbcx_var(:, :, :, i) * dm%fbcx_var(:, :, :, 1) )
    else
      call Get_x_1st_derivative_P2C_3D( -gx_ppc * qy_ppc, acpc, dm, dm%ibcx(:, i), dm%fbcx_var(:, :, :, i) * dm%fbcx_var(:, :, :, 1 + NBC) )
    end if
    fl%my_rhs = fl%my_rhs + acpc
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom convection term (y-c2/3), d(gy * qy)/dy at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_y_1st_derivative_C2P_3D(-qy_ccc_ypencil * qy_ccc_ypencil, acpc_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i) * dm%fbcy_var(:, :, :, 2) )
    else
      call Get_y_1st_derivative_C2P_3D(-gy_ccc_ypencil * qy_ccc_ypencil, acpc_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i) * dm%fbcy_var(:, :, :, 2 + NBC) )
    end if
    my_rhs_ypencil = my_rhs_ypencil + acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Y-mom convection term (y-c3/3), d(<gz>^y * <qy>^z)/dz at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_z_1st_derivative_P2C_3D( -qz_cpp_zpencil * qy_cpp_zpencil, acpc_zpencil, dm, dm%ibcz(:, i), dm%fbcz_var(:, :, :, i) * dm%fbcz_var(:, :, :, 3) )
    else
      call Get_z_1st_derivative_P2C_3D( -gz_cpp_zpencil * qy_cpp_zpencil, acpc_zpencil, dm, dm%ibcz(:, i), dm%fbcz_var(:, :, :, i) * dm%fbcz_var(:, :, :, 3 + NBC) )
    end if
    my_rhs_zpencil = my_rhs_zpencil + acpc_zpencil
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom pressure gradient in y direction, d(sigma_1 p)
!----------------------------------------------------------------------------------------------------------
    call Get_y_1st_derivative_C2P_3D( -pres_ypencil, acpc_ypencil, dm, dm%ibcy(:, 4), dm%fbcy_var(:, :, :, 4) )
    my_rhs_pfc_ypencil =  my_rhs_pfc_ypencil + acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! X-pencil : Y-mom diffusion term (y-v1-1/7), \mu * LL1(uy) at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    !call Get_x_2nd_derivative_C2C_3D(fl%qy, acpc, dm, dm%ibcx(:, 2) )
    appc = ZERO
    acpc = ZERO
    call Get_x_1st_derivative_C2P_3D(fl%qy, appc, dm, dm%ibcx(:, 2), dm%fbcx_var(:, :, :, i) )
    call Get_x_1st_derivative_P2C_3D(appc,  acpc, dm, dm%ibcx(:, 2), dm%fbcx_var(:, :, :, i) )
    if (.not. dm%is_thermo) then
      fl%my_rhs = fl%my_rhs +         fl%rre * acpc
    else
      fl%my_rhs = fl%my_rhs + m_cpc * fl%rre * acpc
    end if
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v1-2/7), \mu * LL2(uy) at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    !call Get_y_2nd_derivative_P2P_3D(qy_ypencil, acpc_ypencil, dm, dm%ibcy(:, 2))
    call Get_y_1st_derivative_P2C_3D(qy_ypencil,   accc_ypencil, dm, dm%ibcy(:, 2), dm%fbcy_var(:, :, :, i))
    call Get_y_1st_derivative_C2P_3D(accc_ypencil, acpc_ypencil, dm, dm%ibcy(:, 2), dm%fbcy_var(:, :, :, i))
    if ( .not. dm%is_thermo ) then
      my_rhs_ypencil = my_rhs_ypencil +                 fl%rre * acpc_ypencil
    else
      my_rhs_ypencil = my_rhs_ypencil + m_cpc_ypencil * fl%rre * acpc_ypencil
    end if
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Y-mom diffusion term (y-v1-3/7), \mu * LL3(uy) at (i, j', k)
!----------------------------------------------------------------------------------------------------------
    !call Get_z_2nd_derivative_C2C_3D(qy_zpencil, acpc_zpencil, dm, dm%ibcz(:, 2))
    call Get_z_1st_derivative_C2P_3D(qy_zpencil,   acpp_zpencil, dm, dm%ibcz(:, 2), dm%fbcz_var(:, :, :, i))
    call Get_z_1st_derivative_P2C_3D(acpp_zpencil, acpc_zpencil, dm, dm%ibcz(:, 2), dm%fbcz_var(:, :, :, i))
    if ( .not. dm%is_thermo ) my_rhs_zpencil = my_rhs_zpencil +                  fl%rre * acpc_zpencil
    if ( dm%is_thermo) my_rhs_zpencil = my_rhs_zpencil + m_cpc_zpencil * fl%rre * acpc_zpencil

    if ( dm%is_thermo) then
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom gravity force in y direction
!----------------------------------------------------------------------------------------------------------
      if(fl%igravity == i .or. fl%igravity == -i) then
        call Get_y_midp_C2P_3D(dDens_ypencil, acpc_ypencil, dm, dm%ibcy(:, 5), dm%ftpbcy_var(:, :, :)%d )
        my_rhs_pfc_ypencil =  my_rhs_pfc_ypencil + fl%fgravity(i) * acpc_ypencil
      end if
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v2/7), \mu^y * 1/3 * d (div)/dy at (i, j', k)
!----------------------------------------------------------------------------------------------------------
      fbcy = ZERO ! check
      call Get_y_1st_derivative_C2P_3D(div_ypencil, acpc_ypencil, dm, dm%ibcy(:, i), fbcx)
      my_rhs_ypencil = my_rhs_ypencil + one_third_rre * m_cpc_ypencil * acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v3/7), 2d\mu/dy * (-1/3 * div(u)) +  2d\mu/dy * dv/dy
!----------------------------------------------------------------------------------------------------------
      fbcy = ZERO ! check
      call Get_y_midp_C2P_3D (div_ypencil, acpc_ypencil, dm, dm%ibcy(:, 2), fbcy )
      my_rhs_ypencil = my_rhs_ypencil - two_third_rre * dmdy_cpc_ypencil * acpc_ypencil
      call Get_y_1st_derivative_P2P_3D(qy_ypencil,  acpc_ypencil, dm, dm%ibcy(:, 2), dm%fbcy_var(:, :, :, 2) )
      my_rhs_ypencil = my_rhs_ypencil + two_rre       * dmdy_cpc_ypencil * acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v4/7), d(mu^y)/dx * d(qx^x)/dy at (i, j', k)
!----------------------------------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2P_3D(qx_ccc_ypencil, acpc_ypencil, dm, dm%ibcy(:, 1), dm%fbcy_var(:, :, :, 1) )
      my_rhs_ypencil =  my_rhs_ypencil + fl%rre * dmdx_cpc_ypencil * acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! X-pencil : Y-mom diffusion term (y-v5/7), d(mu^y)/dx * d(qy^x))/dx at (i, j', k)
!----------------------------------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2C_3D(fl%qy, acpc, dm, dm%ibcx(:, 2), dm%fbcx_var(:, :, :, 2) )
      fl%my_rhs =  fl%my_rhs + fl%rre * dmdx_cpc * acpc
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Y-mom diffusion term (y-v6/7), d(mu^y)/dz * d(qz^z)/dy at (i, j', k)
!----------------------------------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2P_3D(qz_ccc_ypencil, acpc_ypencil, dm, dm%ibcy(:, 3), dm%fbcy_var(:, :, :, 3) )
      my_rhs_ypencil =  my_rhs_ypencil + fl%rre * dmdz_cpc_ypencil * acpc_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Y-mom diffusion term (y-v7/7), d(mu^y)/dz * d(qy)/dz at (i, j', k)
!----------------------------------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2C_3D(qy_zpencil, acpc_zpencil, dm, dm%ibcz(:, 2), dm%fbcz_var(:, :, :, 2) )
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
    i = 3
!----------------------------------------------------------------------------------------------------------
! X-pencil : Z-mom convection term (z-c1/3), d(gx^z * qz^x)/dx at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_x_1st_derivative_P2C_3D( -qx_pcp * qz_pcp, accp, dm, dm%ibcx(:, i), dm%fbcx_var(:, :, :, i) * dm%fbcx_var(:, :, :, 1)  )
    else
      call Get_x_1st_derivative_P2C_3D( -gx_pcp * qz_pcp, accp, dm, dm%ibcx(:, i), dm%fbcx_var(:, :, :, i) * dm%fbcx_var(:, :, :, 1 + NBC)  )
    end if
    fl%mz_rhs = fl%mz_rhs + accp
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Z-mom convection term (z-c2/3), d(gy^z * qz^y)/dy at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_y_1st_derivative_P2C_3D( -qy_cpp_ypencil * qz_cpp_ypencil, accp_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i) * dm%fbcy_var(:, :, :, 2) )
   else
      call Get_y_1st_derivative_P2C_3D( -gy_cpp_ypencil * qz_cpp_ypencil, accp_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i) * dm%fbcy_var(:, :, :, 2 + NBC) )
    end if
    mz_rhs_ypencil = mz_rhs_ypencil + accp_ypencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Z-mom convection term (y-c3/3), d(gz * qz)/dz at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    if ( .not. dm%is_thermo) then
      call Get_z_1st_derivative_C2P_3D(-qz_ccc_zpencil * qz_ccc_zpencil, accp_zpencil, dm, dm%ibcz(:, 3), dm%fbcz_var(:, :, :, i) * dm%fbcz_var(:, :, :, 3) )
    else
      call Get_z_1st_derivative_C2P_3D(-gz_ccc_zpencil * qz_ccc_zpencil, accp_zpencil, dm, dm%ibcz(:, 3), dm%fbcz_var(:, :, :, i) * dm%fbcz_var(:, :, :, 3 + NBC) )
    end if
    mz_rhs_zpencil = mz_rhs_zpencil + accp_zpencil

!----------------------------------------------------------------------------------------------------------
! z-pencil : pressure gradient in z direction, d(sigma_1 p)
!----------------------------------------------------------------------------------------------------------
    call Get_z_1st_derivative_C2P_3D( -pres_zpencil, accp_zpencil, dm, dm%ibcz(:, 4), dm%fbcz_var(:, :, :, 4) )
    mz_rhs_pfc_zpencil =  mz_rhs_pfc_zpencil + accp_zpencil
!----------------------------------------------------------------------------------------------------------
! X-pencil : Z-mom diffusion term (z-v1-1/7), \mu * L11(uz) at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    !call Get_x_2nd_derivative_C2C_3D(fl%qz, accp, dm, dm%ibcx(:, 3) )
    call Get_x_1st_derivative_C2P_3D(fl%qz, apcp, dm, dm%ibcx(:, i),  dm%fbcx_var(:, :, :, i))
    call Get_x_1st_derivative_P2C_3D(apcp,  accp, dm, dm%ibcx(:, i),  dm%fbcx_var(:, :, :, i))
    if ( .not. dm%is_thermo) then
      fl%mz_rhs = fl%mz_rhs +         fl%rre * accp
    else
      fl%mz_rhs = fl%mz_rhs + m_ccp * fl%rre * accp
    end if
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Z-mom diffusion term (z-v1-2/1), \mu * L22(uz) at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    !call Get_y_2nd_derivative_C2C_3D( qz_ypencil, accp_ypencil, dm, dm%ibcy(:, 3), dm%fbcy_var(:, :, :, 3))
    call Get_y_1st_derivative_C2P_3D( qz_ypencil,   acpp_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i))
    call Get_y_1st_derivative_P2C_3D( acpp_ypencil, accp_ypencil, dm, dm%ibcy(:, i), dm%fbcy_var(:, :, :, i))
    if ( .not. dm%is_thermo) then
      mz_rhs_ypencil = mz_rhs_ypencil +                 fl%rre * accp_ypencil
    else 
      mz_rhs_ypencil = mz_rhs_ypencil + m_ccp_ypencil * fl%rre * accp_ypencil
    end if
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v1-3/7), \mu * L33(uz) at (i, j, k')
!----------------------------------------------------------------------------------------------------------
    !call Get_z_2nd_derivative_P2P_3D(qz_zpencil, accp_zpencil, dm, dm%ibcz(:, 3))
    call Get_z_1st_derivative_P2C_3D(qz_zpencil,   accc_zpencil, dm, dm%ibcz(:, 3), dm%fbcz_var(:, :, :, 3))
    call Get_z_1st_derivative_C2P_3D(accc_zpencil, accp_zpencil, dm, dm%ibcz(:, 3), dm%fbcz_var(:, :, :, 3))
    if ( .not. dm%is_thermo) then
      mz_rhs_zpencil = mz_rhs_zpencil +                 fl%rre * accp_zpencil
    else
      mz_rhs_zpencil = mz_rhs_zpencil + m_ccp_zpencil * fl%rre * accp_zpencil
    end if

    if ( dm%is_thermo) then
!----------------------------------------------------------------------------------------------------------
! z-pencil : Z-mom gravity force in z direction
!----------------------------------------------------------------------------------------------------------
      if( fl%igravity == i .or. fl%igravity == -i) then
        call Get_z_midp_C2P_3D(dDens_zpencil, accp_zpencil, dm, dm%ibcz(:, 5), dm%ftpbcz_var(:, :, :)%d )
        mz_rhs_pfc_zpencil =  mz_rhs_pfc_zpencil + fl%fgravity(i) * accp_zpencil
      end if
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v2/7), \mu * 1/3 * d (div)/dz at (i, j, k')
!----------------------------------------------------------------------------------------------------------
      fbcz = ZERO
      call Get_z_1st_derivative_C2P_3D(div_zpencil, accp_zpencil, dm, dm%ibcz(:, 3), fbcz )
      mz_rhs_zpencil = mz_rhs_zpencil + one_third_rre * m_ccp_zpencil * accp_zpencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v3/7), 2d\mu/dz * (-1/3 * div(u)) +  2d\mu/dz * dw/dz
!----------------------------------------------------------------------------------------------------------
      fbcz = ZERO
      call Get_z_midp_C2P_3D(div_zpencil, accp_zpencil, dm,  dm%ibcz(:, 3), fbcz )
      mz_rhs_zpencil = mz_rhs_zpencil - two_third_rre * dmdz_ccp_zpencil * accp_zpencil
      call Get_z_1st_derivative_P2P_3D( qz_zpencil,  accp_zpencil, dm, dm%ibcz(:, 3), dm%fbcz_var(:, :, :, 3) )
      mz_rhs_zpencil = mz_rhs_zpencil + two_rre * dmdz_ccp_zpencil * accp_zpencil
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v4/7), d(mu^z)/dx * d(qx^x)/dz at (i, j, k')
!----------------------------------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2P_3D(qx_ccc_zpencil, accp_zpencil, dm, dm%ibcz(:, 1), dm%fbcz_var(:, :, :, 1) )
      mz_rhs_zpencil = mz_rhs_zpencil + fl%rre * dmdx_ccp_zpencil * accp_zpencil
!----------------------------------------------------------------------------------------------------------
! X-pencil : Z-mom diffusion term (z-v5/7), d(mu^z)/dx * d(qz)/dx at (i, j, k')
!----------------------------------------------------------------------------------------------------------
      call Get_x_1st_derivative_C2C_3D(fl%qz, accp, dm, dm%ibcx(:, 3), dm%fbcx_var(:, :, :, 3) )
      fl%mz_rhs =  fl%mz_rhs + fl%rre * dmdx_ccp * accp
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Z-mom diffusion term (z-v6/7), d(mu^z)/dy * d(qy^y)/dz at (i, j, k')
!----------------------------------------------------------------------------------------------------------
      call Get_z_1st_derivative_C2P_3D(qy_ccc_zpencil, accp_zpencil, dm, dm%ibcz(:, 2), dm%fbcz_var(:, :, :, 2) )
      mz_rhs_zpencil = mz_rhs_zpencil + fl%rre * dmdy_ccp_zpencil * accp_zpencil
!----------------------------------------------------------------------------------------------------------
! Y-pencil : Z-mom diffusion term (z-v7/7), d(mu^z)/dy * d(qz)/dy at (i, j, k')
!----------------------------------------------------------------------------------------------------------
      call Get_y_1st_derivative_C2C_3D(qz_ypencil, accp_ypencil, dm, dm%ibcy(:, 3), dm%fbcy_var(:, :, :, 3) )
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
    
#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_mid_msg("Updating the velocity/mass flux ...")
#endif
!----------------------------------------------------------------------------------------------------------
!   x-pencil, ux = ux - dt * alpha * d(phi_ccc)/dx
!----------------------------------------------------------------------------------------------------------
    dphidx_pcc = ZERO
    call Get_x_1st_derivative_C2P_3D(phi_ccc,  dphidx_pcc, dm, dm%ibcx(:, 4), dm%fbcx_var(:, :, :, 4) )
    ux = ux - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidx_pcc
!----------------------------------------------------------------------------------------------------------
!   y-pencil, uy = uy - dt * alpha * d(phi_ccc)/dy
!----------------------------------------------------------------------------------------------------------
    phi_ccc_ypencil = ZERO
    dphidy_cpc_ypencil = ZERO
    dphidy_cpc = ZERO
    call transpose_x_to_y (phi_ccc, phi_ccc_ypencil, dm%dccc)
    call Get_y_1st_derivative_C2P_3D(phi_ccc_ypencil, dphidy_cpc_ypencil, dm, dm%ibcy(:, 4), dm%fbcy_var(:, :, :, 4) )
    call transpose_y_to_x (dphidy_cpc_ypencil, dphidy_cpc, dm%dcpc)
    uy = uy - dm%dt * dm%tAlpha(isub) * dm%sigma2p * dphidy_cpc
!----------------------------------------------------------------------------------------------------------
!   z-pencil, uz = uz - dt * alpha * d(phi_ccc)/dz
!----------------------------------------------------------------------------------------------------------
    pphi_ccc_zpencil =  ZERO
    dphidz_ccp_zpencil = ZERO
    dphidz_ccp_ypencil = ZERO
    dphidz_ccp = ZERO
    call transpose_y_to_z (phi_ccc_ypencil, pphi_ccc_zpencil, dm%dccc)
    call Get_z_1st_derivative_C2P_3D(pphi_ccc_zpencil, dphidz_ccp_zpencil, dm, dm%ibcz(:, 4), dm%fbcz_var(:, :, :, 4) )
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
    !integer :: i, j, k, jj, ii


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
    !integer :: i, j, k, jj, ii

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
      call Calculate_drhodt(dm, tm%dDens, tm%dDensm1, tm%dDensm2, rhs)
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
    use wtformat_mod
#endif
    implicit none

    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(inout) :: dm
    integer,        intent(in)    :: isub

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

    !in order for a high order spacial accuracy
    ! to use Alternating direction implicit method
    ! ref: Cui2013: Convergence analysis of high-order compact 
    ! alternating direction implicit schemes for the two-dimensional 
    ! time fractional equation

#ifdef DEBUG_STEPS
  call write_snapshot_any3darray(fl%qx, 'qxs_RK'//trim(int2str(isub)), 'debug', dm%dpcc, dm, fl%iteration)
  call write_snapshot_any3darray(fl%qy, 'qys_RK'//trim(int2str(isub)), 'debug', dm%dcpc, dm, fl%iteration)
  call write_snapshot_any3darray(fl%qz, 'qzs_RK'//trim(int2str(isub)), 'debug', dm%dccp, dm, fl%iteration)
#endif
!----------------------------------------------------------------------------------------------------------
! to solve Poisson equation
!----------------------------------------------------------------------------------------------------------
    !if(nrank == 0) call Print_debug_mid_msg("  Solving Poisson Equation ...") 
    !call solve_poisson_x2z(fl, dm, isub) !
    call solve_poisson(fl, dm, isub) ! test show above two methods gave the same results. 
#ifdef DEBUG_STEPS
    call write_snapshot_any3darray(fl%pcor, 'pcor'//trim(int2str(isub)), 'debug', dm%dccc, dm, fl%iteration)
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
!----------------------------------------------------------------------------------------------------------
! to update pressure
!----------------------------------------------------------------------------------------------------------
#ifdef DEBUG_STEPS  
    if(nrank == 0) &
    call Print_debug_mid_msg("Correcting the pressure term ...")
#endif
    fl%pres(:, :, :) = fl%pres(:, :, :) + fl%pcor(:, :, :)

!----------------------------------------------------------------------------------------------------------
! to update velocity from gx gy gz 
!----------------------------------------------------------------------------------------------------------
  if(dm%is_thermo) then
    call Calculate_velocity_from_massflux(dm, fl)
    call update_gxgygz_bc(dm)
  end if

    
#ifdef DEBUG_STEPS
    call Find_maximum_velocity(dm, fl%qx, fl%qy, fl%qz)
    call Check_mass_conservation(fl, dm, "isub"//trim(int2str(isub))) 
#endif

    return
  end subroutine Solve_momentum_eq

end module eq_momentum_mod
