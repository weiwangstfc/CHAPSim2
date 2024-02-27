module continuity_eq_mod
  use operations
  use decomp_2d
  public :: Calculate_drhodt
  public :: Get_divergence
  !public :: Get_divergence_x2z
  public :: Check_mass_conservation
contains
!==========================================================================================================
!==========================================================================================================
!> \brief To calculate d(\rho)/dt in the continuity eq.
!? to do, to check, same as CHAPSim1 method2, but different from CHAPSim1 method 3.
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  dDens            density at the current time step
!> \param[in]  dDensm1          density at the t-1 time step
!> \param[in]  dDensm2          density at the t-2 time step
!> \param[out] drhodt           d(rho)/dt
!> \param[in]  itime            the sub-step in RK3
!_______________________________________________________________________________
  subroutine Calculate_drhodt(dm, dDens, dDensm1, dDensm2, drhodt)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in) :: dm
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent ( in  ) :: dDens, dDensm1, dDensm2
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent ( out ) :: drhodt

    integer :: i

    if(dm%iTimeScheme == ITIME_AB2) then

      drhodt(:, :, :) = HALF * dDens  (:, :, :) - &
                        TWO  * dDensm1(:, :, :) + &
                        HALF * dDensm2(:, :, :)
      drhodt(:, :, :) = drhodt(:, :, :) * dm%dt

    else if (dm%iTimeScheme == ITIME_RK3 .or. dm%iTimeScheme == ITIME_RK3_CN) then

      ! to check this part, is iteration necessary?
      drhodt(:, :, :) = dDens  (:, :, :)
      do i = 1, dm%nsubitr
        drhodt(:, :, :) = drhodt(:, :, :) + dm%tAlpha(i) * &
                          (dDensm1(:, :, :) - dDensm2(:, :, :))  * dm%dt
      end do

    else  

      ! default, Euler 1st order 
      drhodt(:, :, :) = dDens(:, :, :) - dDensm1(:, :, :)
      drhodt(:, :, :) = drhodt(:, :, :) * dm%dt

    end if


    return
  end subroutine Calculate_drhodt

!==========================================================================================================
!==========================================================================================================
!> \brief To calculate divergence of (rho * u) or divergence of (u)
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    div          div(q) or div(g)
!> \param[in]     d            domain
!_______________________________________________________________________________
  subroutine Get_divergence(fl, div, dm, flg)
    use parameters_constant_mod
    use udf_type_mod
    implicit none

    type(t_domain), intent(in) :: dm
    type(t_flow),    intent(in) :: fl
    character(1), optional, intent(in) :: flg
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent (out) :: div

    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)):: qx
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)):: qy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)):: qz

    real(WP), dimension(4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)) :: fbcx_qx
    real(WP), dimension(dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3)) :: fbcy_qy
    real(WP), dimension(dm%dccp%zsz(1), dm%dccp%zsz(2), 4) :: fbcz_qz
    

    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div0
    real(WP), dimension(dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3)) :: div0_ypencil
    real(WP), dimension(dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3)) :: div0_zpencil

    real(WP), dimension(dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3)) :: qy_ypencil
    real(WP), dimension(dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3)) :: qz_ypencil
    real(WP), dimension(dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3)) :: qz_zpencil

    if(dm%is_thermo) then
    
      if(flg == 'g')then
        qx = fl%gx
        qy = fl%gy
        qz = fl%gz
        fbcx_qx = fl%fbcx_gx
        fbcy_qy = fl%fbcy_gy
        fbcz_qz = fl%fbcz_gz
      else if (flg == 'q')then
        qx = fl%qx
        qy = fl%qy
        qz = fl%qz
        fbcx_qx = fl%fbcx_qx
        fbcy_qy = fl%fbcy_qy
        fbcz_qz = fl%fbcz_qz
      else ! default is 'gx' based
        qx = fl%gx
        qy = fl%gy
        qz = fl%gz
        fbcx_qx = fl%fbcx_gx
        fbcy_qy = fl%fbcy_gy
        fbcz_qz = fl%fbcz_gz
      end if

    else 
      qx = fl%qx
      qy = fl%qy
      qz = fl%qz
      fbcx_qx = fl%fbcx_qx
      fbcy_qy = fl%fbcy_qy
      fbcz_qz = fl%fbcz_qz
    end if

    div = ZERO
!----------------------------------------------------------------------------------------------------------
! operation in x pencil, dqx/dx
!----------------------------------------------------------------------------------------------------------
    div0 = ZERO
    call Get_x_1st_derivative_P2C_3D(qx, div0, dm, dm%ibcx(:, 1), fbcx_qx)
    div(:, :, :) = div(:, :, :) + div0(:, :, :)
!----------------------------------------------------------------------------------------------------------
! operation in y pencil, dqy/dy * (1/r)
!----------------------------------------------------------------------------------------------------------
    qy_ypencil = ZERO
    div0_ypencil = ZERO
    div0 = ZERO
    call transpose_x_to_y(qy, qy_ypencil, dm%dcpc)
    call Get_y_1st_derivative_P2C_3D(qy_ypencil, div0_ypencil, dm, dm%ibcy(:, 2), fbcy_qy)
    call transpose_y_to_x(div0_ypencil, div0, dm%dccc)
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(div0, dm%dccc, dm%rci, 1, IPENCIL(1))
    div(:, :, :) = div(:, :, :) + div0(:, :, :)
!----------------------------------------------------------------------------------------------------------
! operation in z pencil, dw/dz * (1/r)^2
!----------------------------------------------------------------------------------------------------------
    qz_ypencil = ZERO
    qz_zpencil = ZERO
    div0_zpencil = ZERO
    div0_ypencil = ZERO
    div0 = ZERO
    call transpose_x_to_y(qz,         qz_ypencil, dm%dccp)
    call transpose_y_to_z(qz_ypencil, qz_zpencil, dm%dccp)
    call Get_z_1st_derivative_P2C_3D(qz_zpencil, div0_zpencil, dm, dm%ibcz(:, 3), fbcz_qz)
    call transpose_z_to_y(div0_zpencil, div0_ypencil, dm%dccc)
    call transpose_y_to_x(div0_ypencil, div0,         dm%dccc)
    if(dm%icoordinate == ICYLINDRICAL) call multiple_cylindrical_rn(div0, dm%dccc, dm%rci, 2, IPENCIL(1))
    div(:, :, :) = div(:, :, :) + div0(:, :, :)
    
    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
!> \brief To calculate divergence of (rho * u) or divergence of (u)
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ux           ux or gx
!> \param[in]     uy           uy or gy
!> \param[in]     uz           uz or gz
!> \param[out]    div          div(u) or div(g)
!> \param[in]     d            domain
!_______________________________________________________________________________
!   subroutine Get_divergence_x2z(fl, div_zpencil_ggg, dm, flg)
!     use parameters_constant_mod
!     use udf_type_mod
!     use decomp_extended_mod

!     implicit none

!     type(t_domain), intent(in) :: dm
!     type(f_domain), intent(in) :: fl
!     character(1),   intent(in) :: flg
!     real(WP), dimension(dm%dccc%zst(1) : dm%dccc%zen(1), &
!                         dm%dccc%zst(2) : dm%dccc%zen(2), &
!                         dm%dccc%zst(3) : dm%dccc%zen(3)), intent (out) :: div_zpencil_ggg

!     real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)) :: qx
!     real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)) :: qy
!     real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)) :: qz
    

!     real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div0
!     real(WP), dimension(dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3)) :: div0_ypencil
!     real(WP), dimension(dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3)) :: div0_zpencil
!     real(WP), dimension(dm%dccc%yst(1) : dm%dccc%yen(1), &
!                         dm%dccc%yst(2) : dm%dccc%yen(2), &
!                         dm%dccc%ysz(3))                  :: div0_ypencil_ggl
!     real(WP), dimension(dm%dccc%yst(1) : dm%dccc%yen(1), &
!                         dm%dccc%yst(2) : dm%dccc%yen(2), &
!                         dm%dccc%ysz(3))                  :: div_ypencil_ggl
!     real(WP), dimension(dm%dccc%zst(1) : dm%dccc%zen(1), &
!                         dm%dccc%zst(2) : dm%dccc%zen(2), &
!                         dm%dccc%zst(3) : dm%dccc%zen(3)) :: div0_zpencil_ggg

!     real(WP), dimension(dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3)) :: qy_ypencil
!     real(WP), dimension(dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3)) :: qz_ypencil
!     real(WP), dimension(dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3)) :: qz_zpencil

!     if(dm%is_thermo) then
    
!       if(flg == 'g')then
!         qx = fl%gx
!         qy = fl%gy
!         qz = fl%gz
!         fbcx_qx = fl%fbcx_gx
!         fbcy_qy = fl%fbcy_gy
!         fbcz_qz = fl%fbcz_gz
!       else if (flg == 'q')then
!         qx = fl%qx
!         qy = fl%qy
!         qz = fl%qz
!         fbcx_qx = fl%fbcx_qx
!         fbcy_qy = fl%fbcy_qy
!         fbcz_qz = fl%fbcz_qz
!       else ! default is 'gx' based
!         qx = fl%gx
!         qy = fl%gy
!         qz = fl%gz
!         fbcx_qx = fl%fbcx_gx
!         fbcy_qy = fl%fbcy_gy
!         fbcz_qz = fl%fbcz_gz
!       end if

!     else 
!       qx = fl%qx
!       qy = fl%qy
!       qz = fl%qz
!       fbcx_qx = fl%fbcx_qx
!       fbcy_qy = fl%fbcy_qy
!       fbcz_qz = fl%fbcz_qz
!     end if

! !----------------------------------------------------------------------------------------------------------
! ! operation in x pencil, dqx/dx
! !----------------------------------------------------------------------------------------------------------
!     div0 = ZERO
!     div0_ypencil_ggl = ZERO
!     div_ypencil_ggl = ZERO
!     call Get_x_1st_derivative_P2C_3D(qx, div0, dm, dm%ibcx(:, 1), fbcx_qx)
!     call transpose_x_to_y(div0, div0_ypencil_lgl, dm%dccc)
!     call ypencil_index_lgl2ggl(div0_ypencil_lgl, div0_ypencil_ggl, dm%dccc)
!     div_ypencil_ggl = div0_ypencil_ggl
! !----------------------------------------------------------------------------------------------------------
! ! operation in y pencil, dqy/dy
! !----------------------------------------------------------------------------------------------------------
!     qy_ypencil = ZERO
!     div0_ypencil = ZERO
!     div0_ypencil_ggl = ZERO
!     call transpose_x_to_y(qy, qy_ypencil, dm%dcpc)
!     call Get_y_1st_derivative_P2C_3D(qy_ypencil, div0_ypencil, dm, dm%ibcy(:, 2), fbcy_qy)
!     call ypencil_index_lgl2ggl(div0_ypencil, div0_ypencil_ggl, dm%dccc)
!     div_ypencil_ggl = div_ypencil_ggl + div0_ypencil_ggl
!     call transpose_y_to_z(div_ypencil_ggl, div_zpencil_ggg, dm%dccc)
! !----------------------------------------------------------------------------------------------------------
! ! operation in z pencil, dqz/dz
! !----------------------------------------------------------------------------------------------------------
!     qz_ypencil = ZERO
!     qz_zpencil = ZERO
!     div0_zpencil = ZERO
!     div0_zpencil_ggg = ZERO
!     call transpose_x_to_y(qz,         qz_ypencil, dm%dccp)
!     call transpose_y_to_z(qz_ypencil, qz_zpencil, dm%dccp)
!     call Get_z_1st_derivative_P2C_3D(qz_zpencil, div0_zpencil, dm, dm%ibcz(:, 3), fbcz_qz)
!     call zpencil_index_llg2ggg(div0_zpencil, div0_zpencil_ggg, dm%dccc)
!     div_zpencil_ggg = div_zpencil_ggg + div0_zpencil_ggg

!     return
!   end subroutine

!==========================================================================================================
!==========================================================================================================
!> \brief To calculate divergence of (rho * u) or divergence of (u)
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ux           ux or gx
!> \param[in]     uy           uy or gy
!> \param[in]     uz           uz or gz
!> \param[out]    div          div(u) or div(g)
!> \param[in]     d            domain
!_______________________________________________________________________________
  subroutine Check_mass_conservation(fl, dm, str0, tm)
    use precision_mod
    use udf_type_mod
    use input_general_mod    
    use parameters_constant_mod
    use math_mod                
    use mpi_mod
    use solver_tools_mod
    use wtformat_mod
    use io_visulisation_mod

    implicit none

    type(t_domain), intent( in    ) :: dm
    type(t_flow),   intent( inout ) :: fl  
    character(*), intent(in), optional :: str0   
    type(t_thermo), intent(in), optional :: tm             

    character(32) :: str

    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div
    !real(WP)   :: divmax 

    if(present(str0)) then
      str = trim(str0)
    else
      str = ''
    end if

    fl%pcor = ZERO
    div(:, :, :)  = ZERO
!----------------------------------------------------------------------------------------------------------
! $d\rho / dt$ at cell centre
!----------------------------------------------------------------------------------------------------------
    if (dm%is_thermo) then
      call Calculate_drhodt(dm, tm%dDens, tm%dDensm1, tm%dDensm2, fl%pcor)
    end if
!----------------------------------------------------------------------------------------------------------
! $d(\rho u_i)) / dx_i $ at cell centre
!----------------------------------------------------------------------------------------------------------
    call Get_divergence(fl, div, dm)

#ifdef DEBUG_STEPS
    call write_snapshot_any3darray(div, 'divU', trim(str), dm%dccc, dm, fl%iteration)
#endif

    call Find_maximum_absvar3d(div, trim(str)//" Check Mass Conservation:", wrtfmt1e)

    ! if(nrank == 0) then
    !   write (*, wrtfmt1e) "  Check Mass Conservation:", divmax
    ! end if

    return
  end subroutine Check_mass_conservation 

end module continuity_eq_mod
