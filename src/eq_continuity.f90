module continuity_eq_mod
  use operations
  use decomp_2d
  public :: Calculate_drhodt
  public :: Get_divergence_vel
  public :: Get_divergence_vel_x2z
  public  :: Check_mass_conservation
contains
!==========================================================================================================
!==========================================================================================================
!> \brief To calculate d(\rho)/dt in the continuity eq.
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
!> \param[in]     ux           ux or gx
!> \param[in]     uy           uy or gy
!> \param[in]     uz           uz or gz
!> \param[out]    div          div(u) or div(g)
!> \param[in]     d            domain
!_______________________________________________________________________________
  subroutine Get_divergence_vel(ux, uy, uz, div, dm)
    use parameters_constant_mod
    use udf_type_mod
#ifdef DEBUG
    use typeconvert_mod
#endif
    implicit none

    type(t_domain), intent (in) :: dm
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent (in ) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent (in ) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent (in ) :: uz
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent (out) :: div

    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div0
    real(WP), dimension(dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3)) :: div0_ypencil
    real(WP), dimension(dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3)) :: div0_zpencil

    real(WP), dimension(dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3)) :: uy_ypencil
    real(WP), dimension(dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3)) :: uz_ypencil
    real(WP), dimension(dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3)) :: uz_zpencil

    integer :: i, j, k
#ifdef DEBUG
    integer :: jj
    type(DECOMP_INFO) :: dtmp
#endif
    div = ZERO
!----------------------------------------------------------------------------------------------------------
! operation in x pencil, du/dx
!----------------------------------------------------------------------------------------------------------
    div0 = ZERO
    call Get_x_1st_derivative_P2C_3D(ux, div0, dm, dm%ibcx(:, 1))
    div(:, :, :) = div(:, :, :) + div0(:, :, :)

#ifdef DEBUG

    dtmp = dm%dccc
    k = 2
    i = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      open(121, file = 'debugy_div_x'//trim(int2str(nrank))//'.dat', position="append")
      do j = 1, dtmp%xsz(2)
        jj =dtmp%xst(2) + j - 1
        write(121, *) jj, ux(i, j, k), div0(i, j, k), div(i, j, k)
      end do
    end if

    k = 2
    j = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
        open(221, file = 'debugx_div_x'//trim(int2str(nrank))//'.dat', position="append")
        do i = 1, dtmp%xsz(1)
          write(221, *) i, ux(i, j, k), div0(i, j, k), div(i, j, k)
        end do
      end if
    end if

    i = 2
    j = 2
    if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
      open(321, file = 'debugz_div_x'//trim(int2str(nrank))//'.dat', position="append")
      do k = 1, dtmp%xsz(3)
        write(321, *) k, ux(i, j, k), div0(i, j, k), div(i, j, k)
      end do
    end if
#endif
!----------------------------------------------------------------------------------------------------------
! operation in y pencil, dv/dy
!----------------------------------------------------------------------------------------------------------
    uy_ypencil = ZERO
    div0_ypencil = ZERO
    div0 = ZERO
    call transpose_x_to_y(uy, uy_ypencil, dm%dcpc)
    call Get_y_1st_derivative_P2C_3D(uy_ypencil, div0_ypencil, dm, dm%ibcy(:, 2))
    call transpose_y_to_x(div0_ypencil, div0, dm%dccc)
    div(:, :, :) = div(:, :, :) + div0(:, :, :)

#ifdef DEBUG

    dtmp = dm%dccc
    k = 2
    i = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      open(122, file = 'debugy_div_y'//trim(int2str(nrank))//'.dat', position="append")
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        write(122, *) jj, uy(i, j, k), div0(i, j, k), div(i, j, k)
      end do
    end if

    k = 2
    j = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
        open(222, file = 'debugx_div_y'//trim(int2str(nrank))//'.dat', position="append")
        do i = 1, dtmp%xsz(1)
          write(222, *) i, uy(i, j, k), div0(i, j, k), div(i, j, k)
        end do
      end if
    end if

    i = 2
    j = 2
    if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
      open(323, file = 'debugz_div_y'//trim(int2str(nrank))//'.dat', position="append")
      do k = 1, dtmp%xsz(3)
        write(323, *) k, uy(i, j, k), div0(i, j, k), div(i, j, k)
      end do
    end if
#endif
!----------------------------------------------------------------------------------------------------------
! operation in z pencil, dw/dz
!----------------------------------------------------------------------------------------------------------
    uz_ypencil = ZERO
    uz_zpencil = ZERO
    div0_zpencil = ZERO
    div0_ypencil = ZERO
    div0 = ZERO
    call transpose_x_to_y(uz,         uz_ypencil, dm%dccp)
    call transpose_y_to_z(uz_ypencil, uz_zpencil, dm%dccp)
    call Get_z_1st_derivative_P2C_3D(uz_zpencil, div0_zpencil, dm, dm%ibcz(:, 3))
    call transpose_z_to_y(div0_zpencil, div0_ypencil, dm%dccc)
    call transpose_y_to_x(div0_ypencil, div0,         dm%dccc)
    div(:, :, :) = div(:, :, :) + div0(:, :, :)
#ifdef DEBUG

    dtmp = dm%dccc
    k = 2
    i = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      open(123, file = 'debugy_div_z'//trim(int2str(nrank))//'.dat', position="append")
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        write(123, *) jj, uz(i, j, k), div0(i, j, k), div(i, j, k)
      end do
    end if

    k = 2
    j = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
        open(223, file = 'debugx_div_z'//trim(int2str(nrank))//'.dat', position="append")
        do i = 1, dtmp%xsz(1)
          write(223, *) i, uz(i, j, k), div0(i, j, k), div(i, j, k)
        end do
      end if
    end if

    j = 2
    i = 2
    if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
      open(323, file = 'debugz_div_z'//trim(int2str(nrank))//'.dat', position="append")
      do k = 1, dtmp%xsz(3)
        write(323, *) k, uz(i, j, k), div0(i, j, k), div(i, j, k)
      end do
    end if

#endif
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
  subroutine Get_divergence_vel_x2z(ux, uy, uz, div_zpencil_ggg, dm)
    use parameters_constant_mod
    use udf_type_mod
    use decomp_extended_mod
#ifdef DEBUG
    use typeconvert_mod
#endif
    implicit none

    type(t_domain), intent (in) :: dm
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent (in) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent (in) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent (in) :: uz
    real(WP), dimension(dm%dccc%zst(1) : dm%dccc%zen(1), &
                        dm%dccc%zst(2) : dm%dccc%zen(2), &
                        dm%dccc%zst(3) : dm%dccc%zen(3)), intent (out) :: div_zpencil_ggg

    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div0
    real(WP), dimension(dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3)) :: div0_ypencil
    real(WP), dimension(dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3)) :: div0_zpencil
    real(WP), dimension(dm%dccc%yst(1) : dm%dccc%yen(1), &
                        dm%dccc%yst(2) : dm%dccc%yen(2), &
                        dm%dccc%ysz(3))                  :: div0_ypencil_ggl
    real(WP), dimension(dm%dccc%yst(1) : dm%dccc%yen(1), &
                        dm%dccc%yst(2) : dm%dccc%yen(2), &
                        dm%dccc%ysz(3))                  :: div_ypencil_ggl
    real(WP), dimension(dm%dccc%zst(1) : dm%dccc%zen(1), &
                        dm%dccc%zst(2) : dm%dccc%zen(2), &
                        dm%dccc%zst(3) : dm%dccc%zen(3)) :: div0_zpencil_ggg

    real(WP), dimension(dm%dcpc%ysz(1),                  dm%dcpc%ysz(2), dm%dcpc%ysz(3)) :: uy_ypencil
    real(WP), dimension(dm%dcpc%yst(1) : dm%dcpc%yen(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3)) :: uy_ypencil_ggl
    real(WP), dimension(dm%dccp%yst(1) : dm%dccp%yen(1), dm%dccp%ysz(2), dm%dccp%ysz(3)) :: uz_ypencil_ggl

    real(WP), dimension(dm%dccp%ysz(1),                  dm%dccp%ysz(2),                  dm%dccp%ysz(3)) :: uz_ypencil
    real(WP), dimension(dm%dccp%zsz(1),                  dm%dccp%zsz(2),                  dm%dccp%zsz(3)) :: uz_zpencil
    real(WP), dimension(dm%dccp%zst(1) : dm%dccp%zen(1), dm%dccp%zst(2) : dm%dccp%zen(2), dm%dccp%zsz(3)) :: uz_zpencil_ggl

#ifdef DEBUG
    integer :: jj, i, j ,k
    type(DECOMP_INFO) :: dtmp
#endif
!----------------------------------------------------------------------------------------------------------
! operation in x pencil, du/dx
!----------------------------------------------------------------------------------------------------------
    div0 = ZERO
    div0_ypencil_ggl = ZERO
    div_ypencil_ggl = ZERO
    call Get_x_1st_derivative_P2C_3D(ux, div0, dm, dm%ibcx(:, 1))
    call transpose_x_to_y(div0, div0_ypencil_ggl, dm%dccc)
    div_ypencil_ggl = div0_ypencil_ggl

#ifdef DEBUG

    dtmp = dm%dccc
    k = 2
    i = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      open(121, file = 'debugy_div_x'//trim(int2str(nrank))//'.dat', position="append")
      do j = 1, dtmp%xsz(2)
        jj =dtmp%xst(2) + j - 1
        write(121, *) jj, ux(i, j, k), div0(i, j, k), div0_ypencil_ggl(i, jj, k)
      end do
    end if

    k = 2
    j = 2
    jj =dtmp%xst(2) + j - 1
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
        open(221, file = 'debugx_div_x'//trim(int2str(nrank))//'.dat', position="append")
        do i = 1, dtmp%xsz(1)
          write(221, *) i, ux(i, j, k), div0(i, j, k), div0_ypencil_ggl(i, jj, k)
        end do
      end if
    end if

    i = 2
    j = 2
    jj =dtmp%xst(2) + j - 1
    if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
      open(321, file = 'debugz_div_x'//trim(int2str(nrank))//'.dat', position="append")
      do k = 1, dtmp%xsz(3)
        write(321, *) k, ux(i, j, k), div0(i, j, k), div0_ypencil_ggl(i, jj, k)
      end do
    end if
#endif
!----------------------------------------------------------------------------------------------------------
! operation in y pencil, dv/dy
!----------------------------------------------------------------------------------------------------------
    uy_ypencil = ZERO
    div0_ypencil = ZERO
    div0_ypencil_ggl = ZERO
    call transpose_x_to_y(uy, uy_ypencil, dm%dcpc)
    call Get_y_1st_derivative_P2C_3D(uy_ypencil, div0_ypencil, dm, dm%ibcy(:, 2))
    call ypencil_index_lgl2ggl(div0_ypencil, div0_ypencil_ggl, dm%dccc)
    div_ypencil_ggl = div_ypencil_ggl + div0_ypencil_ggl
    call transpose_y_to_z(div_ypencil_ggl, div_zpencil_ggg, dm%dccc)

#ifdef DEBUG

    dtmp = dm%dccc
    k = 2
    i = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      open(122, file = 'debugy_div_y'//trim(int2str(nrank))//'.dat', position="append")
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        write(122, *) jj, uy(i, j, k), div0_ypencil(i, j, k), div_ypencil_ggl(i, jj, k)
      end do
    end if

    k = 2
    j = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
        open(222, file = 'debugx_div_y'//trim(int2str(nrank))//'.dat', position="append")
        do i = 1, dtmp%xsz(1)
          write(222, *) i, uy(i, j, k), div0_ypencil(i, j, k), div_ypencil_ggl(i, jj, k)
        end do
      end if
    end if

    i = 2
    j = 2
    if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
      open(323, file = 'debugz_div_y'//trim(int2str(nrank))//'.dat', position="append")
      do k = 1, dtmp%xsz(3)
        write(323, *) k, uy(i, j, k), div0_ypencil(i, j, k), div_ypencil_ggl(i, jj, k)
      end do
    end if
#endif

!----------------------------------------------------------------------------------------------------------
! operation in z pencil, dw/dz
!----------------------------------------------------------------------------------------------------------
    uz_ypencil = ZERO
    uz_zpencil = ZERO
    div0_zpencil = ZERO
    div0_zpencil_ggg = ZERO
    call transpose_x_to_y(uz,         uz_ypencil, dm%dccp)
    call transpose_y_to_z(uz_ypencil, uz_zpencil, dm%dccp)
    call Get_z_1st_derivative_P2C_3D(uz_zpencil, div0_zpencil, dm, dm%ibcz(:, 3))
    call zpencil_index_llg2ggg(div0_zpencil, div0_zpencil_ggg, dm%dccc)
    div_zpencil_ggg = div_zpencil_ggg + div0_zpencil_ggg

#ifdef DEBUG

    dtmp = dm%dccc
    k = 2
    i = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      open(123, file = 'debugy_div_z'//trim(int2str(nrank))//'.dat', position="append")
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        write(123, *) jj, uz(i, j, k), div0_zpencil(i, j, k), div_zpencil_ggg(i, j, k)
      end do
    end if

    k = 2
    j = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
        open(223, file = 'debugx_div_z'//trim(int2str(nrank))//'.dat', position="append")
        do i = 1, dtmp%xsz(1)
          write(223, *) i, uz(i, j, k), div0_zpencil(i, j, k), div_zpencil_ggg(i, j, k)
        end do
      end if
    end if

    j = 2
    i = 2
    if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
      open(323, file = 'debugz_div_z'//trim(int2str(nrank))//'.dat', position="append")
      do k = 1, dtmp%xsz(3)
        write(323, *) k, uz(i, j, k), div0_zpencil(i, j, k), div_zpencil_ggg(i, j, k)
      end do
    end if

#endif

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
  subroutine Check_mass_conservation(fl, dm)
    use precision_mod
    use udf_type_mod
    use input_general_mod    
    use parameters_constant_mod
    use math_mod                
    use mpi_mod
    use solver_tools_mod
    use wtformat_mod
#ifdef DEBUG
    use visulisation_mod
#endif
    implicit none

    type(t_domain), intent( in    ) :: dm
    type(t_flow),   intent( inout ) :: fl                  

    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: div
    real(WP)   :: divmax 

    fl%pcor = ZERO
    div(:, :, :)  = ZERO
!----------------------------------------------------------------------------------------------------------
! $d\rho / dt$ at cell centre
!_______________________________________________________________________________
    if (dm%ithermo == 1) then
      call Calculate_drhodt(dm, fl%dDens, fl%dDensm1, fl%dDensm2, fl%pcor)
    end if
!----------------------------------------------------------------------------------------------------------
! $d(\rho u_i)) / dx_i $ at cell centre
!_______________________________________________________________________________
    if (dm%ithermo == 1) then
      call Get_divergence_vel(fl%gx, fl%gy, fl%gz, div, dm)
    else
      call Get_divergence_vel(fl%qx, fl%qy, fl%qz, div, dm)
    end if

#ifdef DEBUG
    call view_data_in_rank(div,   dm%dccc, dm, 'div', 0)
#endif
    
    call Find_maximum_absvar3d(div, "Check Mass Conservation:")

    ! if(nrank == 0) then
    !   write (*, wrtfmt1e) "  Check Mass Conservation:", divmax
    ! end if

    return
  end subroutine Check_mass_conservation 

end module continuity_eq_mod
