module solver_tools_mod

  ! procedure
  private

 
  public  :: Check_cfl_convection
  public  :: Check_cfl_diffusion

  public  :: Update_Re
  public  :: Update_PrGr
  public  :: Calculate_xz_mean_yprofile
  public  :: Adjust_to_xzmean_zero
  public  :: Get_volumetric_average_3d

  

contains
!==========================================================================================================
!> \brief The main code for initialising flow variables
!> This subroutine is called once in \ref initialise_chapsim.
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!----------------------------------------------------------------------------------------------------------
!> \param[inout]  
!==========================================================================================================
  subroutine Update_Re(iter, fl)
    use parameters_constant_mod
    use thermo_info_mod
    use udf_type_mod
    implicit none
    integer,     intent(in ) :: iter  
    type(t_flow),   intent(inout) :: fl
  !----------------------------------------------------------------------------------------------------------
  !  1/Re                                   
  !----------------------------------------------------------------------------------------------------------
    if(iter < fl%initReTo) then
      fl%rre = ONE / fl%reninit
    else
      fl%rre = ONE / fl%ren
    end if

    return
  end subroutine Update_Re

  subroutine Update_PrGr(fl, tm)
    use parameters_constant_mod
    use thermo_info_mod
    use udf_type_mod
    implicit none
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm
    

    real(WP) :: u0, rtmp
  
!----------------------------------------------------------------------------------------------------------
!  1/(Re*Pr)                                   
!----------------------------------------------------------------------------------------------------------
    tm%rPrRen = fl%rre * fluidparam%ftp0ref%k / fluidparam%ftp0ref%m / fluidparam%ftp0ref%cp
!----------------------------------------------------------------------------------------------------------
!  gravity force                          
!----------------------------------------------------------------------------------------------------------  
    u0 = ONE / fl%rre * fluidparam%ftp0ref%m / fluidparam%ftp0ref%d / tm%ref_l0
    rtmp = tm%ref_l0 / u0 / u0 * GRAVITY
    fl%fgravity = ZERO
    if (fl%igravity == 1 ) then ! flow/gravity same dirction - x
      fl%fgravity(1) =  rtmp
    else if (fl%igravity == 2 ) then ! flow/gravity same dirction - y
      fl%fgravity(2) =  rtmp
    else if (fl%igravity == 3 ) then ! flow/gravity same dirction - z
      fl%fgravity(3) =  rtmp
    else if (fl%igravity == -1 ) then ! flow/gravity opposite dirction - x
      fl%fgravity(1) =  - rtmp
    else if (fl%igravity == -2 ) then ! flow/gravity opposite dirction - y
      fl%fgravity(2) =  - rtmp
    else if (fl%igravity == -3 ) then ! flow/gravity opposite dirction - z
      fl%fgravity(3) =  - rtmp
    else ! no gravity
      fl%fgravity = ZERO
    end if
    
    return
  end subroutine Update_PrGr
!==========================================================================================================
!> \brief The main code for initialising flow variables
!>
!> not changing storage position, exclude b.c. values, for example, developing
!> flow.
!> MPI : x-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!----------------------------------------------------------------------------------------------------------
!> \param[inout]  none          NA
!==========================================================================================================
  subroutine Calculate_xz_mean_yprofile(var, dtmp, n, varxz_work1)
    use mpi_mod
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), dimension(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3)), intent(in)  :: var ! x-pencil default
    integer,  intent(in)  :: n
    real(WP), dimension(n), optional, intent(out) :: varxz_work1

    real(wp) :: varxz( n )
    integer :: jj, i, j, k
    integer :: nk, ni!, nk_work, ni_work
    real(WP) :: varxz_work(n)
    !----------------------------------------------------------------------------------------------------------
    !   Default X-pencil
    !----------------------------------------------------------------------------------------------------------
    varxz = ZERO
    varxz_work = ZERO
    do j = 1, dtmp%xsz(2)
      nk = 0
      ni = 0
      jj = j - 1 + dtmp%xst(2)
      do k = 1, dtmp%xsz(3)
        nk = nk + 1
        do i = 1, dtmp%xsz(1)
          ni = ni + 1
          varxz(jj) = varxz(jj) + var(i, j, k) !
        end do
      end do
      varxz(jj) = varxz(jj) / real(nk * ni, wp)
    end do
    

    call mpi_barrier(MPI_COMM_WORLD, ierror)
    !call mpi_allreduce(ni, ni_work, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
    !call mpi_allreduce(nk, nk_work, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
    call mpi_allreduce(varxz, varxz_work, n, MPI_REAL_WP, MPI_SUM, MPI_COMM_WORLD, ierror)
    varxz_work = varxz_work / real(p_col * p_col, wp)
    if(PRESENT(varxz_work1)) varxz_work1 = varxz_work

#ifdef DEBUG_STEPS
    if (nrank == 0) then
      open(121, file = 'check_calculate_xz_mean_yprofile.dat', position="append")
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        write(121, *) jj, varxz_work(jj)
      end do
    end if
#endif
    
    
    return
  end subroutine
!==========================================================================================================
!> \brief : 
!> MPI : x-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!----------------------------------------------------------------------------------------------------------
!> \param[inout]            
!==========================================================================================================
  subroutine Adjust_to_xzmean_zero(var, dtmp, n, varxz)
    use mpi_mod
    use udf_type_mod
    implicit none
    type(DECOMP_INFO),  intent(in) :: dtmp
    integer,            intent(in) :: n
    real(WP), dimension(n), intent(in) :: varxz
    real(WP), dimension(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3)), intent(inout) :: var
    
    integer :: jj, i, j, k

    do j = 1, dtmp%xsz(2)
      jj = j - 1 + dtmp%xst(2)
      do k = 1, dtmp%xsz(3)
        do i = 1, dtmp%xsz(1)
          var(:, j, :) = var(:, j, :) - varxz(jj)
        end do
      end do
    end do

#ifdef DEBUG_STEPS
    open(121, file = 'check_adjust_to_xzmean_zero.dat', position="append")
    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
        do i = 1, dtmp%xsz(1)
          write(121, *) k, j, i, var(i, j, k)
        end do
      end do
    end do
    close(121)
#endif

    return
  end subroutine
!==========================================================================================================
!> \brief : 
!> MPI : x-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!----------------------------------------------------------------------------------------------------------
!> \param[inout]         
!==========================================================================================================
  subroutine Check_cfl_diffusion(x2r, rre, dt)
    use parameters_constant_mod
    !use iso_fortran_env
    use mpi_mod
    use wtformat_mod
    use print_msg_mod
    implicit none
    real(WP), intent(in) :: x2r(3)
    real(WP), intent(in) :: rre
    real(WP), intent(in) :: dt
    real(WP) :: cfl_diff

    ! check, ours is two times of the one in xcompact3d.
    cfl_diff = sum(x2r) * TWO * dt * rre
    if(nrank == 0) then
      if(cfl_diff > ONE) call Print_warning_msg("Warning: Diffusion number is larger than 1.")
      write (*, wrtfmt1e) "  Diffusion number :", cfl_diff
    end if
    
    return
  end subroutine
!==========================================================================================================
!> \brief : to check CFL for convection terms
!> CFL = u^x/dx + v^y/dy + w^z/dz < limit
!> MPI : x-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!----------------------------------------------------------------------------------------------------------
!> \param[inout]         
!==========================================================================================================
  subroutine Check_cfl_convection(u, v, w, dm)
    use parameters_constant_mod
    use udf_type_mod
    use operations
    use decomp_2d
    use wtformat_mod
    use find_max_min_ave_mod
    implicit none

    type(t_domain),               intent(in) :: dm
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent(in) :: u
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent(in) :: v
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent(in) :: w

    real(WP) :: var_xpencil (dm%dccc%xsz(1), &
                             dm%dccc%xsz(2), &
                             dm%dccc%xsz(3))
    real(WP) :: var_ypencil (dm%dccc%ysz(1), &
                             dm%dccc%ysz(2), &
                             dm%dccc%ysz(3))
    real(WP) :: var_zpencil (dm%dccc%zsz(1), &
                             dm%dccc%zsz(2), &
                             dm%dccc%zsz(3))
    real(WP) :: accc_xpencil (dm%dccc%xsz(1), &
                             dm%dccc%xsz(2), &
                             dm%dccc%xsz(3))
    real(WP) :: accc_ypencil (dm%dccc%ysz(1), &
                             dm%dccc%ysz(2), &
                             dm%dccc%ysz(3))
    real(WP) :: accc_zpencil (dm%dccc%zsz(1), &
                             dm%dccc%zsz(2), &
                             dm%dccc%zsz(3))
    real(WP) ::   v_ypencil (dm%dcpc%ysz(1), &
                             dm%dcpc%ysz(2), &
                             dm%dcpc%ysz(3))
    real(WP) ::   w_ypencil (dm%dccp%ysz(1), &
                             dm%dccp%ysz(2), &
                             dm%dccp%ysz(3))
    real(WP) ::   w_zpencil (dm%dccp%zsz(1), &
                             dm%dccp%zsz(2), &
                             dm%dccp%zsz(3))
    !real(WP)   :: cfl_convection, cfl_convection_work
    real(wp) :: dummy
!----------------------------------------------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------------------------------------------
    var_xpencil = ZERO
    var_ypencil = ZERO
    var_zpencil = ZERO
    accc_xpencil = ZERO
    accc_ypencil = ZERO
    accc_zpencil = ZERO
!----------------------------------------------------------------------------------------------------------
! X-pencil : u_ccc / dx * dt
!----------------------------------------------------------------------------------------------------------
    call Get_x_midp_P2C_3D(u, accc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_qx)
    var_xpencil = accc_xpencil * dm%h1r(1) * dm%dt
!----------------------------------------------------------------------------------------------------------
! Y-pencil : v_ccc / dy * dt
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(var_xpencil, var_ypencil, dm%dccc)
    call transpose_x_to_y(v,             v_ypencil, dm%dcpc)
    call Get_y_midp_P2C_3D(v_ypencil, accc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)
    var_ypencil = var_ypencil +  accc_ypencil * dm%h1r(2) * dm%dt
!----------------------------------------------------------------------------------------------------------
! Z-pencil : \overline{w}^z/dz at cell centre
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_z(var_ypencil, var_zpencil, dm%dccc)
    call transpose_x_to_y(w,             w_ypencil, dm%dccp)
    call transpose_y_to_z(w_ypencil,     w_zpencil, dm%dccp)
    call Get_z_midp_P2C_3D(w_zpencil, accc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)
    var_zpencil = var_zpencil +  accc_zpencil * dm%h1r(3) * dm%dt
!----------------------------------------------------------------------------------------------------------
! Z-pencil : Find the maximum 
!----------------------------------------------------------------------------------------------------------
    call Find_maximum_absvar3d(var_zpencil, dummy, dm%dccc, "CFL (convection) :", wrtfmt1e)

    ! if(nrank == 0) then
    !   if(cfl_convection_work > ONE) call Print_warning_msg("Warning: CFL is larger than 1.")
    !   write (*, wrtfmt1r) "  CFL (convection) :", cfl_convection_work
    ! end if
    
    return
  end subroutine
!==========================================================================================================
!>\brief : to calculate:
!>         fo = \int_1^nx \int_
!> This is based only y-direction stretching.
!> \todo Here is 2nd order Trapezoid Method. Need to improve! Check!
!---------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    needed         specified   pubic
!----------------------------------------------------------------------------------------------------------
!> MPI : 
!>     default x-pencil
!>     working in : y-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!> Y: index arrangment
!>      j'-1   j'-1  j'    j'+1  j'+2
!>      _|__.__|__.__|__.__|__.__|__.__
!>         j-2   j-1   j     j+1    j+2
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!----------------------------------------------------------------------------------------------------------
!> \param[inout]         
!==========================================================================================================
  subroutine Get_volumetric_average_3d(is_ynp, ibcy, fbcy, dm, dtmp, var, fo_work)
    use mpi_mod
    use udf_type_mod
    use parameters_constant_mod
    use operations
    use decomp_2d
    use wtformat_mod
    implicit none
    type(t_domain),  intent(in) :: dm
    logical,           intent(in) :: is_ynp
    integer,           intent(in) :: ibcy(2)
    real(WP),          intent(in) :: fbcy(:, :, :)
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP),          intent(in) :: var(:, :, :)
    real(WP),          intent(out):: fo_work
 
    real(WP), dimension( dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3) )  :: var_ypencil
    real(WP), allocatable   :: vcp_ypencil(:, :, :)
    real(WP)   :: vol, fo, vol_work
    integer :: i, j, k, noy, jp

! #ifdef DEBUG_STEPS  
!     if(nrank == 0) then
!       if(present(str)) then
!         call Print_debug_mid_msg("Calculating volumeric average of "//trim(str)//" in 3-D ...")
!       else
!         call Print_debug_mid_msg("Calculating volumeric average in 3-D ...")
!       end if
!     end if
! #endif

    if(.not. dm%is_stretching(2) ) then 
      vol = ZERO
      fo  = ZERO
      do k = 1, dtmp%xsz(3)
        do j = 1, dtmp%xsz(2)
          do i = 1, dtmp%xsz(1)
            fo = fo + var(i, j, k)
            vol = vol + ONE
          end do
        end do
      end do
      
    else
    !----------------------------------------------------------------------------------------------------------
    !   transpose to y pencil. Default is x-pencil.
    !----------------------------------------------------------------------------------------------------------
      var_ypencil = ZERO

      call transpose_x_to_y(var, var_ypencil, dtmp)
      !----------------------------------------------------------------------------------------------------------
      !   In Y-pencil now
      !----------------------------------------------------------------------------------------------------------
      if( is_ynp )then
        !----------------------------------------------------------------------------------------------------------
        !   if variable is stored in y-nodes, extend them to y-cell centres (P2C)
        !   for example, uy.
        !----------------------------------------------------------------------------------------------------------
        if( dm%is_periodic(2) ) then
          noy = dtmp%ysz(2)
        else
          noy = dtmp%ysz(2) - 1
        end if

        allocate( vcp_ypencil(dtmp%ysz(1), noy, dtmp%ysz(3)) )
        vcp_ypencil = ZERO

        call Get_y_midp_P2C_3D(var_ypencil, vcp_ypencil, dm, dm%iAccuracy, ibcy)

        fo = ZERO
        vol = ZERO
        do k = 1, dtmp%ysz(3)
          do i = 1, dtmp%ysz(1)
            do j = 1, noy
              !----------------------------------------------------------------------------------------------------------
              !       j'    j'+1
              !      _|__.__|_
              !         j     
              !----------------------------------------------------------------------------------------------------------
              jp = j + 1
              if( dm%is_periodic(2) .and. jp > dtmp%ysz(2)) jp = 1
              fo = fo + &      
                  ( var_ypencil(i, jp, k) + vcp_ypencil(i, j, k) ) * &
                  ( dm%yp(j + 1) - dm%yc(j) ) * HALF + &
                  ( var_ypencil(i, j,     k) + vcp_ypencil(i, j, k) ) * &
                  ( dm%yc(j    ) - dm%yp(j) ) * HALF
              vol = vol + ( dm%yp(j + 1) - dm%yp(j) )
            end do
          end do
        end do
        deallocate(vcp_ypencil)
      else
        !----------------------------------------------------------------------------------------------------------
        !   if variable is not stored in y-nodes, extends them to y-nodes. C2P
        !   for example, ux, density, etc.
        !----------------------------------------------------------------------------------------------------------
        if( dm%is_periodic(2) ) then
          noy = dtmp%ysz(2)
        else
          noy = dtmp%ysz(2) + 1
        end if
        allocate( vcp_ypencil(dtmp%ysz(1), noy, dtmp%ysz(3)) )
        vcp_ypencil = ZERO
        call Get_y_midp_C2P_3D(var_ypencil, vcp_ypencil, dm, dm%iAccuracy, ibcy, fbcy)

        fo = ZERO
        vol = ZERO
        do k = 1, dtmp%ysz(3)
          do i = 1, dtmp%ysz(1)
            do j = 1, dtmp%ysz(2)
              !----------------------------------------------------------------------------------------------------------
              !      j'    j'+1
              !      _|__.__|_
              !         j
              !----------------------------------------------------------------------------------------------------------
              jp = j + 1
              if( dm%is_periodic(2) .and. jp > noy) jp = 1
              ! method 1: 2nd order
              ! fo = fo + &
              !     ( vcp_ypencil(i, jp, k) + var_ypencil(i, j, k) ) * &
              !     ( dm%yp(j + 1) - dm%yc(j) ) * HALF + &
              !     ( var_ypencil(i, j,     k) + var_ypencil(i, j, k) ) * &
              !     ( dm%yc(j    ) - dm%yp(j) ) * HALF
              ! method 2: 1st order, same as CHAPSim1
              fo = fo + vcp_ypencil(i, j, k)*(dm%yp(j + 1) - dm%yp(j))
              vol = vol + ( dm%yp(j + 1) - dm%yp(j) )
            end do
          end do
        end do
        deallocate(vcp_ypencil)
      end if

    end if


    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce( fo,  fo_work, 1, MPI_REAL_WP, MPI_SUM, MPI_COMM_WORLD, ierror)
    call mpi_allreduce(vol, vol_work, 1, MPI_REAL_WP, MPI_SUM, MPI_COMM_WORLD, ierror)
    fo_work = fo_work / vol_work

#ifdef DEBUG_STEPS  
    if(nrank == 0 ) then
      write (*, wrtfmt1e) " volumetric average is ", fo_work
    end if
#endif

    return 
  end subroutine Get_volumetric_average_3d

end module
