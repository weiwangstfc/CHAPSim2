module solver_tools_mod

  ! procedure
  private

  public  :: Calculate_massflux_from_velocity
  public  :: Update_thermal_properties
  public  :: Check_cfl_convection
  public  :: Check_cfl_diffusion

  public  :: Update_Re
  public  :: Update_PrGr
  public  :: Calculate_xz_mean_yprofile
  public  :: Adjust_to_xzmean_zero
  public  :: Get_volumetric_average_3d
  public  :: Find_maximum_absvar3d

contains
!=============================================================================================================================================
!> \brief The main code for initializing flow variables
!> This subroutine is called once in \ref Initialize_chapsim.
!>
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[inout]  
!=============================================================================================================================================
  subroutine Update_Re(iter, fl)
    use parameters_constant_mod
    use thermo_info_mod
    use udf_type_mod
    implicit none
    integer,     intent(in ) :: iter  
    type(t_flow),   intent(inout) :: fl

  !---------------------------------------------------------------------------------------------------------------------------------------------
  !  1/Re                                   
  !---------------------------------------------------------------------------------------------------------------------------------------------
    if(iter < fl%nIterIniRen) then
      fl%rre = ONE / fl%renIni
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
  
!---------------------------------------------------------------------------------------------------------------------------------------------
!  1/(Re*Pr)                                   
!---------------------------------------------------------------------------------------------------------------------------------------------
    tm%rPrRen = fl%rre * fluidparam%ftp0ref%k / fluidparam%ftp0ref%m / fluidparam%ftp0ref%cp
!---------------------------------------------------------------------------------------------------------------------------------------------
!  gravity force                          
!---------------------------------------------------------------------------------------------------------------------------------------------  
    u0 = ONE / fl%rre * fluidparam%ftp0ref%m / fluidparam%ftp0ref%d / tm%lenRef
    rtmp = tm%lenRef / u0 / u0 * GRAVITY
    fl%fgravity(:) = ZERO
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
      fl%fgravity(:) = ZERO
    end if
    
    return
  end subroutine Update_PrGr
!=============================================================================================================================================
!> \brief The main code for initializing flow variables
!>
!> not changing storage position, exclude b.c. values, for example, developing
!> flow.
!> MPI : x-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[inout]  none          NA
!=============================================================================================================================================
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
    integer :: nk, ni, nk_work, ni_work
    real(WP) :: varxz_work(n)
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   Default X-pencil
    !---------------------------------------------------------------------------------------------------------------------------------------------
    varxz(:) = ZERO
    varxz_work(:) = ZERO
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
    call mpi_allreduce(varxz, varxz_work, n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    varxz_work(:) = varxz_work(:) / real(p_col, wp)
    if(PRESENT(varxz_work1)) varxz_work1 = varxz_work

#ifdef DEBUG
    if (nrank == 0) then
      call Print_debug_mid_msg("xz plane averaged data:")
      do j = 1, n
        write(*, *) n, varxz_work(j)
      end do
    end if
#endif
    
    
    return
  end subroutine
!=============================================================================================================================================
!> \brief : 
!> MPI : x-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[inout]            
!=============================================================================================================================================
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

    return
  end subroutine
!=============================================================================================================================================
!> \brief Calculate the conservative variables from primary variable.     
!> This subroutine is called to update $\rho u_i$ from $u_i$.
!--------------------------------------------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    needed         specified   pubic
!---------------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     dm             domain
!> \param[in]     fm             flow
!=============================================================================================================================================
  subroutine Calculate_massflux_from_velocity(fl, dm)
    use udf_type_mod
    use operations
    implicit none
    type(t_domain), intent(in )   :: dm
    type(t_flow  ), intent(inout) :: fl
    real(WP), dimension( dm%nc(1) ) :: fix
    real(WP), dimension( dm%np(1) ) :: fox
    real(WP), dimension( dm%nc(2) ) :: fiy
    real(WP), dimension( dm%np(2) ) :: foy
    real(WP), dimension( dm%nc(3) ) :: fiz
    real(WP), dimension( dm%np(3) ) :: foz
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) ::  uy_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: duy_ypencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) ::  uz_ypencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) :: duz_ypencil
    real(WP), dimension( dm%dccc%ysz(1), &
                         dm%dccc%ysz(2), &
                         dm%dccc%ysz(3)) ::   d_ypencil
    
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) ::  uz_zpencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: duz_zpencil
    real(WP), dimension( dm%dccc%zsz(1), &
                         dm%dccc%zsz(2), &
                         dm%dccc%zsz(3)) ::   d_zpencil
    
    integer :: i, j, k
    type(DECOMP_INFO) :: dtmp

    integer  :: ibc(2)
    real(WP) :: fbc(2)

!---------------------------------------------------------------------------------------------------------------------------------------------
! x-pencil : u1 -> g1
!---------------------------------------------------------------------------------------------------------------------------------------------
    dtmp = dm%dpcc
    ibc(:) = dm%ibcx(:, 5)
    fbc(:) = dm%fbc_dend(:, 1)
    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
        fix(:) = fl%dDens(:, j, k)
        call Get_x_midp_C2P_1D (fix, fox, dm, ibc, fbc)
        fl%gx(:, j, k) = fox(:) * fl%qx(:, j, k)
      end do
    end do
!---------------------------------------------------------------------------------------------------------------------------------------------
! x-pencil --> y-pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qy,    uy_ypencil, dm%dcpc)
    call transpose_x_to_y(fl%dDens,  d_ypencil, dm%dccc)
    call transpose_x_to_y(fl%qz,    uz_ypencil, dm%dccp)
!---------------------------------------------------------------------------------------------------------------------------------------------
! y-pencil : u2 -> g2
!---------------------------------------------------------------------------------------------------------------------------------------------
    dtmp = dm%dcpc
    ibc(:) = dm%ibcy(:, 5)
    fbc(:) = dm%fbc_dend(:, 2)
    do k = 1, dtmp%ysz(3)
      do i = 1, dtmp%ysz(1)
        fiy(:) = d_ypencil(i, :, k)
        call Get_y_midp_C2P_1D (fiy, foy, dm, ibc, fbc)
        duy_ypencil(i, :, k) = foy(:) * uy_ypencil(i, :, k)
      end do
    end do
!---------------------------------------------------------------------------------------------------------------------------------------------
! y-pencil --> z-pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
    call transpose_y_to_z( d_ypencil,  d_zpencil, dm%dccc)
    call transpose_y_to_z(uz_ypencil, uz_zpencil, dm%dccp)
!---------------------------------------------------------------------------------------------------------------------------------------------
! Z-pencil : u3 -> g3
!---------------------------------------------------------------------------------------------------------------------------------------------
    dtmp = dm%dcpc
    ibc(:) = dm%ibcz(:, 5)
    fbc(:) = dm%fbc_dend(:, 3)
    do j = 1, dtmp%zsz(2)
      do i = 1, dtmp%zsz(1)
        fiz(:) = d_zpencil(i, j, :)
        call Get_z_midp_C2P_1D (fiz, foz, dm, ibc, fbc)
        duz_zpencil(i, j, :) = foz(:) * uz_zpencil(i, j, :)
      end do
    end do
!---------------------------------------------------------------------------------------------------------------------------------------------
! z-pencil --> y-pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
    call transpose_z_to_y(duz_zpencil, duz_ypencil, dm%dccp)
!---------------------------------------------------------------------------------------------------------------------------------------------
! y-pencil --> x-pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
    call transpose_y_to_x(duz_ypencil, fl%gz, dm%dccp)
    call transpose_y_to_x(duy_ypencil, fl%gy, dm%dcpc)

    return
  end subroutine Calculate_massflux_from_velocity
!=============================================================================================================================================
!> \brief : 
!> MPI : x-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[inout]         
!=============================================================================================================================================
  subroutine Check_cfl_diffusion(x2r, rre, dt)
    use parameters_constant_mod
    !use iso_fortran_env
    use mpi_mod
    use wtformat_mod
    implicit none
    real(WP), intent(in) :: x2r(3)
    real(WP), intent(in) :: rre
    real(WP), intent(in) :: dt
    real(WP) :: cfl_diff

    ! check, ours is two times of the one in xcompact3d.
    cfl_diff = sum(x2r) * TWO * dt * rre
    if(nrank == 0) then
      if(cfl_diff > ONE) call Print_warning_msg("Warning: Diffusion number is larger than 1.")
      write (*, wrtfmt1r) "  Diffusion number :", cfl_diff
    end if
    
    return
  end subroutine
!=============================================================================================================================================
!> \brief : to check CFL for convection terms
!> CFL = u^x/dx + v^y/dy + w^z/dz < limit
!> MPI : x-pencil
!>  (y) ^_____ _____ ______
!>      |_____|_____|______|
!>      |_____|_____|______|__> (z)
!>
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[inout]         
!=============================================================================================================================================
  subroutine Check_cfl_convection(u, v, w, dm)
    use parameters_constant_mod
    use udf_type_mod
    use operations
    use decomp_2d
    use wtformat_mod
    implicit none

    type(t_domain),               intent(in) :: dm
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent(in) :: u
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent(in) :: v
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent(in) :: w

    real(WP), dimension( dm%np(1) ) :: fix
    real(WP), dimension( dm%np(2) ) :: fiy
    real(WP), dimension( dm%np(3) ) :: fiz
    real(WP), dimension( dm%nc(1) ) :: fox
    real(WP), dimension( dm%nc(2) ) :: foy
    real(WP), dimension( dm%nc(3) ) :: foz

    real(WP) :: udx_xpencil (dm%dccc%xsz(1), &
                             dm%dccc%xsz(2), &
                             dm%dccc%xsz(3))
    real(WP) :: udx_ypencil (dm%dccc%ysz(1), &
                             dm%dccc%ysz(2), &
                             dm%dccc%ysz(3))
    real(WP) :: udx_zpencil (dm%dccc%zsz(1), &
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
    real(WP)   :: cfl_convection, cfl_convection_work
    integer :: i, j, k
    type(DECOMP_INFO) :: dtmp

!---------------------------------------------------------------------------------------------------------------------------------------------
! Initialisation
!---------------------------------------------------------------------------------------------------------------------------------------------
    udx_xpencil(:, :, :) = ZERO
    udx_ypencil(:, :, :) = ZERO
    udx_zpencil(:, :, :) = ZERO
      v_ypencil(:, :, :) = ZERO
      w_ypencil(:, :, :) = ZERO
      w_zpencil(:, :, :) = ZERO
!---------------------------------------------------------------------------------------------------------------------------------------------
! X-pencil : u_ccc / dx * dt
!---------------------------------------------------------------------------------------------------------------------------------------------
    dtmp = dm%dpcc
    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
        fix(:) = u(:, j, k)
        call Get_x_midp_P2C_1D (fix, fox, dm, dm%ibcx(:, 1))
        udx_xpencil(:, j, k) = fox(:) * dm%h1r(1) * dm%dt
      end do
    end do
!---------------------------------------------------------------------------------------------------------------------------------------------
! Convert X-pencil to Y-Pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(udx_xpencil, udx_ypencil, dm%dccc)
    call transpose_x_to_y(v,             v_ypencil, dm%dcpc)
    call transpose_x_to_y(w,             w_ypencil, dm%dccp)
!---------------------------------------------------------------------------------------------------------------------------------------------
! Y-pencil : v_ccc / dy * dt
!---------------------------------------------------------------------------------------------------------------------------------------------
    dtmp = dm%dcpc
    do k = 1, dtmp%ysz(3)
      do i = 1, dtmp%ysz(1)
        fiy(:) = v_ypencil(i, :, k)
        call Get_y_midp_P2C_1D (fiy, foy, dm, dm%ibcy(:, 2))
        udx_ypencil(i, :, k) = udx_ypencil(i, :, k) + foy(:) * dm%h1r(2) * dm%dt
      end do
    end do
!---------------------------------------------------------------------------------------------------------------------------------------------
! Convert Y-pencil to Z-Pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
    call transpose_y_to_z(udx_ypencil, udx_zpencil, dm%dccc)
    call transpose_y_to_z(  w_ypencil,   w_zpencil, dm%dccp)
!---------------------------------------------------------------------------------------------------------------------------------------------
! Z-pencil : \overline{w}^z/dz at cell centre
!---------------------------------------------------------------------------------------------------------------------------------------------
    dtmp = dm%dccp
    do j = 1, dtmp%zsz(2)
      do i = 1, dtmp%zsz(1)
        fiz(:) = w_zpencil(i, j, :)
        call Get_z_midp_P2C_1D (fiz, foz, dm, dm%ibcz(:, 3))
        udx_zpencil(i, j, :) = udx_zpencil(i, j, :) + foz(:) * dm%h1r(3) * dm%dt
      end do
    end do
!---------------------------------------------------------------------------------------------------------------------------------------------
! Z-pencil : Find the maximum 
!---------------------------------------------------------------------------------------------------------------------------------------------
    call Find_maximum_absvar3d(udx_zpencil, "CFL (convection) :")

    ! if(nrank == 0) then
    !   if(cfl_convection_work > ONE) call Print_warning_msg("Warning: CFL is larger than 1.")
    !   write (*, wrtfmt1r) "  CFL (convection) :", cfl_convection_work
    ! end if
    
    return
  end subroutine
!=============================================================================================================================================
!>\brief : to calculate:
!>         fo = \int_1^nx \int_
!> This is based only y-direction stretching.
!> \todo Here is 2nd order Trapezoid Method. Need to improve! Check!
!--------------------------------------------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    needed         specified   pubic
!---------------------------------------------------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[inout]         
!=============================================================================================================================================
  subroutine Get_volumetric_average_3d(is_ynp, ibcy, fbcy, dm, dtmp, var, fo_work)
    use mpi_mod
    use udf_type_mod
    use parameters_constant_mod
    use operations
    implicit none
    type(t_domain),  intent(in) :: dm
    logical,           intent(in) :: is_ynp
    integer,           intent(in) :: ibcy(2)
    real(WP),          intent(in) :: fbcy(2)
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP),          intent(in) :: var(:, :, :)
    real(WP),          intent(out):: fo_work
 
    real(WP), dimension( dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3) )  :: var_ypencil
    real(WP), allocatable   :: vcp_ypencil(:, :, :)
    real(WP)   :: vol, fo, vol_work
    integer :: i, j, k, noy, jp

  !---------------------------------------------------------------------------------------------------------------------------------------------
  !   transpose to y pencil. Default is x-pencil.
  !---------------------------------------------------------------------------------------------------------------------------------------------
    var_ypencil = ZERO

    call transpose_x_to_y(var, var_ypencil, dtmp)
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   In Y-pencil now
    !---------------------------------------------------------------------------------------------------------------------------------------------
    if( is_ynp )then
      !---------------------------------------------------------------------------------------------------------------------------------------------
      !   if variable is stored in y-nodes, extend them to y-cell centres (P2C)
      !   for example, uy.
      !---------------------------------------------------------------------------------------------------------------------------------------------
      if( dm%is_periodic(2) ) then
        noy = dtmp%ysz(2)
      else
        noy = dtmp%ysz(2) - 1
      end if

      allocate( vcp_ypencil(dtmp%ysz(1), noy, dtmp%ysz(3)) )
      vcp_ypencil = ZERO

      call Get_y_midp_P2C_3D(var_ypencil, vcp_ypencil, dm, ibcy)

      fo = ZERO
      vol = ZERO
      do k = 1, dtmp%ysz(3)
        do i = 1, dtmp%ysz(1)
          do j = 1, noy
            !---------------------------------------------------------------------------------------------------------------------------------------------
            !       j'    j'+1
            !      _|__.__|_
            !         j     
            !---------------------------------------------------------------------------------------------------------------------------------------------
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
      !---------------------------------------------------------------------------------------------------------------------------------------------
      !   if variable is not stored in y-nodes, extends them to y-nodes. C2P
      !   for example, ux, density, etc.
      !---------------------------------------------------------------------------------------------------------------------------------------------
      if( dm%is_periodic(2) ) then
        noy = dtmp%ysz(2)
      else
        noy = dtmp%ysz(2) + 1
      end if
      allocate( vcp_ypencil(dtmp%ysz(1), noy, dtmp%ysz(3)) )
      vcp_ypencil = ZERO
      call Get_y_midp_C2P_3D(var_ypencil, vcp_ypencil, dm, ibcy, fbcy)

      fo = ZERO
      vol = ZERO
      do k = 1, dtmp%ysz(3)
        do i = 1, dtmp%ysz(1)
          do j = 1, dtmp%ysz(2)
            !---------------------------------------------------------------------------------------------------------------------------------------------
            !      j'    j'+1
            !      _|__.__|_
            !         j
            !---------------------------------------------------------------------------------------------------------------------------------------------
            jp = j + 1
            if( dm%is_periodic(2) .and. jp > noy) jp = 1
            fo = fo + &
                ( vcp_ypencil(i, jp, k) + var_ypencil(i, j, k) ) * &
                ( dm%yp(j + 1) - dm%yc(j) ) * HALF + &
                ( var_ypencil(i, j,     k) + var_ypencil(i, j, k) ) * &
                ( dm%yc(j    ) - dm%yp(j) ) * HALF
            vol = vol + ( dm%yp(j + 1) - dm%yp(j) )
          end do
        end do
      end do
      deallocate(vcp_ypencil)
    end if
    
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce( fo,  fo_work, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    call mpi_allreduce(vol, vol_work, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    fo_work = fo_work / vol_work

    return 
  end subroutine Get_volumetric_average_3d

!=============================================================================================================================================
  subroutine Find_maximum_absvar3d(var,  str)
    use precision_mod
    use math_mod
    use mpi_mod
    use wtformat_mod
    implicit none

    real(WP), intent(in)  :: var(:, :, :)
    character(len = *), intent(in) :: str
    
    real(WP):: varmax_work
    real(WP)   :: varmax

    integer :: i, j, k, nx, ny, nz
    nx = size(var, 1)
    ny = size(var, 2)
    nz = size(var, 3)

    varmax = ZERO
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if(abs_wp(var(i, j, k)) > varmax) varmax = abs_wp(var(i, j, k))
        end do
      end do
    end do

    !varmax = MAXVAL( abs_wp( var(:, :, :) ) ) 
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(varmax, varmax_work, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)

    if(nrank == 0) then
      write (*, wrtfmt1e) str, varmax_work
    end if

    return
  end subroutine
!=============================================================================================================================================
  subroutine Update_thermal_properties(fl, tm, dm)
    use udf_type_mod
    use thermo_info_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm

    integer :: i, j, k
    type(t_fluidThermoProperty) :: ftp
!---------------------------------------------------------------------------------------------------------------------------------------------
!   x-pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
    do k = dm%dccc%xst(3), dm%dccc%xen(3)
      do j = dm%dccc%xst(2), dm%dccc%xen(2)
        do i = dm%dccc%xst(1), dm%dccc%xen(1)
          ftp%dh = tm%dh(i, j, k)
          call ftp_refresh_thermal_properties_from_DH(ftp)
          tm%hEnth(i, j, k) = ftp%h
          tm%tTemp(i, j, k) = ftp%T
          tm%kCond(i, j, k) = ftp%k
          fl%dDens(i, j, k) = ftp%d
          fl%mVisc(i, j, k) = ftp%m
        end do
      end do
    end do

    fl%dDensm1(:, :, :) = fl%dDens(:, :, :)
    fl%dDensm2(:, :, :) = fl%dDensm1(:, :, :)

  return
  end subroutine Update_thermal_properties

end module
