!----------------------------------------------------------------------------------------------------------
!                      CHAPSim version 2.0.0
!                      --------------------------
! This file is part of CHAPSim, a general-purpose CFD tool.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 3 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.
!----------------------------------------------------------------------------------------------------------
module flow_thermo_initialiasation
  use vars_df_mod
  use solver_tools_mod
  implicit none

  private :: Allocate_flow_variables
  private :: Allocate_thermo_variables
  private :: Generate_poiseuille_flow_profile
  private :: Generate_random_field

  private :: Initialize_poiseuille_flow
  private :: Initialize_flow_from_given_values
  private :: Initialize_vortexgreen_2dflow
  private :: Initialize_vortexgreen_3dflow

  public  :: Validate_TGV2D_error
  public  :: Initialize_flow_fields
  public  :: Initialize_thermo_fields

contains
!==========================================================================================================
!> \brief Allocate flow and thermal variables.     
!---------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   private
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                         
!----------------------------------------------------------------------------------------------------------
!> \param[in]     none          NA
!> \param[out]    none          NA
!==========================================================================================================
  subroutine Allocate_flow_variables (dm, fl)
    use parameters_constant_mod
    use mpi_mod
    implicit none

    type(t_domain), intent(in)    :: dm
    type(t_flow),   intent(inout) :: fl

    if(nrank == 0) call Print_debug_start_msg("Allocating flow variables ...")
    !----------------------------------------------------------------------------------------------------------
    ! default : x pencil. 
    ! varaible index is LOCAL. means 1:xsize(1)
    !----------------------------------------------------------------------------------------------------------
    call alloc_x(fl%qx,      dm%dpcc) ; fl%qx = ZERO ! qx = ux
    call alloc_x(fl%qy,      dm%dcpc) ; fl%qy = ZERO ! qy = ur * rp
    call alloc_x(fl%qz,      dm%dccp) ; fl%qz = ZERO ! qz = u_theta * rc

    call alloc_x(fl%pres,    dm%dccc) ; fl%pres = ZERO
    call alloc_x(fl%pcor,    dm%dccc) ; fl%pcor = ZERO
    call alloc_z(fl%pcor_zpencil_ggg,    dm%dccc, .true.) ; fl%pcor_zpencil_ggg = ZERO

    call alloc_x(fl%mx_rhs,  dm%dpcc) ; fl%mx_rhs = ZERO
    call alloc_x(fl%my_rhs,  dm%dcpc) ; fl%my_rhs = ZERO
    call alloc_x(fl%mz_rhs,  dm%dccp) ; fl%mz_rhs = ZERO

    call alloc_x(fl%mx_rhs0, dm%dpcc) ; fl%mx_rhs0 = ZERO
    call alloc_x(fl%my_rhs0, dm%dcpc) ; fl%my_rhs0 = ZERO
    call alloc_x(fl%mz_rhs0, dm%dccp) ; fl%mz_rhs0 = ZERO

    allocate (fl%fbcx_qx(4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)))
    allocate (fl%fbcx_qy(4, dm%dcpc%xsz(2), dm%dcpc%xsz(3)))
    allocate (fl%fbcx_qz(4, dm%dccp%xsz(2), dm%dccp%xsz(3)))

    allocate (fl%fbcy_qx(dm%dpcc%ysz(1), 4, dm%dpcc%ysz(3)))
    allocate (fl%fbcy_qy(dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3)))
    allocate (fl%fbcy_qz(dm%dccp%ysz(1), 4, dm%dccp%ysz(3)))

    allocate (fl%fbcz_qx(dm%dpcc%zsz(1), dm%dpcc%zsz(2), 4))
    allocate (fl%fbcz_qy(dm%dcpc%zsz(1), dm%dcpc%zsz(2), 4))
    allocate (fl%fbcz_qz(dm%dccp%zsz(1), dm%dccp%zsz(2), 4))

    allocate (fl%fbcx_pr(4, dm%dccc%xsz(2), dm%dccc%xsz(3)))
    allocate (fl%fbcy_pr(dm%dccc%ysz(1), 4, dm%dccc%ysz(3)))
    allocate (fl%fbcz_pr(dm%dccc%zsz(1), dm%dccc%zsz(2), 4))


    if(dm%icoordinate == ICYLINDRICAL) then
      allocate (fl%fbcy_qyr(dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3)))
      allocate (fl%fbcz_qyr(dm%dcpc%zsz(1), dm%dcpc%zsz(2), 4))
      allocate (fl%fbcy_qzr(dm%dccp%ysz(1), 4, dm%dccp%ysz(3)))
      allocate (fl%fbcz_qzr(dm%dccp%zsz(1), dm%dccp%zsz(2), 4))
    end if

    if(dm%is_thermo) then
      call alloc_x(fl%gx,      dm%dpcc) ; fl%gx = ZERO ! gx = rho * qx = rho * ux
      call alloc_x(fl%gy,      dm%dcpc) ; fl%gy = ZERO ! gy = rho * qy = rho * ur * rp
      call alloc_x(fl%gz,      dm%dccp) ; fl%gz = ZERO ! gz = rho * qz = rho * u_theta * rc

      allocate (fl%fbcx_gx(4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)))
      allocate (fl%fbcx_gy(4, dm%dcpc%xsz(2), dm%dcpc%xsz(3)))
      allocate (fl%fbcx_gz(4, dm%dccp%xsz(2), dm%dccp%xsz(3)))
  
      allocate (fl%fbcy_gx(dm%dpcc%ysz(1), 4, dm%dpcc%ysz(3)))
      allocate (fl%fbcy_gy(dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3)))
      allocate (fl%fbcy_gz(dm%dccp%ysz(1), 4, dm%dccp%ysz(3)))
  
      allocate (fl%fbcz_gx(dm%dpcc%zsz(1), dm%dpcc%zsz(2), 4))
      allocate (fl%fbcz_gy(dm%dcpc%zsz(1), dm%dcpc%zsz(2), 4))
      allocate (fl%fbcz_gz(dm%dccp%zsz(1), dm%dccp%zsz(2), 4))
  
      if(dm%icoordinate == ICYLINDRICAL) then
        allocate (fl%fbcy_gyr(dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3)))
        allocate (fl%fbcz_gyr(dm%dcpc%zsz(1), dm%dcpc%zsz(2), 4))
        allocate (fl%fbcy_gzr(dm%dccp%ysz(1), 4, dm%dccp%ysz(3)))
        allocate (fl%fbcz_gzr(dm%dccp%zsz(1), dm%dccp%zsz(2), 4))
      end if

    end if

    if(nrank == 0) call Print_debug_end_msg
    return

  end subroutine Allocate_flow_variables
  !==========================================================================================================
  subroutine Allocate_thermo_variables (dm, tm)
    use parameters_constant_mod
    use mpi_mod
    use udf_type_mod
    use thermo_info_mod
    implicit none

    type(t_domain), intent(in)    :: dm
    type(t_thermo), intent(inout) :: tm

    if(.not. dm%is_thermo) return
    if(nrank == 0) call Print_debug_start_msg("Allocating thermal variables ...")
    !----------------------------------------------------------------------------------------------------------
    ! default : x pencil. 
    ! varaible index is LOCAL. means 1:xsize(1)
    !----------------------------------------------------------------------------------------------------------
    call alloc_x(tm%dh,    dm%dccc) ; tm%dh    = ZERO
    call alloc_x(tm%hEnth, dm%dccc) ; tm%hEnth = ZERO
    call alloc_x(tm%kCond, dm%dccc) ; tm%kCond = ONE
    call alloc_x(tm%tTemp, dm%dccc) ; tm%tTemp = ONE

    call alloc_x(tm%dDens,   dm%dccc) ; tm%dDens = ONE
    call alloc_x(tm%mVisc,   dm%dccc) ; tm%mVisc = ONE

    call alloc_x(tm%dDensm1, dm%dccc) ; tm%dDensm1 = ONE
    call alloc_x(tm%dDensm2, dm%dccc) ; tm%dDensm2 = ONE

    if(nrank == 0) call Print_debug_end_msg
    return

  end subroutine Allocate_thermo_variables
  !==========================================================================================================
  !> \brief Generate a flow profile for Poiseuille flow in channel or pipe.     
  !---------------------------------------------------------------------------------------------------------- 
  !> Scope:  mpi    called-freq    xdomain     module
  !>         all    once           specified   private
  !----------------------------------------------------------------------------------------------------------
  ! Arguments
  !----------------------------------------------------------------------------------------------------------
  !  mode           name          role                                           
  !----------------------------------------------------------------------------------------------------------
  !> \param[in]     
  !> \param[out]    
  !==========================================================================================================
  subroutine Generate_random_field(dm, fl)
    use random_number_generation_mod
    use parameters_constant_mod
    use mpi_mod
    use math_mod
    use boundary_conditions_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow), intent(inout) :: fl
    
    integer :: seed
    integer :: i, j, k! local id
    integer :: n, nsz  
    integer :: ii, jj, kk ! global id
    real(WP) :: rd, lownoise
    type(DECOMP_INFO) :: dtmp

    if(nrank == 0) call Print_debug_mid_msg("Generating random field ...")
    !----------------------------------------------------------------------------------------------------------
    !   Initialisation in x pencil
    !----------------------------------------------------------------------------------------------------------
    seed = 0
    fl%pres(:, :, :) = ZERO
    fl%qx(:, :, :) = ZERO
    fl%qy(:, :, :) = ZERO
    fl%qz(:, :, :) = ZERO
    nsz = dm%np(1) * dm%np(2) * dm%np(3)

    do n = 1, NDIM
      if(n == 1) then
        dtmp = dm%dpcc
      else if(n == 2) then
        dtmp = dm%dcpc
      else if(n == 3) then
        dtmp = dm%dccp
      else
      end if

      do i = 1, dtmp%xsz(1)
        ii = i
        do j = 1, dtmp%xsz(2)
          jj = dtmp%xst(2) + j - 1
          ! if( ( 1.0_WP - abs_wp(dm%yp(jj)) ) < 0.250_WP) then
          !   lownoise = lnoise * lnoise
          ! else
             lownoise = fl%noiselevel
          ! end if
          do k = 1, dtmp%xsz(3)
            kk = dtmp%xst(3) + k - 1
            seed = (nrank + 1) * nsz * n + ii * 313 + jj * 571 + kk * 937
            !seed = ii + jj + kk + nsz * (n - 1)
            !write(*,*) ii, jj, kk, seed
            call Initialize_random_number ( seed )
            call Generate_r_random( -ONE, ONE, rd)
            if(n == 1) fl%qx(i, j, k) = lownoise * rd
            if(n == 2) fl%qy(i, j, k) = lownoise * rd / dm%rpi(jj)
            if(n == 3) fl%qz(i, j, k) = lownoise * rd / dm%rci(jj)

          end do
        end do
      end do
    end do

    !----------------------------------------------------------------------------------------------------------
    ! Correct dirichlet bc for random 
    !----------------------------------------------------------------------------------------------------------
    if(dm%ibcx(1, 1) == IBC_DIRICHLET) fl%qx(1,              :, :) = ZERO
    if(dm%ibcx(2, 1) == IBC_DIRICHLET) fl%qx(dm%dpcc%xsz(1), :, :) = ZERO

    if(dm%ibcy(1, 2) == IBC_DIRICHLET .and. dm%dcpc%xst(2) == 1)        fl%qy(:, 1,              :) = ZERO
    if(dm%ibcy(2, 2) == IBC_DIRICHLET .and. dm%dcpc%xen(2) == dm%np(2)) fl%qy(:, dm%dcpc%xsz(2), :) = ZERO

    if(dm%ibcz(1, 3) == IBC_DIRICHLET .and. dm%dccp%xst(3) == 1)        fl%qz(:, :, 1             ) = ZERO
    if(dm%ibcz(2, 3) == IBC_DIRICHLET .and. dm%dccp%xen(3) == dm%np(3)) fl%qz(:, :, dm%dccp%xsz(3)) = ZERO
    

    return
  end subroutine

  !==========================================================================================================
  !> \brief Generate a flow profile for Poiseuille flow in channel or pipe.     
  !---------------------------------------------------------------------------------------------------------- 
  !> Scope:  mpi    called-freq    xdomain     module
  !>         all    once           specified   private
  !----------------------------------------------------------------------------------------------------------
  ! Arguments
  !----------------------------------------------------------------------------------------------------------
  !  mode           name          role                                           
  !----------------------------------------------------------------------------------------------------------
  !> \param[in]     d             domain
  !> \param[out]    ux_1c1          u(yc), velocity profile along wall-normal direction
  !==========================================================================================================
  subroutine Generate_poiseuille_flow_profile(dm, ux_1c1)
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none

    type(t_domain), intent(in)  :: dm
    real(WP),       intent(out) :: ux_1c1(:)
    
    real(WP)   :: a, b, c, yy, ymax, ymin
    integer :: j
    
    if(nrank == 0) call Print_debug_mid_msg("Generate poiseuille flow profile ...")

    ux_1c1 (:) = ZERO

    ymax = dm%yp( dm%np_geo(2) )
    ymin = dm%yp( 1 )
    if (dm%icase == ICASE_CHANNEL) then
      a = (ymax - ymin) * HALF
      b = ZERO
      c = ONEPFIVE
    else if (dm%icase == ICASE_PIPE) then
      a = (ymax - ymin)
      b = ZERO
      c = TWO
    else if (dm%icase == ICASE_ANNUAL) then
      a = (ymax - ymin) * HALF
      b = (ymax + ymin) * HALF
      c = TWO
    else 
      a = (ymax - ymin) * HALF
      b = ZERO
      c = ONEPFIVE
    end if

    do j = 1, dm%nc(2)
      yy = dm%yc(j)
      ux_1c1(j) = ( ONE - ( (yy - b)**2 ) / a / a ) * c
    end do

    return
  end subroutine Generate_poiseuille_flow_profile

  !==========================================================================================================
  !> \brief Initialize Poiseuille flow in channel or pipe.     
  !---------------------------------------------------------------------------------------------------------- 
  !> Scope:  mpi    called-freq    xdomain     module
  !>         all    once           specified   private
  !----------------------------------------------------------------------------------------------------------
  ! Arguments
  !----------------------------------------------------------------------------------------------------------
  !  mode           name          role                                           
  !----------------------------------------------------------------------------------------------------------
  !> \param[in]     d             domain
  !> \param[out]    f             flow
  !==========================================================================================================
  subroutine Initialize_poiseuille_flow(dm, fl)
    use input_general_mod
    use udf_type_mod
    use boundary_conditions_mod
    use parameters_constant_mod
    use wtformat_mod
    use files_io_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow), intent(inout) :: fl
    integer :: pf_unit
    integer :: i, j, k, jj
    real(WP) :: ubulk
    real(WP) :: ux_1c1(dm%nc(2))
    ! real(WP) :: uxxza (dm%nc(2))
    ! real(WP) :: uyxza (dm%np(2))
    ! real(WP) :: uzxza (dm%nc(2))
    real(WP) :: ux_ypencil(dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3))

    type(DECOMP_INFO) :: dtmp

    if(nrank == 0) call Print_debug_start_msg("Initializing Poiseuille flow field ...")
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : initial
    !----------------------------------------------------------------------------------------------------------
    fl%pres(:, :, :) = ZERO
    fl%qx(:, :, :) = ZERO
    fl%qy(:, :, :) = ZERO
    fl%qz(:, :, :) = ZERO
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : to get random fields [-1,1] for qx, qy, qz
    !----------------------------------------------------------------------------------------------------------
    call Generate_random_field(dm, fl)
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : to get Poiseuille profile for all ranks
    !----------------------------------------------------------------------------------------------------------
    ux_1c1(:) = ZERO
    call Generate_poiseuille_flow_profile (dm, ux_1c1)
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : to add profile to ux (default: x streamwise)
    !----------------------------------------------------------------------------------------------------------
    dtmp = dm%dpcc
    do i = 1, dtmp%xsz(1)
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        do k = 1, dtmp%xsz(3)
          fl%qx(i, j, k) =  fl%qx(i, j, k) + ux_1c1(jj)
        end do
      end do
    end do

    if(nrank == 0) Call Print_debug_mid_msg(" Maximum velocity for random velocities + given profile")
    call Find_maximum_velocity(dm, fl%qx, fl%qy, fl%qz)
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : Ensure the mass flow rate is 1.
    !----------------------------------------------------------------------------------------------------------
    !if(nrank == 0) call Print_debug_mid_msg("Ensure u, v, w, averaged in x and z direction is zero...")
    !call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), fl%fbcy_var(:, :, :, 1), dm, dm%dpcc, ux, ubulk, "ux")
    call Get_volumetric_average_3d_for_var_xcx(dm, dm%dpcc, fl%qx, ubulk, "ux")
    if(nrank == 0) then
      Call Print_debug_mid_msg("  The initial mass flux is:")
      write (*, wrtfmt1r) ' average[u(x,y,z)]_[x,y,z]: ', ubulk
    end if

    fl%qx(:, :, :) = fl%qx(:, :, :) / ubulk

    !call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), fl%fbcy_var(:, :, :, 1), dm, dm%dpcc, ux, ubulk, "ux")
    call Get_volumetric_average_3d_for_var_xcx(dm, dm%dpcc, fl%qx, ubulk, "ux")
    if(nrank == 0) then
      Call Print_debug_mid_msg("  The scaled mass flux is:")
      write (*, wrtfmt1r) ' average[u(x,y,z)]_[x,y,z]: ', ubulk
    end if

    ! to do : to add a scaling for turbulence generator inlet scaling, u = u * m / rho, check

    !----------------------------------------------------------------------------------------------------------
    !   some checking
    !----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qx, ux_ypencil, dm%dpcc)
    !----------------------------------------------------------------------------------------------------------
    !   Y-pencil : write out velocity profile
    !----------------------------------------------------------------------------------------------------------
    if(nrank == 0) then
      open ( newunit = pf_unit,     &
              file    = trim(dir_chkp)//'/check_poiseuille_profile.dat', &
              status  = 'replace',         &
              action  = 'write')
      write(pf_unit, '(A)') "# yc, ux_laminar, ux_real"
      do j = 1, dm%nc(2)
        write(pf_unit, '(5ES13.5)') dm%yc(j), ux_1c1(j), ux_ypencil(dm%dpcc%yen(1)/2, j, dm%dpcc%yen(3)/2)
      end do
      close(pf_unit)
    end if
    
    if(nrank == 0) call Print_debug_end_msg

    return
  end subroutine  Initialize_poiseuille_flow
  !==========================================================================================================
  !==========================================================================================================
  subroutine Initialize_flow_from_given_values(dm, fl)
    use udf_type_mod, only: t_domain
    use precision_mod, only: WP
    use parameters_constant_mod, only: ZERO
    use boundary_conditions_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow), intent(inout) :: fl
    
    if(nrank == 0) call Print_debug_mid_msg("Initializing flow field with given values...")
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : initial
    !----------------------------------------------------------------------------------------------------------
    fl%pres (:, :, :) = ZERO
    fl%qx(:, :, :) = ZERO
    fl%qy(:, :, :) = ZERO
    fl%qz(:, :, :) = ZERO
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : to get random fields [-1,1] for ux, uy, uz
    !----------------------------------------------------------------------------------------------------------
    call Generate_random_field(dm, fl)
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : update values
    !----------------------------------------------------------------------------------------------------------
    fl%qx(:, :, :) = fl%qx(:, :, :) + fl%init_velo3d(1)
    fl%qy(:, :, :) = fl%qy(:, :, :) + fl%init_velo3d(2)
    fl%qz(:, :, :) = fl%qz(:, :, :) + fl%init_velo3d(3)
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : apply b.c.
    !----------------------------------------------------------------------------------------------------------

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine

!==========================================================================================================
  !==========================================================================================================
  subroutine Initialize_flow_from_given_inlet(dm, fl)
    use udf_type_mod, only: t_domain
    use precision_mod, only: WP
    use parameters_constant_mod, only: ZERO
    use boundary_conditions_mod
    implicit none
    type(t_domain), intent(in)  :: dm
    type(t_flow), intent(inout) :: fl

    integer :: i, j, k, ii, jj, kk
    
    if(nrank == 0) call Print_debug_mid_msg("Initializing flow field with given profile...")
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : initial
    !----------------------------------------------------------------------------------------------------------
    fl%pres (:, :, :) = ZERO
    fl%qx(:, :, :) = ZERO
    fl%qy(:, :, :) = ZERO
    fl%qz(:, :, :) = ZERO
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : to get random fields [-1,1] for ux, uy, uz
    !----------------------------------------------------------------------------------------------------------
    call Generate_random_field(dm, fl)
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : update values
    !----------------------------------------------------------------------------------------------------------
    do k = 1, dm%dpcc%xsz(3)
      kk = dm%dpcc%xst(3) + k - 1
      do j = 1, dm%dpcc%xsz(2)
        jj = dm%dpcc%xst(2) + j - 1
        do i = 1, dm%dpcc%xsz(1)
          ii = dm%dpcc%xst(1) + i - 1
          fl%qx(i, j, k) = fl%qx(i, j, k) + fl%fbcx_qx(1, j, k)
        end do
      end do
    end do

    do k = 1, dm%dcpc%xsz(3)
      kk = dm%dcpc%xst(3) + k - 1
      do j = 1, dm%dcpc%xsz(2)
        jj = dm%dcpc%xst(2) + j - 1
        do i = 1, dm%dcpc%xsz(1)
          ii = dm%dcpc%xst(1) + i - 1
          fl%qy(i, j, k) = fl%qy(i, j, k) + fl%fbcx_qy(1, j, k)
        end do
      end do
    end do

    do k = 1, dm%dccp%xsz(3)
      kk = dm%dccp%xst(3) + k - 1
      do j = 1, dm%dccp%xsz(2)
        jj = dm%dccp%xst(2) + j - 1
        do i = 1, dm%dccp%xsz(1)
          ii = dm%dccp%xst(1) + i - 1
          fl%qz(i, j, k) = fl%qz(i, j, k) + fl%fbcx_qz(1, j, k)
        end do
      end do
    end do
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : apply b.c.
    !----------------------------------------------------------------------------------------------------------

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine

  !==========================================================================================================
  subroutine Initialize_flow_fields(dm, fl)
    use udf_type_mod
    use parameters_constant_mod
    use io_restart_mod
    use burgers_eq_mod
    use io_visulisation_mod
    use wtformat_mod
    use solver_tools_mod
    use continuity_eq_mod
    use boundary_conditions_mod
    implicit none

    type(t_domain), intent(inout) :: dm
    type(t_flow), intent(inout)   :: fl

    real(WP) :: velo(3)

    if(nrank == 0) call Print_debug_start_msg("Initialize flow fields ...")

  !----------------------------------------------------------------------------------------------------------
  ! to allocate flow variables
  !----------------------------------------------------------------------------------------------------------
    call Allocate_flow_variables (dm, fl)
  !----------------------------------------------------------------------------------------------------------
  ! to set up Re
  !----------------------------------------------------------------------------------------------------------
    call Update_Re(fl%iterfrom, fl)
  !----------------------------------------------------------------------------------------------------------
  ! initialize primary variables
  !----------------------------------------------------------------------------------------------------------
    fl%time = ZERO
    fl%iteration = 0

    if(fl%inittype == INIT_RESTART) then
      call read_instantanous_flow_raw_data(dm, fl)
      call restore_flow_variables_from_restart(dm, fl)

    else if (fl%inittype == INIT_INTERPL) then

    else if (fl%inittype == INIT_RANDOM) then
      call Generate_random_field(dm, fl)

    else if (fl%inittype == INIT_INLET) then
      call Initialize_flow_from_given_inlet(dm, fl)

    else if (fl%inittype == INIT_GVCONST) then
      call Initialize_flow_from_given_values(dm, fl)

    else if (fl%inittype == INIT_POISEUILLE) then
      call Initialize_poiseuille_flow(dm, fl)

    else if (fl%inittype == INIT_FUNCTION) then
      if (dm%icase == ICASE_TGV2D) then
        call Initialize_vortexgreen_2dflow (dm, fl)
      else if (dm%icase == ICASE_TGV3D) then
        call Initialize_vortexgreen_3dflow (dm, fl)
      else if (dm%icase == ICASE_BURGERS) then
        call Initialize_burgers_flow      (dm, fl)
      else
      end if
    else
    end if
!----------------------------------------------------------------------------------------------------------
! to initialize pressure correction term
!----------------------------------------------------------------------------------------------------------
    fl%pcor(:, :, :) = ZERO

    if(nrank == 0) call Print_debug_end_msg

    return
  end subroutine

  !==========================================================================================================
  subroutine Initialize_thermo_fields(dm, fl, tm)
    use udf_type_mod
    use parameters_constant_mod
    use eq_energy_mod
    use thermo_info_mod
    use io_restart_mod
    implicit none

    type(t_domain), intent(inout) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm

    integer :: i

    if(.not. dm%is_thermo) return
    if(nrank == 0) call Print_debug_start_msg("Initialize thermo fields ...")

!----------------------------------------------------------------------------------------------------------
! to allocate thermal variables
!----------------------------------------------------------------------------------------------------------
    call Allocate_thermo_variables (dm, tm)
!----------------------------------------------------------------------------------------------------------
! to set up Fr etc, require update flow Re first
!----------------------------------------------------------------------------------------------------------
    call Update_Re(fl%iterfrom, fl)
    call Update_PrGr(fl, tm) 
!----------------------------------------------------------------------------------------------------------
! initialize primary variables
!----------------------------------------------------------------------------------------------------------
    if(tm%inittype == INIT_RESTART) then
      call read_instantanous_thermo_raw_data  (dm, tm)
      call restore_thermo_variables_from_restart(dm, fl, tm)

    else if (tm%inittype == INIT_INTERPL) then
    else
      call Initialize_thermal_properties (tm)
      tm%time = ZERO
      tm%iteration = 0
    end if

    call Calculate_massflux_from_velocity (dm, fl, tm)
    ! to do, to check, scaling of gx? to be unified?
    !----------------------------------------------------------------------------------------------------------
    ! to set up old arrays 
    !----------------------------------------------------------------------------------------------------------
    tm%dDensm1(:, :, :) = tm%dDens(:, :, :)
    tm%dDensm2(:, :, :) = tm%dDens(:, :, :)
    
    return
  end subroutine




!==========================================================================================================
!==========================================================================================================
!> \brief Initialize Vortex Green flow
!>
!> This subroutine is called locally once.
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!> \param[out]    f             flow
!_______________________________________________________________________________
  subroutine  Initialize_vortexgreen_2dflow(dm, fl)
    use parameters_constant_mod!, only : HALF, ZERO, SIXTEEN, TWO
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in ) :: dm
    type(t_flow), intent(inout) :: fl
    real(WP) :: xc, yc
    real(WP) :: xp, yp
    integer :: i, j, ii, jj
    type(DECOMP_INFO) :: dtmp

    if(nrank == 0) call Print_debug_mid_msg("Initializing vortexgreen 2dflow ...")
!----------------------------------------------------------------------------------------------------------
!   ux in x-pencil
!----------------------------------------------------------------------------------------------------------
    dtmp = dm%dpcc
    do j = 1, dtmp%xsz(2)
      jj = dtmp%xst(2) + j - 1
      yc = dm%yc(jj)
      do i = 1, dtmp%xsz(1)
        ii = dtmp%xst(1) + i - 1
        xp = dm%h(1) * real(ii - 1, WP)
        fl%qx(i, j, :) =  sin_wp ( xp ) * cos_wp ( yc )
      end do
    end do
!----------------------------------------------------------------------------------------------------------
!   uy in x-pencil
!---------------------------------------------------------------------------------------------------------- 
    dtmp = dm%dcpc
    do j = 1, dtmp%xsz(2)
      jj = dtmp%xst(2) + j - 1
      yp = dm%yp(jj)
      do i = 1, dtmp%xsz(1)
        ii = dtmp%xst(1) + i - 1
        xc = dm%h(1) * (real(ii - 1, WP) + HALF)
        fl%qy(i, j, :) = -cos_wp ( xc ) * sin_wp ( yp )
      end do
    end do
!----------------------------------------------------------------------------------------------------------
!   uz in x-pencil
!---------------------------------------------------------------------------------------------------------- 
    fl%qz(:, :, :) =  ZERO
!----------------------------------------------------------------------------------------------------------
!   p in x-pencil
!----------------------------------------------------------------------------------------------------------
    fl%pres(:, :, :) =  ZERO
    ! dtmp = dm%dccc
    ! do j = 1, dtmp%xsz(2)
    !   jj = dtmp%xst(2) + j - 1
    !   yc = dm%yc(jj)
    !   do i = 1, dtmp%xsz(1)
    !     ii = dtmp%xst(1) + i - 1
    !     xc = dm%h(1) * (real(ii - 1, WP) + HALF)
    !     p(i, j, :)= ( cos_wp(TWO * xc) + sin(TWO * yc) ) * QUARTER
    !   end do
    ! end do
    
    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine Initialize_vortexgreen_2dflow
!==========================================================================================================
!==========================================================================================================
  subroutine  Validate_TGV2D_error(dm, fl)
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    !use iso_fortran_env
    implicit none

    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(in) :: fl
    integer :: k, i, j, ii, jj!, kk
    real(wp) :: uerr, ue, uc, verr, perr
    real(wp) :: xc, yc, xp, yp
    real(wp) :: uerrmax, verrmax, perrmax
    real(wp) :: perr_work, perrmax_work
    real(wp) :: uerr_work, uerrmax_work
    real(wp) :: verr_work, verrmax_work

    type(DECOMP_INFO) :: dtmp
    character( len = 128) :: filename
    logical :: file_exists = .FALSE.
    integer :: outputunit

!----------------------------------------------------------------------------------------------------------
!   X-pencil : Find Max. error of ux
!----------------------------------------------------------------------------------------------------------
    if(nrank == 0) call Print_debug_mid_msg("Validat TGV2D error ...")

    dtmp = dm%dpcc
    uerr = ZERO
    uerrmax = ZERO
    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yc = dm%yc(jj)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xp = dm%h(1) * real(ii - 1, WP)
          uc = fl%qx(i, j, k)
          ue = sin_wp ( xp ) * cos_wp ( yc ) * exp(- TWO * fl%rre * fl%time)
          uerr = uerr + (uc - ue)**2
          if(abs_wp(uc - ue) > uerrmax) uerrmax = abs_wp(uc - ue)
        end do
      end do
    end do
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(uerr,    uerr_work,    1, MPI_REAL_WP, MPI_SUM, MPI_COMM_WORLD, ierror)
    call mpi_allreduce(uerrmax, uerrmax_work, 1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
    uerr_work = uerr_work / real(dm%np(1), wp) / real(dm%nc(2), wp) / real(dm%nc(3), wp)
    uerr_work = sqrt_wp(uerr_work)
!----------------------------------------------------------------------------------------------------------
!   X-pencil : Find Max. error of uy
!----------------------------------------------------------------------------------------------------------
    dtmp = dm%dcpc
    verr = ZERO
    verrmax = ZERO
    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yp = dm%yp(jj)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xc = dm%h(1) * (real(ii - 1, WP) + HALF)
          uc = fl%qy(i, j, k)
          ue = - cos_wp ( xc ) * sin_wp ( yp ) * exp(- TWO * fl%rre * fl%time)
          verr = verr + (uc - ue)**2
          if(abs_wp(uc - ue) > verrmax) verrmax = abs_wp(uc - ue)
        end do
      end do
    end do
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(verr,    verr_work,    1, MPI_REAL_WP, MPI_SUM, MPI_COMM_WORLD, ierror)
    call mpi_allreduce(verrmax, verrmax_work, 1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
    verr_work = verr_work / real(dm%nc(1), wp) / real(dm%np(2), wp) / real(dm%nc(3), wp)
    verr_work = sqrt_wp(verr_work)
!----------------------------------------------------------------------------------------------------------
!   X-pencil : Find Max. error of p
!----------------------------------------------------------------------------------------------------------
    dtmp = dm%dccc
    perr = ZERO
    perrmax = ZERO
    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yc = dm%yc(jj)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xc = dm%h(1) * (real(ii - 1, WP) + HALF)
          uc = fl%pres(i, j, k)
          ue = ( cos_wp ( TWO * xc ) + sin_wp ( TWO * yc ) ) * QUARTER * (exp(- TWO * fl%rre * fl%time))**2
          perr = perr + (uc - ue)**2
          if(abs_wp(uc - ue) > perrmax) perrmax = abs_wp(uc - ue)
        end do
      end do
    end do
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(perr,    perr_work,    1, MPI_REAL_WP, MPI_SUM, MPI_COMM_WORLD, ierror)
    call mpi_allreduce(perrmax, perrmax_work, 1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
    perr_work = perr_work / real(dm%nc(1), wp) / real(dm%nc(2), wp) / real(dm%nc(3), wp)
    perr_work = sqrt_wp(perr_work)
!----------------------------------------------------------------------------------------------------------
!   X-pencil : write data in rank=0
!----------------------------------------------------------------------------------------------------------
    if(nrank == 0) then
      filename = 'Validation_TGV2d.dat'
      INQUIRE(FILE = trim(filename), exist = file_exists)
      if(.not.file_exists) then
        open(newunit = outputunit, file = trim(filename), action = "write", status = "new")
        write(outputunit, '(A)') 'Time, SD(u), SD(v), SD(p)'
      else
        open(newunit = outputunit, file = trim(filename), action = "write", status = "old", position="append")
      end if
      write(outputunit, '(1F10.4, 6ES17.7E3)') fl%time, uerr_work, verr_work, perr_work, &
            uerrmax_work, verrmax_work, perrmax_work
      close(outputunit)
    end if

    if(nrank == 0) call Print_debug_end_msg

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
!> \brief Initialize Vortex Green flow
!>
!> This subroutine is called locally once.
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!> \param[out]    f             flow
!_______________________________________________________________________________
  subroutine  Initialize_vortexgreen_3dflow(dm, fl)
    use parameters_constant_mod!, only : HALF, ZERO, SIXTEEN, TWO, PI
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in ) :: dm
    type(t_flow), intent(inout) :: fl
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    integer :: i, j, k, ii, jj, kk
    type(DECOMP_INFO) :: dtmp

    if(nrank == 0) call Print_debug_mid_msg("Initializing Taylor Green Vortex flow field ...")
!----------------------------------------------------------------------------------------------------------
!   ux in x-pencil
!---------------------------------------------------------------------------------------------------------- 
    dtmp = dm%dpcc
    do k = 1, dtmp%xsz(3)
      kk = dtmp%xst(3) + k - 1
      zc = dm%h(3) * (real(kk - 1, WP) + HALF)
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yc = dm%yc(jj)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xp = dm%h(1) * real(ii - 1, WP)
          fl%qx(i, j, k) =  sin_wp ( xp ) * cos_wp ( yc ) * cos_wp ( zc )
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
!   uy in x-pencil
!----------------------------------------------------------------------------------------------------------
    dtmp = dm%dcpc
    do k = 1, dtmp%xsz(3)
      kk = dtmp%xst(3) + k - 1
      zc = dm%h(3) * (real(kk - 1, WP) + HALF)
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yp = dm%yp(jj)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xc = dm%h(1) * (real(ii - 1, WP) + HALF)
          fl%qy(i, j, k) = -cos_wp ( xc ) * sin_wp ( yp ) * cos_wp ( zc )
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
!   uz in x-pencil
!---------------------------------------------------------------------------------------------------------- 
    !uz(:, :, :) =  ZERO
    dtmp = dm%dccp
    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
        do i = 1, dtmp%xsz(1)
          fl%qz(i, j, k) = zero
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
!   p in x-pencil
!----------------------------------------------------------------------------------------------------------
    dtmp = dm%dccc
    do k = 1, dtmp%xsz(3)
      kk = dtmp%xst(3) + k - 1
      zc = dm%h(3) * (real(kk - 1, WP) + HALF)
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yc = dm%yc(jj)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xc = dm%h(1) * (real(ii - 1, WP) + HALF)
          fl%pres(i, j, k)= ONE / SIXTEEN * ( cos(TWO * xc) + cos(TWO * yc) ) * &
                      (cos(TWO * zc) + TWO)
        end do
      end do
    end do

    if(nrank == 0) call Print_debug_end_msg
    
    return
  end subroutine Initialize_vortexgreen_3dflow

end module flow_thermo_initialiasation