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
module flow_thermo_initialisation
  use vars_df_mod
  use solver_tools_mod
  implicit none

  private  :: Allocate_flow_variables
  private  :: Allocate_thermo_variables
  private  :: Allocate_vof_variables
  private :: Generate_poiseuille_flow_profile
  private :: Generate_random_field

  private :: Initialize_poiseuille_flow
  private :: Initialize_flow_from_given_values
  private :: Initialize_vortexgreen_2dflow
  private :: Initialize_vortexgreen_3dflow
  private :: Initialize_rotational_flow
  private :: Initialize_vof_2d
  private :: Initialize_vof_3d

  public  :: Validate_TGV2D_error
  public  :: Initialize_flow_fields
  public  :: Initialize_thermo_fields
  public  :: Initialize_vof_fields
  public  :: Initialize_update_vortex_flow

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
  subroutine Allocate_flow_variables (fl, dm)
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
    call alloc_x(fl%qx,      dm%dpcc) ; fl%qx = ZERO
    call alloc_x(fl%qy,      dm%dcpc) ; fl%qy = ZERO
    call alloc_x(fl%qz,      dm%dccp) ; fl%qz = ZERO

    call alloc_x(fl%pres,    dm%dccc) ; fl%pres = ZERO
    call alloc_x(fl%pres0,    dm%dccc) ; fl%pres0 = ZERO
    call alloc_x(fl%pcor,    dm%dccc) ; fl%pcor = ZERO
    call alloc_z(fl%pcor_zpencil_ggg,    dm%dccc, .true.) ; fl%pcor_zpencil_ggg = ZERO

    call alloc_x(fl%mx_rhs,  dm%dpcc) ; fl%mx_rhs = ZERO
    call alloc_x(fl%my_rhs,  dm%dcpc) ; fl%my_rhs = ZERO
    call alloc_x(fl%mz_rhs,  dm%dccp) ; fl%mz_rhs = ZERO

    call alloc_x(fl%mx_rhs0, dm%dpcc) ; fl%mx_rhs0 = ZERO
    call alloc_x(fl%my_rhs0, dm%dcpc) ; fl%my_rhs0 = ZERO
    call alloc_x(fl%mz_rhs0, dm%dccp) ; fl%mz_rhs0 = ZERO

    if(dm%is_thermo) then
      call alloc_x(fl%gx,      dm%dpcc) ; fl%gx = ZERO
      call alloc_x(fl%gy,      dm%dcpc) ; fl%gy = ZERO
      call alloc_x(fl%gz,      dm%dccp) ; fl%gz = ZERO
      call alloc_x(fl%dDens,   dm%dccc) ; fl%dDens = ONE
      call alloc_x(fl%mVisc,   dm%dccc) ; fl%mVisc = ONE

      call alloc_x(fl%dDensm1, dm%dccc) ; fl%dDensm1 = ONE
      call alloc_x(fl%dDensm2, dm%dccc) ; fl%dDensm2 = ONE
    end if

    if(nrank == 0) call Print_debug_end_msg
    return

  end subroutine Allocate_flow_variables
  !==========================================================================================================
  subroutine Allocate_thermo_variables (tm, dm)
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

    if(nrank == 0) call Print_debug_end_msg
    return

  end subroutine Allocate_thermo_variables
  !==========================================================================================================
  subroutine Allocate_vof_variables (vf, dm)
    use parameters_constant_mod
    use mpi_mod
    implicit none

    type(t_domain), intent(in)    :: dm
    type(t_vof), intent(inout) :: vf

    if(.not. dm%is_vof) return
    if(nrank == 0) call Print_debug_start_msg("Allocating vof variables ...")
    !----------------------------------------------------------------------------------------------------------
    ! default : x pencil. 
    ! varaible index is LOCAL. means 1:xsize(1)
    !----------------------------------------------------------------------------------------------------------
    call alloc_x(vf%phi,     dm%dccc) ; vf%phi    = ZERO
    call alloc_x(vf%phi1,    dm%dccc) ; vf%phi1   = ZERO
    call alloc_x(vf%phi2,    dm%dccc) ; vf%phi2   = ZERO
    call alloc_x(vf%phi3,    dm%dccc) ; vf%phi3   = ZERO

    call alloc_x(vf%rhostar,  dm%dccc) ; vf%rhostar = ONE
    call alloc_x(vf%mustar,   dm%dccc) ; vf%mustar  = ONE

    call alloc_x(vf%lnx,     dm%dccc) ; vf%lnx    = ZERO
    call alloc_x(vf%lny,     dm%dccc) ; vf%lny    = ZERO
    call alloc_x(vf%lnz,     dm%dccc) ; vf%lnz    = ZERO

    call alloc_x(vf%llxx,    dm%dccc) ; vf%llxx   = ZERO
    call alloc_x(vf%llyy,    dm%dccc) ; vf%llyy   = ZERO
    call alloc_x(vf%llzz,    dm%dccc) ; vf%llzz   = ZERO
    call alloc_x(vf%llxy,    dm%dccc) ; vf%llxy   = ZERO
    call alloc_x(vf%llyz,    dm%dccc) ; vf%llyz   = ZERO
    call alloc_x(vf%llxz,    dm%dccc) ; vf%llxz   = ZERO

    call alloc_x(vf%a200,    dm%dccc) ; vf%a200   = ZERO
    call alloc_x(vf%a020,    dm%dccc) ; vf%a020   = ZERO
    call alloc_x(vf%a002,    dm%dccc) ; vf%a002   = ZERO
    call alloc_x(vf%a110,    dm%dccc) ; vf%a110   = ZERO
    call alloc_x(vf%a011,    dm%dccc) ; vf%a011   = ZERO
    call alloc_x(vf%a101,    dm%dccc) ; vf%a101   = ZERO
    call alloc_x(vf%a100,    dm%dccc) ; vf%a100   = ZERO
    call alloc_x(vf%a010,    dm%dccc) ; vf%a010   = ZERO
    call alloc_x(vf%a001,    dm%dccc) ; vf%a001   = ZERO
    call alloc_x(vf%dd,      dm%dccc) ; vf%dd     = ZERO

    call alloc_x(vf%kappa,   dm%dccc) ; vf%kappa  = ZERO

    call alloc_x(vf%flx,     dm%dpcc) ; vf%flx    = ZERO
    call alloc_x(vf%gly,     dm%dcpc) ; vf%gly    = ZERO
    call alloc_x(vf%hlz,     dm%dccp) ; vf%hlz    = ZERO

    if(nrank == 0) call Print_debug_end_msg
    return

  end subroutine Allocate_vof_variables
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
  subroutine Generate_random_field(dm, ux, uy, uz, p, lnoise)
    use random_number_generation_mod
    use parameters_constant_mod
    use mpi_mod
    use math_mod
    use boundary_conditions_mod
    implicit none
    type(t_domain), intent(in) :: dm
    real(WP),       intent(in) :: lnoise
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent(inout) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent(inout) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent(inout) :: uz
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent(inout) :: p
    
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
    p(:, :, :) = ZERO
    ux(:, :, :) = ZERO
    uy(:, :, :) = ZERO
    uz(:, :, :) = ZERO
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
             lownoise = lnoise
          ! end if
          do k = 1, dtmp%xsz(3)
            kk = dtmp%xst(3) + k - 1
            seed = (nrank + 1) * nsz * n + ii * 313 + jj * 571 + kk * 937
            !seed = ii + jj + kk + nsz * (n - 1)
            !write(*,*) ii, jj, kk, seed
            call Initialize_random_number ( seed )
            call Generate_r_random( -ONE, ONE, rd)
            if(n == 1) ux(i, j, k) = lownoise * rd
            if(n == 2) uy(i, j, k) = lownoise * rd
            if(n == 3) uz(i, j, k) = lownoise * rd

          end do
        end do
      end do
    end do

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
  subroutine Initialize_poiseuille_flow(dm, ux, uy, uz, p, lnoise)
    use input_general_mod
    use udf_type_mod
    use boundary_conditions_mod
    use parameters_constant_mod
    use wtformat_mod
    use files_io_mod
    implicit none
    type(t_domain), intent(in) :: dm
    real(WP),       intent(in) :: lnoise   
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent(out) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent(out) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent(out) :: uz
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent(out) ::  p
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
    p (:, :, :) = ZERO
    ux(:, :, :) = ZERO
    uy(:, :, :) = ZERO
    uz(:, :, :) = ZERO
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : to get random fields [-1,1] for ux, uy, uz
    !----------------------------------------------------------------------------------------------------------
    call Generate_random_field(dm, ux, uy, uz, p, lnoise)
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
          ux(i, j, k) =  ux(i, j, k) + ux_1c1(jj)
        end do
      end do
    end do

    if(nrank == 0) Call Print_debug_mid_msg(" Maximum velocity for random velocities + given profile")
    call Find_maximum_absvar3d(ux, "maximum ux:", wrtfmt1e)
    call Find_maximum_absvar3d(uy, "maximum uy:", wrtfmt1e)
    call Find_maximum_absvar3d(uz, "maximum uz:", wrtfmt1e)
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : Ensure the mass flow rate is 1.
    !----------------------------------------------------------------------------------------------------------
    !if(nrank == 0) call Print_debug_mid_msg("Ensure u, v, w, averaged in x and z direction is zero...")
    call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), dm%fbcy_qx(:, :, :), dm, dm%dpcc, ux, ubulk, "ux")
    if(nrank == 0) then
      Call Print_debug_mid_msg("  The initial mass flux is:")
      write (*, wrtfmt1r) ' average[u(x,y,z)]_[x,y,z]: ', ubulk
    end if

    ux(:, :, :) = ux(:, :, :) / ubulk

    call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), dm%fbcy_qx(:, :, :), dm, dm%dpcc, ux, ubulk, "ux")
    if(nrank == 0) then
      Call Print_debug_mid_msg("  The scaled mass flux is:")
      write (*, wrtfmt1r) ' average[u(x,y,z)]_[x,y,z]: ', ubulk
    end if
    if(nrank == 0) Call Print_debug_mid_msg(" Maximum velocity for after scaling unit mass flux")
    ! to do : to add a scaling for turbulence generator inlet scaling, u = u * m / rho

    !----------------------------------------------------------------------------------------------------------
    !   some checking
    !----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(ux, ux_ypencil, dm%dpcc)
    !----------------------------------------------------------------------------------------------------------
    !   Y-pencil : write out velocity profile
    !----------------------------------------------------------------------------------------------------------
    if(nrank == 0) then
      open ( newunit = pf_unit,     &
              file    = trim(dir_chkp)//'/check_poiseuille_profile.dat', &
              status  = 'replace',         &
              action  = 'write')
      write(pf_unit, '(A)') "#id,  yc, ux_laminar, ux_real"
      do j = 1, dm%nc(2)
        write(pf_unit, '(1I3.1, 5ES15.7)') j, dm%yc(j), ux_1c1(j), ux_ypencil(1, j, 1)
      end do
      close(pf_unit)
    end if
    
    if(nrank == 0) call Print_debug_end_msg

    return
  end subroutine  Initialize_poiseuille_flow
  !==========================================================================================================
  !==========================================================================================================
  subroutine Initialize_flow_from_given_values(dm, ux, uy, uz, p, lnoise, velo)
    use udf_type_mod, only: t_domain
    use precision_mod, only: WP
    use parameters_constant_mod, only: ZERO
    use boundary_conditions_mod
    implicit none
    type(t_domain), intent(in) :: dm
    real(WP),       intent(in) :: lnoise
    real(WP),       intent(in) :: velo(3)   
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent(out) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent(out) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent(out) :: uz
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent(out) ::  p
    
    if(nrank == 0) call Print_debug_mid_msg("Initializing flow field with given values...")
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : initial
    !----------------------------------------------------------------------------------------------------------
    p (:, :, :) = ZERO
    ux(:, :, :) = ZERO
    uy(:, :, :) = ZERO
    uz(:, :, :) = ZERO
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : to get random fields [-1,1] for ux, uy, uz
    !----------------------------------------------------------------------------------------------------------
    call Generate_random_field(dm, ux, uy, uz, p, lnoise)
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : update values
    !----------------------------------------------------------------------------------------------------------
    ux(:, :, :) = ux(:, :, :) + velo(1)
    uy(:, :, :) = uy(:, :, :) + velo(2)
    uz(:, :, :) = uz(:, :, :) + velo(3)
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : apply b.c.
    !----------------------------------------------------------------------------------------------------------

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine

!==========================================================================================================
  !==========================================================================================================
  subroutine Initialize_flow_from_given_inlet(dm, ux, uy, uz, p, lnoise)
    use udf_type_mod, only: t_domain
    use precision_mod, only: WP
    use parameters_constant_mod, only: ZERO
    use boundary_conditions_mod
    implicit none
    type(t_domain), intent(in) :: dm
    real(WP),       intent(in) :: lnoise
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent(out) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent(out) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent(out) :: uz
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent(out) ::  p

    integer :: i, j, k, ii, jj, kk
    
    if(nrank == 0) call Print_debug_mid_msg("Initializing flow field with given profile...")
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : initial
    !----------------------------------------------------------------------------------------------------------
    p (:, :, :) = ZERO
    ux(:, :, :) = ZERO
    uy(:, :, :) = ZERO
    uz(:, :, :) = ZERO
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : to get random fields [-1,1] for ux, uy, uz
    !----------------------------------------------------------------------------------------------------------
    call Generate_random_field(dm, ux, uy, uz, p, lnoise)
    !----------------------------------------------------------------------------------------------------------
    !   x-pencil : update values
    !----------------------------------------------------------------------------------------------------------
    do k = 1, dm%dpcc%xsz(3)
      kk = dm%dpcc%xst(3) + k - 1
      do j = 1, dm%dpcc%xsz(2)
        jj = dm%dpcc%xst(2) + j - 1
        do i = 1, dm%dpcc%xsz(1)
          ii = dm%dpcc%xst(1) + i - 1
          ux(i, j, k) = ux(i, j, k) + dm%fbcx_qx(1, j, k)
        end do
      end do
    end do

    do k = 1, dm%dcpc%xsz(3)
      kk = dm%dcpc%xst(3) + k - 1
      do j = 1, dm%dcpc%xsz(2)
        jj = dm%dcpc%xst(2) + j - 1
        do i = 1, dm%dcpc%xsz(1)
          ii = dm%dcpc%xst(1) + i - 1
          uy(i, j, k) = uy(i, j, k) + dm%fbcx_qy(1, j, k)
        end do
      end do
    end do

    do k = 1, dm%dccp%xsz(3)
      kk = dm%dccp%xst(3) + k - 1
      do j = 1, dm%dccp%xsz(2)
        jj = dm%dccp%xst(2) + j - 1
        do i = 1, dm%dccp%xsz(1)
          ii = dm%dccp%xst(1) + i - 1
          uz(i, j, k) = uz(i, j, k) + dm%fbcx_qz(1, j, k)
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
  subroutine Initialize_flow_fields(fl, dm)
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
    call Allocate_flow_variables (fl, dm)
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
      call read_instantanous_flow_raw_data(fl, dm)
      call restore_flow_variables_from_restart(fl, dm)

    else if (fl%inittype == INIT_INTERPL) then

    else if (fl%inittype == INIT_RANDOM) then
      call Generate_random_field(dm, fl%qx, fl%qy, fl%qz, fl%pres, fl%noiselevel)

    else if (fl%inittype == INIT_INLET) then
      call Initialize_flow_from_given_inlet(dm, fl%qx, fl%qy, fl%qz, fl%pres, fl%noiselevel)

    else if (fl%inittype == INIT_GVCONST) then
      velo(:) = fl%init_velo3d(:)
      call Initialize_flow_from_given_values(dm, fl%qx, fl%qy, fl%qz, fl%pres, fl%noiselevel, velo(:))

    else if (fl%inittype == INIT_POISEUILLE) then
      call Initialize_poiseuille_flow(dm, fl%qx, fl%qy, fl%qz, fl%pres, fl%noiselevel)

    else if (fl%inittype == INIT_FUNCTION) then
      if (dm%icase == ICASE_TGV2D) then
        call Initialize_vortexgreen_2dflow (dm, fl%qx, fl%qy, fl%qz, fl%pres)
      else if (dm%icase == ICASE_TGV3D) then
        call Initialize_vortexgreen_3dflow (dm, fl%qx, fl%qy, fl%qz, fl%pres)
      else if (dm%icase == ICASE_BURGERS) then
        call Initialize_burgers_flow      (dm, fl%qx, fl%qy, fl%qz, fl%pres)
      else if (dm%icase == ICASE_ROTATE ) then
        call Initialize_rotational_flow   (dm, fl%qx, fl%qy, fl%qz, fl%pres)
      else if (dm%icase == ICASE_VORTEX ) then
        call Initialize_update_vortex_flow (dm, fl%time, fl%qx, fl%qy, fl%qz, fl%pres)
      else
      end if
    else
    end if
!----------------------------------------------------------------------------------------------------------
! to initialize pressure correction term
!----------------------------------------------------------------------------------------------------------
    fl%pcor(:, :, :) = ZERO

    !==========================================================================================================
    !  validation for each time step
    !==========================================================================================================
    call Find_maximum_absvar3d(fl%qx, "init maximum ux:", wrtfmt1e)
    call Find_maximum_absvar3d(fl%qy, "init maximum uy:", wrtfmt1e)
    call Find_maximum_absvar3d(fl%qz, "init maximum uz:", wrtfmt1e)
    call Check_mass_conservation(fl, dm, 'initialization') 

    call write_snapshot_flow(fl, dm)

    if(nrank == 0) call Print_debug_end_msg

    return
  end subroutine

  !==========================================================================================================
  subroutine Initialize_thermo_fields(tm, fl, dm)
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
    call Allocate_thermo_variables (tm, dm)
!----------------------------------------------------------------------------------------------------------
! to set up Fr etc, require update flow Re first
!----------------------------------------------------------------------------------------------------------
    call Update_Re(fl%iterfrom, fl)
    call Update_PrGr(fl, tm) 
!----------------------------------------------------------------------------------------------------------
! initialize primary variables
!----------------------------------------------------------------------------------------------------------
    if(tm%inittype == INIT_RESTART) then
      call read_instantanous_thermo_raw_data  (tm, dm)
      call restore_thermo_variables_from_restart(fl, tm, dm)

    else if (tm%inittype == INIT_INTERPL) then
    else
      call Initialize_thermal_properties (fl, tm)
      tm%time = ZERO
      tm%iteration = 0
    end if

    call Calculate_massflux_from_velocity (fl, dm)
    !----------------------------------------------------------------------------------------------------------
    ! to set up old arrays 
    !----------------------------------------------------------------------------------------------------------
    fl%dDensm1(:, :, :) = fl%dDens(:, :, :)
    fl%dDensm2(:, :, :) = fl%dDens(:, :, :)
    
    return
  end subroutine

  !==========================================================================================================
  subroutine Initialize_vof_fields(vf, fl, dm)
    use udf_type_mod
    use parameters_constant_mod
    use io_restart_mod
    use io_visulisation_mod
    use vof_info_mod
    use eq_vof_mod
    implicit none

    type(t_domain), intent(inout) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_vof), intent(inout) :: vf

    integer :: i

    if(.not. dm%is_vof) return
    if(nrank == 0) call Print_debug_start_msg("Initialize vof fields ...")

!----------------------------------------------------------------------------------------------------------
! to allocate vof variables
!----------------------------------------------------------------------------------------------------------
    call Allocate_vof_variables (vf, dm)
!----------------------------------------------------------------------------------------------------------
! initialize primary variables
!----------------------------------------------------------------------------------------------------------
    vf%time = ZERO
    vf%iteration = 0
    vf%voflim = 1.e-8_WP

    if(vf%inittype == INIT_RESTART) then
      call read_instantanous_vof_raw_data  (dm, vf)
      call Initialize_vof_properties (fl, vf)
      call restore_vof_variables_from_restart(dm, fl, vf)

    else if (vf%inittype == INIT_INTERPL) then
    
    else if (vf%inittype == INIT_SLOTDISK) then
      call Initialize_vof_2d      (dm, vf%inittype, vf%phi, vf%ibeta)
    else if (vf%inittype == INIT_SLOTSPHERE) then
      call Initialize_vof_3d    (dm, vf%inittype, vf%phi, vf%ibeta)
    else if (vf%inittype == INIT_2DBUB) then
      call Initialize_vof_2d  (dm, vf%inittype, vf%phi, vf%ibeta)
      call Initialize_vof_properties (fl, vf)
      vf%time = ZERO
      vf%iteration = 0
    else if (vf%inittype == INIT_3DBUB) then
      call Initialize_vof_3d  (dm, vf%inittype, vf%phi, vf%ibeta)
      call Initialize_vof_properties (fl, vf)
      vf%time = ZERO
      vf%iteration = 0
    end if

    call Update_Re(fl%iterfrom, fl)
    call Update_gravity_sigma(fl, vf)
    call cal_colour_function(dm, vf)

!    call Calculate_massflux_from_velocity (fl, dm)
    call write_snapshot_vof(vf, dm)

!    fl%dDensm1(:, :, :) = fl%dDens(:, :, :)
!    fl%dDensm2(:, :, :) = fl%dDens(:, :, :)

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
  subroutine  Initialize_vortexgreen_2dflow(dm, ux, uy, uz, p)
    use parameters_constant_mod!, only : HALF, ZERO, SIXTEEN, TWO
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in ) :: dm
    real(WP),       intent(inout) :: ux(:, :, :) , &
                                     uy(:, :, :) , &
                                     uz(:, :, :) , &
                                     p (:, :, :)
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
        ux(i, j, :) =  sin_wp ( xp ) * cos_wp ( yc )
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
        uy(i, j, :) = -cos_wp ( xc ) * sin_wp ( yp )
      end do
    end do
!----------------------------------------------------------------------------------------------------------
!   uz in x-pencil
!---------------------------------------------------------------------------------------------------------- 
    uz(:, :, :) =  ZERO
!----------------------------------------------------------------------------------------------------------
!   p in x-pencil
!----------------------------------------------------------------------------------------------------------
    p(:, :, :) =  ZERO
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
!> \brief Initialize rotating velocity field (for the slotted disk/sphere test)
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
  subroutine  Initialize_rotational_flow (dm, ux, uy, uz, p)
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in ) :: dm
    real(WP),       intent(inout) :: ux(:, :, :) , &
                                     uy(:, :, :) , &
                                     uz(:, :, :) , &
                                     p (:, :, :)
    real(WP) :: xc, yc
    real(WP) :: xp, yp
    integer :: i, j, ii, jj
    type(DECOMP_INFO) :: dtmp

    if(nrank == 0) call Print_debug_mid_msg("Initializing rotational flow ...")
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
        ux(i, j, :) =  real(0.5, WP) - yc
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
        uy(i, j, :) = xc - real(0.5, WP)
      end do
    end do
!----------------------------------------------------------------------------------------------------------
!   uz in x-pencil
!---------------------------------------------------------------------------------------------------------- 
    uz(:, :, :) =  ZERO
!----------------------------------------------------------------------------------------------------------
!   p in x-pencil
!----------------------------------------------------------------------------------------------------------
    p(:, :, :) =  ZERO

    if(nrank == 0) call Print_debug_end_msg

    return

  end subroutine Initialize_rotational_flow
!==========================================================================================================
!==========================================================================================================
!> \brief Initialize a time dependent vortex flow (for the stretched sphere test)
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
  subroutine  Initialize_update_vortex_flow (dm, time, ux, uy, uz, p)
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in ) :: dm
    real(WP),       intent(in ) :: time
    real(WP),       intent(inout) :: ux(:, :, :) , &
                                     uy(:, :, :) , &
                                     uz(:, :, :) , &
                                     p (:, :, :)
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    real(WP), parameter :: ct = 3.
    integer :: i, j, k, ii, jj, kk
    type(DECOMP_INFO) :: dtmp

    if(nrank == 0) call Print_debug_mid_msg("Initializing vortex flow ...")
!----------------------------------------------------------------------------------------------------------
!   ux in x-pencil
!----------------------------------------------------------------------------------------------------------
    dtmp = dm%dpcc
    do k = 1, dtmp%xsz(3)
      kk = dtmp%xst(3) + k - 1
      zc = dm%h(3) * (real(kk-1, WP) + HALF)
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yc = dm%yc(jj)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xp = dm%h(1) * real(ii-1, WP)
          ux(i, j, k) =  TWO*(sin(PI*xp)**2)*sin(TWO*PI*yc)*sin(TWO*PI*zc)*cos(PI*time/ct)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
!   uy in x-pencil
!---------------------------------------------------------------------------------------------------------- 
    dtmp = dm%dcpc
    do k = 1, dtmp%xsz(3)
      kk = dtmp%xst(3) + k - 1
      zc = dm%h(3) * (real(kk-1, WP) + HALF)
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yp = dm%yp(jj)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xc = dm%h(1) * (real(ii-1, WP) + HALF)
          uy(i, j, k) = -sin(TWO*PI*xc)*(sin(PI*yp)**2)*sin(TWO*PI*zc)*cos(PI*time/ct)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
!   uz in x-pencil
!---------------------------------------------------------------------------------------------------------- 
    dtmp = dm%dccp
    do k = 1, dtmp%xsz(3)
      kk = dtmp%xst(3) + k - 1
      zp = dm%h(3) * real(kk-1, WP)
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yc = dm%yc(jj)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xc = dm%h(1) * (real(ii - 1, WP) + HALF)
          uz(i, j, k) = -sin(TWO*PI*xc)*sin(TWO*PI*yc)*(sin(PI*zp)**2)*cos(PI*time/ct)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
!   p in x-pencil
!----------------------------------------------------------------------------------------------------------
    p(:, :, :) =  ZERO

    if(nrank == 0) call Print_debug_end_msg

    return

  end subroutine Initialize_update_vortex_flow
!==========================================================================================================
!==========================================================================================================
!> \brief Initialize vof field of a rotating slotted disk 
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
  subroutine  Initialize_vof_2d(dm, icase, phi, ibeta)
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(IN) :: dm
    integer, intent(IN) :: icase
    integer, intent(IN) :: ibeta
    real(WP),       intent(INOUT) :: phi(:, :, :)
    real(WP) :: xp1, yp1, xp2, yp2
    integer :: i, j, ii, jj
    type(DECOMP_INFO) :: dtmp

    real(WP) :: hh, beta
    real(WP) :: dx, dy
    integer :: ires_x, ires_y

    real(WP), parameter :: eps = 1.d-30

    if(nrank == 0) call Print_debug_mid_msg("Initializing 2d vof field ...")
!----------------------------------------------------------------------------------------------------------
!   vof in x-pencil
!----------------------------------------------------------------------------------------------------------

    beta = real(ibeta, WP)

    ires_x = 20
    ires_y = 20

    dtmp = dm%dccc
    do j = 1, dtmp%xsz(2)
      jj = dtmp%xst(2) + j - 1
      yp1 = dm%yp(jj)
      yp2 = dm%yp(jj+1)
      do i = 1, dtmp%xsz(1)
        ii = dtmp%xst(1) + i - 1
        xp1 = dm%h(1) * real(ii - 1, WP)
        xp2 = dm%h(1) * real(ii, WP)
        hh = xp2-xp1
        dx = (xp2-xp1)/real(ires_x, WP)
        dy = (yp2-yp1)/real(ires_y, WP)
        phi(i, j, :) =  CAL_VOF(xp1, yp1, xp2, yp2, ires_x, ires_y)
      end do
    end do

    if(nrank == 0) call Print_debug_end_msg

    return

    contains

    function CAL_VOF(x1, y1, x2, y2, ires_x, ires_y) result (vof)

      implicit none

      real(WP), intent(IN) :: x1, y1, x2, y2
      integer, intent(IN) :: ires_x, ires_y
      real(WP) :: vof

      real(WP) :: vol
      real(WP), allocatable, dimension(:,:) :: xx, yy
      real(WP), allocatable, dimension(:,:) :: psi, colorfunc
      real(WP), allocatable, dimension(:) :: temp
      real(WP), allocatable, dimension(:) :: x,y

      integer :: i, j

      allocate (xx(ires_x+1,ires_y+1))
      allocate (yy(ires_x+1,ires_y+1))
      allocate (psi(ires_x+1,ires_y+1))
      allocate (colorfunc(ires_x+1,ires_y+1))
      allocate (temp(ires_y+1))

      xx(1,:) = x1
      yy(:,1) = y1
      do i=1, ires_x
        xx(i+1,:) = xx(i,:)+dx
      enddo
      do j=1, ires_y
        yy(:,j+1) = yy(:,j)+dy
      enddo
      do j=1, ires_y+1
        do i=1, ires_x+1
          if(icase==INIT_SLOTDISK) then
            psi(i,j) = CAL_PSI_SLOTDISK(xx(i,j), yy(i,j))
          else if(icase==INIT_2DBUB)then
            psi(i,j) = CAL_PSI_2DBUB(xx(i,j), yy(i,j))
          end if
          colorfunc(i,j) = 0.5d0*(1.d0+tanh(beta/hh*psi(i,j)))
        enddo
      enddo

      vol = abs((y2-y1)*(x2-x1))

      do j=1, ires_y+1
        temp(j) = SIMPSON(xx(:,j), colorfunc(:,j))
      enddo

      vof = SIMPSON(yy(1,:), temp)/vol

      deallocate(xx, yy, psi, colorfunc, temp)

    end function CAL_VOF

    function CAL_PSI_2DBUB(x, y) result (psi)

      implicit none

      real(WP), intent(IN) :: x, y
      real(WP) :: psi
      real(WP) :: radi, xc0, yc0
      real(WP) :: x1, y1

      radi = 0.25d0
      xc0 = 0.5d0
      yc0 = 0.5d0

      psi = radi-sqrt((x-xc0)**2+(y-yc0)**2)

    end function CAL_PSI_2DBUB

    function CAL_PSI_SLOTDISK(x, y) result (psi)

      implicit none

      real(WP), intent(IN) :: x, y
      real(WP) :: psi

      real(WP) :: radi, xo, yo, xc0, yc0, x110, y110, x120, y120, x210, y210,  &
                  x220, y220
      real(WP) :: psi1, psi2, psi3, psi4
      real(WP) :: x1, y1
      real(WP) :: d11, d12, d13, d21, d22, d23, d31, d32, d33

      radi = 0.15d0
      xo = 0.5d0
      yo = 0.5d0
      xc0 = 0.5d0
      yc0 = 0.75d0
      x110 = xc0-0.025d0
      y110 = -sqrt(radi**2-0.025d0**2)+yc0
      x120 = xc0-0.025d0
      y120 = 0.85d0
      x210 = xc0+0.025d0
      y210 = -sqrt(radi**2-0.025d0**2)+yc0
      x220 = xc0+0.025d0
      y220 = 0.85d0

      x1 = x-xo
      y1 = y-yo
      psi1 = radi-sqrt((x-xc0)**2+(y-yc0)**2)
      if (x1>=x120-xo .and. x1<=x220-xo .and. y1<=y120-yo) psi1 = 1.d0
      d11 = (x-x110)*(y120-y110)-(x120-x110)*(y-y110)
      d11 = d11/sqrt((x120-x110)**2+(y120-y110)**2+eps)
      d12 = sqrt((x-x110)**2+(y-y110)**2)
      d13 = sqrt((x-x120)**2+(y-y120)**2)
      if (y1>=y110-yo .and. y1<=y120-yo) then
        psi2 = d11
      else
        psi2 = min(d12,d13)
      endif
      d21 = (x-x210)*(y220-y210)-(x220-x210)*(y-y210)
      d21 = d21/sqrt((x220-x210)**2+(y220-y210)**2+eps)
      d22 = sqrt((x-x210)**2+(y-y210)**2)
      d23 = sqrt((x-x220)**2+(y-y220)**2)
      if (y1>=y210-yo .and. y1<=y220-yo)then
        psi3 = d21
      else
        psi3 = min(d22,d23)
      endif
      d31 = (x-x120)*(y220-y120)-(x220-x120)*(y-y120)
      d31 = d31/sqrt((x220-x120)**2+(y220-y120)**2+eps)
      d32 = sqrt((x-x120)**2+(y-y120)**2)
      d33 = sqrt((x-x220)**2+(y-y220)**2)
      if (x1>=x120-xo .and. x1<=x220-xo)then
        psi4 = d31
      else
        psi4 = min(d32,d33)
      endif

      if (psi1>=0.d0 .and. (d11<=0.d0 .or. d21>=0.d0 .or. d31<=0.d0))then
        psi = min(abs(psi1),abs(psi2),abs(psi3),abs(psi4))
      else
        psi = -min(abs(psi1),abs(psi2),abs(psi3),abs(psi4))
      endif

    end function CAL_PSI_SLOTDISK

    function SIMPSON(x, func) result(integral)

      implicit none

      real(WP), dimension(:), intent(IN) :: x, func
      real(WP) :: integral

      integer :: i, n
      real(WP) :: h

      n = size(x)-1
      if(mod(n,2)/=0) then
        if(nrank == 0) then
          write (*,*) 'Error: n is not a even number in SIMPSON, simulation &&
                       terminated'
        end if
        stop
      endif
      h = (x(n+1)-x(1))/real(n,kind=8)

      integral = 0.d0
      do i=3, n-1, 2
        integral = integral+2.d0*func(i)
      enddo

      integral = integral+4.d0*func(2)+4.d0*func(n)

      do i=4, n-2, 2
        integral = integral+4.d0*func(i)
      enddo

      integral = (h/3.d0)*(func(1)+func(n+1)+integral)

    end function SIMPSON

  end subroutine Initialize_vof_2d
!==========================================================================================================
!==========================================================================================================
!> \brief Initialize vof field of a rotating slotted sphere 
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
  subroutine  Initialize_vof_3d(dm, icase, phi, ibeta)
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(IN) :: dm
    integer, intent(IN) :: icase
    integer, intent(IN) :: ibeta
    real(WP),       intent(INOUT) :: phi(:, :, :)
    real(WP) :: xp1, yp1, xp2, yp2, zp1, zp2
    integer :: i, j, k, ii, jj, kk
    type(DECOMP_INFO) :: dtmp

    real(WP) :: hh, beta
    real(WP) :: dx, dy, dz
    integer :: ires_x, ires_y, ires_z

    real(WP), parameter :: eps = 1.e-30_WP

    if(nrank == 0) call Print_debug_mid_msg("Initializing 3d vof field ...")
!----------------------------------------------------------------------------------------------------------
!   vof in x-pencil
!----------------------------------------------------------------------------------------------------------

    beta = real(ibeta, WP)

    ires_x = 20
    ires_y = 20
    ires_z = 20

    dtmp = dm%dccc
    do k = 1, dtmp%xsz(3)
      kk = dtmp%xst(3) + k - 1
      zp1 = dm%h(3) * real(kk - 1, WP)
      zp2 = dm%h(3) * real(kk, WP)
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yp1 = dm%yp(jj)
        yp2 = dm%yp(jj+1)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xp1 = dm%h(1) * real(ii - 1, WP)
          xp2 = dm%h(1) * real(ii, WP)
          hh = xp2-xp1
          dx = (xp2-xp1)/real(ires_x, WP)
          dy = (yp2-yp1)/real(ires_y, WP)
          dz = (zp2-zp1)/real(ires_z, WP)
          phi(i, j, k) =  CAL_VOF(xp1, yp1, zp1, xp2, yp2, zp2, &
                                  ires_x, ires_y, ires_z)
        end do
      end do
    end do

    if(nrank == 0) call Print_debug_end_msg

    return

    contains

    function CAL_VOF(x1, y1, z1, x2, y2, z2, ires_x, ires_y, ires_z) result (vof)

      implicit none

      real(WP), intent(IN) :: x1, y1, z1, x2, y2, z2
      integer, intent(IN) :: ires_x, ires_y, ires_z
      real(WP) :: vof

      real(WP) :: vol
      real(WP), allocatable, dimension(:,:,:) :: xx, yy, zz
      real(WP), allocatable, dimension(:,:,:) :: psi, colorfunc
      real(WP), allocatable, dimension(:,:) :: temp1
      real(WP), allocatable, dimension(:) :: temp2
      real(WP), allocatable, dimension(:) :: x,y

      integer :: i, j, k

      allocate (xx(ires_x+1,ires_y+1,ires_z+1))
      allocate (yy(ires_x+1,ires_y+1,ires_z+1))
      allocate (zz(ires_x+1,ires_y+1,ires_z+1))
      allocate (psi(ires_x+1,ires_y+1,ires_z+1))
      allocate (colorfunc(ires_x+1,ires_y+1,ires_z+1))
      allocate (temp1(ires_y+1,ires_z+1))
      allocate (temp2(ires_z+1))

      xx(1,:,:) = x1
      yy(:,1,:) = y1
      zz(:,:,1) = z1
      do i=1, ires_x
        xx(i+1,:,:) = xx(i,:,:)+dx
      enddo
      do j=1, ires_y
        yy(:,j+1,:) = yy(:,j,:)+dy
      enddo
      do k=1, ires_z
        zz(:,:,k+1) = zz(:,:,k)+dz
      enddo
      do k=1, ires_z+1
        do j=1, ires_y+1
          do i=1, ires_x+1
            if(icase==INIT_SLOTSPHERE) then
              psi(i,j,k) = CAL_PSI_SLOTSPHERE(xx(i,j,k), yy(i,j,k), zz(i,j,k))
            else if(icase==INIT_3DBUB) then
              psi(i,j,k) = CAL_PSI_3DBUB(xx(i,j,k), yy(i,j,k), zz(i,j,k))
            end if
            colorfunc(i,j,k) = 0.5_WP*(1._WP+tanh(beta/hh*psi(i,j,k)))
          enddo
        enddo
      enddo

      vol = abs((z2-z1)*(y2-y1)*(x2-x1))

      do k=1, ires_z+1
        do j=1, ires_y+1
          temp1(j,k) = SIMPSON(xx(:,j,k), colorfunc(:,j,k))
        enddo
      enddo

      do k=1, ires_z+1
        temp2(k) = SIMPSON(yy(1,:,k), temp1(:,k))
      enddo

      vof = SIMPSON(zz(1,1,:), temp2)/vol

      deallocate(xx, yy, zz, psi, colorfunc, temp1, temp2)

    end function CAL_VOF

    function CAL_PSI_3DBUB(x, y, z) result (psi)

      implicit none

      real(WP), intent(IN) :: x, y, z
      real(WP) :: psi

      real(WP) :: radi, xc0, yc0, zc0

!      radi = 1.0_WP
!      xc0 = 4.0_WP
!      yc0 = 4.0_WP
!      zc0 = 2.0_WP

      radi = 0.25_WP
      xc0 = 0.5_WP
      yc0 = 0.5_WP
      zc0 = 0.5_WP

      psi = radi-sqrt((x-xc0)**2+(y-yc0)**2+(z-zc0)**2)

    end function CAL_PSI_3DBUB

    function CAL_PSI_SLOTSPHERE(x, y, z) result (psi)

      implicit none

      real(WP), intent(IN) :: x, y, z
      real(WP) :: psi

      real(WP) :: radi, r1, r2, xc0, yc0, zc0,           &
                  x110, y110, x120, y120, x210, y210,    &
                  x220, y220, z121, z122
      real(WP) :: psi1, psi2, psi3, psi4
      real(WP) :: d1, d2, d3, d4, k

      radi = 0.15_WP
      xc0 = 0.5_WP
      yc0 = 0.75_WP
      zc0 = 0.25_WP
      r1 = sqrt(radi**2-0.025_WP**2)
      r2 = sqrt(radi**2-0.1_WP**2)
      x110 = xc0-0.025_WP
      y110 = -r1+yc0
      x120 = xc0-0.025_WP
      y120 = 0.85_WP
      z121 = -sqrt(r1**2-(y120-yc0)**2)+zc0
      z122 = sqrt(r1**2-(y120-yc0)**2)+zc0
      x210 = xc0+0.025_WP
      y210 = -r1+yc0
      x220 = xc0+0.025_WP
      y220 = 0.85_WP

      d1 = radi-sqrt((x-xc0)**2+(y-yc0)**2+(z-zc0)**2)
      psi1 = d1

      if (x>=x120 .and. x<=x220 .and. y<=y120) psi1 = 1._WP

      k = (z122-zc0)/(y120-yc0)
      if ((z-zc0)**2+(y-yc0)**2<=r1**2 .and. y<=y120) then
        d2 = (x-x110)*(y120-y110)-(x120-x110)*(y-y110)
        d2 = d2/sqrt((x120-x110)**2+(y120-y110)**2+eps)
        psi2 = d2
      else if ((y-yc0)>=k*(z-zc0) .and. (y-yc0)>=-k*(z-zc0) .and. y>y120) then
        psi2 = CAL_LINE_DIST(z,y,x-x110,z121,y120,z122,y120)
      else
        psi2 = CAL_CIRC_DIST(z,y,x-x110,zc0,yc0,r1)
      endif

      if ((z-zc0)**2+(y-yc0)**2<=r1**2 .and. y<=y120) then
        d3 = (x-x210)*(y220-y210)-(x220-x210)*(y-y210)
        d3 = d3/sqrt((x220-x210)**2+(y220-y210)**2+eps)
        psi3 = d3
      else if ((y-yc0)>=(z-zc0) .and. (y-yc0)>=-(z-zc0) .and. y>y120) then
        psi3 = CAL_LINE_DIST(z,y,x-x210,z121,y120,z122,y120)
      else
        psi3 = CAL_CIRC_DIST(z,y,x-x210,zc0,yc0,r1)
      endif

      k = (z122-zc0)/(x210-xc0)
      if ((x-xc0)**2+(z-zc0)**2<=r2**2 .and. x<=x210 .and. x>=x110) then
        d4 = (x-x120)*(y220-y120)-(x220-x120)*(y-y120)
        d4 = d4/sqrt((x220-x120)**2+(y220-y120)**2+eps)
        psi4 = d4 
      else if ((z-zc0)<=k*(x-xc0) .and. (z-zc0)>=-k*(x-xc0) .and. x>x210) then
        psi4 = CAL_LINE_DIST(x,z,y-y120,x210,z121,x210,z122)
      else if ((z-zc0)<=-k*(x-xc0) .and. (z-zc0)>=k*(x-xc0) .and. x<x110) then
        psi4 = CAL_LINE_DIST(x,z,y-y120,x110,z121,x110,z122)
      else
        psi4 = CAL_CIRC_DIST(x,z,y-y120,xc0,yc0,r2)
      endif   

!      if (d1>=0._WP .and. (d2<=0._WP .or. d3>=0._WP .or. d4<=0._WP))then
      if (d1>=0._WP .and. (x<=x120 .or. x>=x220 .or. y>=y120))then
        psi = min(abs(psi1),abs(psi2),abs(psi3),abs(psi4))
      else
        psi = -min(abs(psi1),abs(psi2),abs(psi3),abs(psi4))
      endif

    end function CAL_PSI_SLOTSPHERE

    function CAL_CIRC_DIST(x1,x2,norm,o1,o2,r) result(dist)

      implicit none

      real(WP), intent(IN) :: x1,x2,norm,o1,o2,r
      real(WP) :: dist

      dist = sqrt((x1-o1)**2+(x2-o2)**2)-r
      dist = sqrt(norm**2+dist**2)

    end function CAL_CIRC_DIST

    function CAL_LINE_DIST(x1,x2,norm,p11,p12,p21,p22) result(dist)

      implicit none

      real(WP), intent(IN) :: x1,x2,norm,p11,p12,p21,p22
      real(WP) :: dist

      real(WP) :: d1,d2

      d1 = 0._WP
      d2 = 0._WP
      if(x1>=p11 .and. x1<=p21)then
        dist = (x1-p11)*(p22-p12)-(p21-p11)*(x2-p12)
        dist = dist/sqrt((p21-p11)**2+(p22-p12)**2+eps)
        dist = sqrt(norm**2+dist**2)
      else
        d1 = sqrt((x1-p11)**2+(x2-p12)**2+norm**2)
        d2 = sqrt((x1-p21)**2+(x2-p22)**2+norm**2)
        dist = min(d1,d2)
      end if

    end function CAL_LINE_DIST

    function SIMPSON(x, func) result(integral)

      implicit none

      real(WP), dimension(:), intent(IN) :: x, func
      real(WP) :: integral

      integer :: i, n
      real(WP) :: h

      n = size(x)-1
      if(mod(n,2)/=0) then
        if(nrank == 0) then
          write (*,*) 'Error: n is not a even number in SIMPSON, simulation &&
                       terminated'
        end if
        stop
      endif
      h = (x(n+1)-x(1))/real(n,kind=8)

      integral = 0._WP
      do i=3, n-1, 2
        integral = integral+2._WP*func(i)
      enddo

      integral = integral+4._WP*func(2)+4._WP*func(n)

      do i=4, n-2, 2
        integral = integral+4._WP*func(i)
      enddo

      integral = (h/3._WP)*(func(1)+func(n+1)+integral)

    end function SIMPSON

  end subroutine Initialize_vof_3d
!==========================================================================================================
!==========================================================================================================
  subroutine  Validate_TGV2D_error(fl, dm)
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
  subroutine  Initialize_vortexgreen_3dflow(dm, ux, uy, uz, p)
    use parameters_constant_mod!, only : HALF, ZERO, SIXTEEN, TWO, PI
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in ) :: dm
    real(WP),       intent(inout) :: ux(:, :, :) , &
                                     uy(:, :, :) , &
                                     uz(:, :, :) , &
                                     p (:, :, :)
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
          ux(i, j, k) =  sin_wp ( xp ) * cos_wp ( yc ) * cos_wp ( zc )
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
          uy(i, j, k) = -cos_wp ( xc ) * sin_wp ( yp ) * cos_wp ( zc )
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
          uz(i, j, k) = zero
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
          p(i, j, k)= ONE / SIXTEEN * ( cos(TWO * xc) + cos(TWO * yc) ) * &
                      (cos(TWO * zc) + TWO)
        end do
      end do
    end do

    if(nrank == 0) call Print_debug_end_msg
    
    return
  end subroutine Initialize_vortexgreen_3dflow

end module flow_thermo_initialisation
