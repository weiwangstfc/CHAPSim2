!---------------------------------------------------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------------------------
module flow_thermo_initialiasation
  use vars_df_mod
  use solver_tools_mod
  implicit none

  
  private :: Generate_poiseuille_flow_profile
  private :: Initialize_poiseuille_flow
  private :: Initialize_vortexgreen_2dflow
  private :: Initialize_vortexgreen_3dflow
  private :: Generate_random_field

  private  :: Allocate_flow_variables
  private  :: Allocate_thermo_variables
  private  :: Initialize_flow_variables
  private  :: Initialize_thermo_variables

  public  :: Validate_TGV2D_error
  public  :: Initialize_flow_thermal_fields

contains
!=============================================================================================================================================
!> \brief Initialisation and preprocessing of the flow field
!--------------------------------------------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain   module
!>         all    once           all       public
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[in]     none          NA
!> \param[out]    none          NA
!=============================================================================================================================================
  subroutine Initialize_flow_thermal_fields
    use mpi_mod
    use vars_df_mod
    use thermo_info_mod
    use boundary_conditions_mod
    use decomp_2d_io
    use typeconvert_mod
    use restart_mod
    use mpi_mod
    use eq_momentum_mod
    use continuity_eq_mod
    implicit none

    integer :: i, j, k, l, s, iter
    type(t_fluidThermoProperty) :: ftpx, ftpy, ftpz
    integer :: ipencil
    type(DECOMP_INFO) :: dtmp

    if(nrank == 0) call Print_debug_start_msg("Initialize flow and thermal fields ...")
    iter = 0
    do l = 1, nxdomain

      !---------------------------------------------------------------------------------------------------------------------------------------------
      ! initialize common thermal variables
      !---------------------------------------------------------------------------------------------------------------------------------------------
      if(domain(l)%ithermo == 0) then
        domain(l)%fbc_dend(:, :) = ONE
        domain(l)%fbc_vism(:, :) = ONE
      else 
        do s = 1, 2
          ftpx = thermo(l)%ftpbcx(s)
          ftpy = thermo(l)%ftpbcy(s)
          ftpz = thermo(l)%ftpbcz(s)

          domain(l)%fbc_dend(s, 1) = ftpx%d
          domain(l)%fbc_dend(s, 2) = ftpy%d
          domain(l)%fbc_dend(s, 3) = ftpz%d
          domain(l)%fbc_vism(s, 1) = ftpx%m
          domain(l)%fbc_vism(s, 2) = ftpy%m
          domain(l)%fbc_vism(s, 3) = ftpz%m
        end do
      end if
      !---------------------------------------------------------------------------------------------------------------------------------------------
      ! to allocate flow variables
      !---------------------------------------------------------------------------------------------------------------------------------------------
      call Allocate_flow_variables (flow(l), domain(l))
      !---------------------------------------------------------------------------------------------------------------------------------------------
      ! to allocate thermal variables
      !---------------------------------------------------------------------------------------------------------------------------------------------
      if(domain(l)%ithermo == 1) call Allocate_thermo_variables (thermo(l), domain(l))
      !---------------------------------------------------------------------------------------------------------------------------------------------
      ! to set up Re, Fr etc 
      !---------------------------------------------------------------------------------------------------------------------------------------------
      call Update_Re(flow(l)%nrsttckpt, flow(l))
      if(domain(l)%ithermo == 1) &
      call Update_PrGr(flow(l), thermo(l))
      !---------------------------------------------------------------------------------------------------------------------------------------------
      ! to intialize primary variable
      !---------------------------------------------------------------------------------------------------------------------------------------------
      if (flow(l)%irestart == INITIAL_RANDOM) then

        call Initialize_flow_variables ( flow(l), domain(l) )
        if(domain(l)%ithermo == 1) &
        call Initialize_thermo_variables ( flow(l), thermo(l), domain(l) )

      else if (flow(l)%irestart == INITIAL_RESTART) then
        
        call read_instantanous_flow_raw_data(flow(l), domain(l))
        call restore_flow_variables_from_restart(flow(l), domain(l))
        if(domain(l)%ithermo == 1) then
            call read_instantanous_thermo_raw_data  (thermo(l), domain(l) )
            call restore_flow_variables_from_restart(flow(l),   domain(l))
        end if

      else if (flow(l)%irestart == INITIAL_INTERPL) then

        ! to add ...

      else
        call Print_error_msg("Error in flow initialisation flag.")
      end if

      call Test_algorithms()

!---------------------------------------------------------------------------------------------------------------------------------------------
! update initial results
!---------------------------------------------------------------------------------------------------------------------------------------------
      call Solve_momentum_eq(flow(l), domain(l), 0) ! check, necessary?
!---------------------------------------------------------------------------------------------------------------------------------------------
! to write out data for check
!---------------------------------------------------------------------------------------------------------------------------------------------
      call write_instantanous_flow_data(flow(l), domain(l))
      if(domain(l)%ithermo == 1) &
      call write_instantanous_thermo_data(thermo(l), domain(l)) 

    end do
    stop
!---------------------------------------------------------------------------------------------------------------------------------------------
! to test algorithms based on given values.
!---------------------------------------------------------------------------------------------------------------------------------------------
    
!---------------------------------------------------------------------------------------------------------------------------------------------
! test
!---------------------------------------------------------------------------------------------------------------------------------------------
    !if(nrank == 0) call Print_debug_end_msg
     
    return
  end subroutine Initialize_flow_thermal_fields
!=============================================================================================================================================
!> \brief Allocate flow and thermal variables.     
!--------------------------------------------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   private
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                         
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[in]     none          NA
!> \param[out]    none          NA
!=============================================================================================================================================
  subroutine Allocate_flow_variables (fl, dm)
    use parameters_constant_mod
    use mpi_mod
    implicit none

    type(t_domain), intent(in)    :: dm
    type(t_flow),   intent(inout) :: fl

    if(nrank == 0) call Print_debug_start_msg("Allocating flow variables ...")
    !---------------------------------------------------------------------------------------------------------------------------------------------
    ! default : x pencil. 
    ! varaible index is LOCAL. means 1:xsize(1)
    !---------------------------------------------------------------------------------------------------------------------------------------------
    call alloc_x(fl%qx,      dm%dpcc) ; fl%qx = ZERO
    call alloc_x(fl%qy,      dm%dcpc) ; fl%qy = ZERO
    call alloc_x(fl%qz,      dm%dccp) ; fl%qz = ZERO

    call alloc_x(fl%pres,    dm%dccc) ; fl%pres = ZERO
    call alloc_x(fl%pcor,    dm%dccc) ; fl%pcor = ZERO
    call alloc_x(fl%pcor_zpencil_ggg,    dm%dccc, .true.) ; fl%pcor_zpencil_ggg = ZERO

    call alloc_x(fl%mx_rhs,  dm%dpcc) ; fl%mx_rhs = ZERO
    call alloc_x(fl%my_rhs,  dm%dcpc) ; fl%my_rhs = ZERO
    call alloc_x(fl%mz_rhs,  dm%dccp) ; fl%mz_rhs = ZERO

    call alloc_x(fl%mx_rhs0, dm%dpcc) ; fl%mx_rhs0 = ZERO
    call alloc_x(fl%my_rhs0, dm%dcpc) ; fl%my_rhs0 = ZERO
    call alloc_x(fl%mz_rhs0, dm%dccp) ; fl%mz_rhs0 = ZERO

    if(dm%ithermo == 1) then
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
  !=============================================================================================================================================
  subroutine Allocate_thermo_variables (tm, dm)
    use parameters_constant_mod
    use mpi_mod
    use udf_type_mod
    use thermo_info_mod
    implicit none

    type(t_domain), intent(in)    :: dm
    type(t_thermo), intent(inout) :: tm

    if(dm%ithermo == 0) return
    if(nrank == 0) call Print_debug_start_msg("Allocating thermal variables ...")
    !---------------------------------------------------------------------------------------------------------------------------------------------
    ! default : x pencil. 
    ! varaible index is LOCAL. means 1:xsize(1)
    !---------------------------------------------------------------------------------------------------------------------------------------------
    call alloc_x(tm%dh,    dm%dccc) ; tm%dh    = ZERO
    call alloc_x(tm%hEnth, dm%dccc) ; tm%hEnth = ZERO
    call alloc_x(tm%kCond, dm%dccc) ; tm%kCond = ONE
    call alloc_x(tm%tTemp, dm%dccc) ; tm%tTemp = ONE

    if(nrank == 0) call Print_debug_end_msg
    return

  end subroutine Allocate_thermo_variables
!=============================================================================================================================================
!> \brief The main code for initializing flow variables
!--------------------------------------------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   private
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[in]     dm
!=============================================================================================================================================
  subroutine Initialize_flow_variables(fl, dm)
    use udf_type_mod
    use solver_tools_mod
    use parameters_constant_mod
    use burgers_eq_mod
    use visulisation_mod
#ifdef DEBUG
    use typeconvert_mod
#endif
    implicit none
    type(t_domain), intent(inout) :: dm
    type(t_flow),   intent(inout) :: fl
#ifdef DEBUG
    integer :: i, j, k, jj
    type(DECOMP_INFO) :: dtmp
#endif

    if(nrank == 0) call Print_debug_start_msg("Initializing flow variables ...")

    !---------------------------------------------------------------------------------------------------------------------------------------------
    ! to initialize flow velocity and pressure
    !---------------------------------------------------------------------------------------------------------------------------------------------
    if ( (dm%icase == ICASE_CHANNEL) .or. &
         (dm%icase == ICASE_PIPE) .or. &
         (dm%icase == ICASE_ANNUAL) ) then
      call Initialize_poiseuille_flow    (dm, fl%qx, fl%qy, fl%qz, fl%pres, fl%initNoise)
    else if (dm%icase == ICASE_TGV2D) then
      call Initialize_vortexgreen_2dflow (dm, fl%qx, fl%qy, fl%qz, fl%pres)
    else if (dm%icase == ICASE_TGV3D) then
      call Initialize_vortexgreen_3dflow (dm, fl%qx, fl%qy, fl%qz, fl%pres)
    else if (dm%icase == ICASE_BURGERS) then
      call Initialize_burgers_flow      (dm, fl%qx, fl%qy, fl%qz, fl%pres)
    else 
      if(nrank == 0) call Print_error_msg("No such case defined" )
    end if
    !---------------------------------------------------------------------------------------------------------------------------------------------
    ! to initialize pressure correction term
    !---------------------------------------------------------------------------------------------------------------------------------------------
    fl%pcor(:, :, :) = ZERO
#ifdef DEBUG
    call view_data_in_rank(fl%qx,   dm%dpcc, dm, 'ux', 0)
    call view_data_in_rank(fl%qy,   dm%dcpc, dm, 'uy', 0)
    call view_data_in_rank(fl%qz,   dm%dccp, dm, 'uz', 0)
    call view_data_in_rank(fl%pres, dm%dccc, dm, 'pr', 0)
#endif

#ifdef DEBUG
    dtmp = dm%dccc
    k = 2
    i = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      open(121, file = 'debugy_init_uvwp_'//trim(int2str(nrank))//'.dat', position="append")
      do j = 1, dm%dccc%xsz(2) 
        jj = dm%dpcc%xst(2) + j - 1
        write(121, *) jj, fl%qx(i+1, j, k), fl%qy(i, j, k), fl%qz(i, j, k+1), fl%pres(i, j, k)
      end do
    end if

    k = 2
    j = 2
    if( k >= dtmp%xst(3) .and. k <= dtmp%xen(3)) then
      if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
        open(221, file = 'debugx_init_uvwp_'//trim(int2str(nrank))//'.dat', position="append")
        do i = 1, dm%dccc%xsz(1)
          write(221, *) i, fl%qx(i, j, k), fl%qy(i, j+1, k), fl%qz(i, j, k+1), fl%pres(i, j, k)
        end do
      end if
    end if

    i = 2
    j = 2
    if( j >= dtmp%xst(2) .and. j <= dtmp%xen(2)) then
      open(321, file = 'debugz_init_uvwp_'//trim(int2str(nrank))//'.dat', position="append")
      do k = 1, dm%dccc%xsz(3)
        write(321, *) k, fl%qx(i+1, j, k), fl%qy(i, j+1, k), fl%qz(i, j, k), fl%pres(i, j, k)
      end do
    end if

#endif

    !---------------------------------------------------------------------------------------------------------------------------------------------
    ! to check maximum velocity
    !---------------------------------------------------------------------------------------------------------------------------------------------
    call Find_maximum_absvar3d(fl%qx, "maximum ux:")
    call Find_maximum_absvar3d(fl%qy, "maximum uy:")
    call Find_maximum_absvar3d(fl%qz, "maximum uz:")
    !call Check_mass_conservation(fl, dm) 
    !---------------------------------------------------------------------------------------------------------------------------------------------
    ! to set up flow iterations 
    !---------------------------------------------------------------------------------------------------------------------------------------------
    fl%time = ZERO
    fl%iteration = 0

    !if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
  !=============================================================================================================================================
  subroutine Initialize_thermo_variables( fl, tm, dm )
    use udf_type_mod
    use solver_tools_mod
    use thermo_info_mod
    use boundary_conditions_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm

    if (dm%ithermo /= 1) return
    !---------------------------------------------------------------------------------------------------------------------------------------------
    ! to initialize thermal variables 
    !---------------------------------------------------------------------------------------------------------------------------------------------
    if(nrank == 0) call Print_debug_mid_msg("Initializing thermal field ...")
    call Initialize_thermal_properties (fl, tm)

    call Calculate_massflux_from_velocity (fl, dm)
    !---------------------------------------------------------------------------------------------------------------------------------------------
    ! to set up old arrays 
    !---------------------------------------------------------------------------------------------------------------------------------------------
    fl%dDensm1(:, :, :) = fl%dDens(:, :, :)
    fl%dDensm2(:, :, :) = fl%dDens(:, :, :)
    
    tm%time = ZERO
    tm%iteration = 0
    
    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine

!=============================================================================================================================================
!> \brief Generate a flow profile for Poiseuille flow in channel or pipe.     
!--------------------------------------------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   private
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[in]     d             domain
!> \param[out]    ux_1c1          u(yc), velocity profile along wall-normal direction
!=============================================================================================================================================
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
      a = MAXP
      b = ZERO
      c = ONE
    end if

    do j = 1, dm%nc(2)
      yy = dm%yc(j)
      ux_1c1(j) = ( ONE - ( (yy - b)**2 ) / a / a ) * c
    end do

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine Generate_poiseuille_flow_profile
!=============================================================================================================================================
!> \brief Generate a flow profile for Poiseuille flow in channel or pipe.     
!--------------------------------------------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   private
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[in]     
!> \param[out]    
!=============================================================================================================================================
  subroutine Generate_random_field(dm, ux, uy, uz, lnoise)
    use random_number_generation_mod
    use parameters_constant_mod
    use mpi_mod
    implicit none
    type(t_domain), intent(in) :: dm
    real(WP),       intent(in) :: lnoise
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent(inout) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent(inout) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent(inout) :: uz
    
    integer :: seed
    integer :: i, j, k! local id
    integer :: n, nsz  
    integer :: ii, jj, kk ! global id
    real(WP) :: rd
    type(DECOMP_INFO) :: dtmp

    if(nrank == 0) call Print_debug_start_msg("Generating random field ...")
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   Initialisation in x pencil
    !---------------------------------------------------------------------------------------------------------------------------------------------
    seed = 0
    ux(:, :, :) = ZERO
    uy(:, :, :) = ZERO
    uz(:, :, :) = ZERO
    nsz = dm%np(1) * dm%np(2) * dm%np(3)

    do n = 1, NVD
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
          do k = 1, dtmp%xsz(3)
            kk = dtmp%xst(3) + k - 1
            seed = ii + jj + kk + nsz * (n - 1)
            call Initialize_random_number ( seed )
            call Generate_r_random( -ONE, ONE, rd)

            if(n == 1) ux(i, j, k) = lnoise * rd
            if(n == 2) uy(i, j, k) = lnoise * rd
            if(n == 3) uz(i, j, k) = lnoise * rd
          end do
        end do
      end do
    end do

    if(nrank==0) call Print_debug_end_msg
    return
  end subroutine
!=============================================================================================================================================
!> \brief Initialize Poiseuille flow in channel or pipe.     
!--------------------------------------------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   private
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[in]     d             domain
!> \param[out]    f             flow
!=============================================================================================================================================
  subroutine Initialize_poiseuille_flow(dm, ux, uy, uz, p, lnoise)
    use input_general_mod
    use udf_type_mod
    use boundary_conditions_mod
    use parameters_constant_mod
    use wtformat_mod
    implicit none
    type(t_domain), intent(in ) :: dm
    real(WP),       intent(in ) :: lnoise   
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent(out) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent(out) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent(out) :: uz
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)), intent(out) ::  p
    integer :: pf_unit
    integer :: i, j, k, jj
    real(WP) :: ubulk
    real(WP) :: ux_1c1(dm%nc(2))
    real(WP) :: uxxza (dm%nc(2))
    real(WP) :: uyxza (dm%np(2))
    real(WP) :: uzxza (dm%nc(2))
    real(WP) :: ux_ypencil(dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3))

    type(DECOMP_INFO) :: dtmp

    if(nrank == 0) call Print_debug_mid_msg("Initializing Poiseuille flow field ...")
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   x-pencil : initial
    !---------------------------------------------------------------------------------------------------------------------------------------------
    p (:, :, :) = ZERO
    ux(:, :, :) = ZERO
    uy(:, :, :) = ZERO
    uz(:, :, :) = ZERO
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   x-pencil : to get random fields [-1,1] for ux, uy, uz
    !---------------------------------------------------------------------------------------------------------------------------------------------
    call Generate_random_field(dm, ux, uy, uz, lnoise)
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   x-pencil : Get the averaged u, v, and w in both the x and z directions
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !if(nrank == 0) call Print_debug_mid_msg("Get the u, v, w, averged in x and z directions ...")
    uxxza = ZERO
    uyxza = ZERO
    uzxza = ZERO
    call Calculate_xz_mean_yprofile(ux, dm%dpcc, dm%nc(2), uxxza)
    call Calculate_xz_mean_yprofile(uy, dm%dcpc, dm%np(2), uyxza)
    call Calculate_xz_mean_yprofile(uz, dm%dpcc, dm%nc(2), uzxza)
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   x-pencil : Ensure u-u_given, v, w, averaged in x and z direction is zero.
    !              added perturbation is zero in mean. 
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !if(nrank == 0) call Print_debug_mid_msg("Calculate xzmean perturbation...")
    call Adjust_to_xzmean_zero(ux, dm%dpcc, dm%nc(1), uxxza)
    call Adjust_to_xzmean_zero(uy, dm%dcpc, dm%np(2), uyxza)
    call Adjust_to_xzmean_zero(uz, dm%dccp, dm%nc(3), uzxza)
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   x-pencil : to get Poiseuille profile for all ranks
    !---------------------------------------------------------------------------------------------------------------------------------------------
    ux_1c1(:) = ZERO
    call Generate_poiseuille_flow_profile (dm, ux_1c1)
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   x-pencil : to add profile to ux (default: x streamwise)
    !---------------------------------------------------------------------------------------------------------------------------------------------
    dtmp = dm%dpcc
    do i = 1, dtmp%xsz(1)
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        do k = 1, dtmp%xsz(3)
          ux(i, j, k) =  ux(i, j, k) + ux_1c1(jj)
        end do
      end do
    end do
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   x-pencil : build up boundary
    !---------------------------------------------------------------------------------------------------------------------------------------------
    call Apply_BC_velocity(dm, ux, uy, uz)
    call Find_maximum_absvar3d(ux, "maximum ux:")
    call Find_maximum_absvar3d(uy, "maximum uy:")
    call Find_maximum_absvar3d(uz, "maximum uz:")
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   x-pencil : Ensure the mass flow rate is 1.
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !if(nrank == 0) call Print_debug_mid_msg("Ensure u, v, w, averaged in x and z direction is zero...")
    call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), dm%fbcy(:, 1), dm, dm%dpcc, ux, ubulk)
    if(nrank == 0) then
      Call Print_debug_mid_msg("  The initial mass flux is:")
      write (*, wrtfmt1r) ' average[u(x,y,z)]_[x,y,z]: ', ubulk
    end if

    ux(:, :, :) = ux(:, :, :) / ubulk

    call Apply_BC_velocity(dm, ux, uy, uz)
    call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), dm%fbcy(:, 1), dm, dm%dpcc, ux, ubulk)
    if(nrank == 0) then
      Call Print_debug_mid_msg("  The scaled mass flux is:")
      write (*, wrtfmt1r) ' average[u(x,y,z)]_[x,y,z]: ', ubulk
    end if
    ! to do : to add a scaling for turbulence generator inlet scaling, u = u * m / rho

    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   some checking
    !---------------------------------------------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(ux, ux_ypencil, dm%dpcc)
    !---------------------------------------------------------------------------------------------------------------------------------------------
    !   Y-pencil : write out velocity profile
    !---------------------------------------------------------------------------------------------------------------------------------------------
    if(nrank == 0) then
      open ( newunit = pf_unit,     &
              file    = 'check_poiseuille_profile.dat', &
              status  = 'replace',         &
              action  = 'write')
      write(pf_unit, '(A)') "# yc, ux_laminar, ux_real"
      do j = 1, dm%nc(2)
        write(pf_unit, '(5ES13.5)') dm%yc(j), ux_1c1(j), ux_ypencil(dm%dpcc%yen(1)/2, j, dm%dpcc%yen(3)/2)
      end do
      close(pf_unit)
    end if
    
    !if(nrank == 0) call Print_debug_end_msg

    return
  end subroutine  Initialize_poiseuille_flow
!=============================================================================================================================================
!=============================================================================================================================================
!> \brief Initialize Vortex Green flow
!>
!> This subroutine is called locally once.
!>
!---------------------------------------------------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------------------------
!   ux in x-pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------------------------
!   uy in x-pencil
!--------------------------------------------------------------------------------------------------------------------------------------------- 
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
!---------------------------------------------------------------------------------------------------------------------------------------------
!   uz in x-pencil
!--------------------------------------------------------------------------------------------------------------------------------------------- 
    uz(:, :, :) =  ZERO
!---------------------------------------------------------------------------------------------------------------------------------------------
!   p in x-pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
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
!=============================================================================================================================================
!=============================================================================================================================================
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

!---------------------------------------------------------------------------------------------------------------------------------------------
!   X-pencil : Find Max. error of ux
!---------------------------------------------------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------------------------
!   X-pencil : Find Max. error of uy
!---------------------------------------------------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------------------------
!   X-pencil : Find Max. error of p
!---------------------------------------------------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------------------------
!   X-pencil : write data in rank=0
!---------------------------------------------------------------------------------------------------------------------------------------------
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
!=============================================================================================================================================
!=============================================================================================================================================
!> \brief Initialize Vortex Green flow
!>
!> This subroutine is called locally once.
!>
!---------------------------------------------------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------------------------
!   ux in x-pencil
!--------------------------------------------------------------------------------------------------------------------------------------------- 
    dtmp = dm%dpcc
    do k = 1, dtmp%xsz(3)
      kk = dtmp%xst(3) + k - 1
      zc = dm%h(3) * (real(kk - 1, WP) + HALF) - PI
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yc = dm%yc(jj)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xp = dm%h(1) * real(ii - 1, WP) - PI
          ux(i, j, k) =  sin_wp ( xp ) * cos_wp ( yc ) * cos_wp ( zc )
          ! if(i<=dtmp%xsz(1)/2 + 1) then
          !   ux(i, j, k) =  xp
          ! else
          !   ux(i, j, k) = -xp + twopi
          !end if
        end do
      end do
    end do
!---------------------------------------------------------------------------------------------------------------------------------------------
!   uy in x-pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
    dtmp = dm%dcpc
    do k = 1, dtmp%xsz(3)
      kk = dtmp%xst(3) + k - 1
      zc = dm%h(3) * (real(kk - 1, WP) + HALF) - PI
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yp = dm%yp(jj)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xc = dm%h(1) * (real(ii - 1, WP) + HALF) - PI
          !uy(i, j, k) = sin_wp(yp)
          uy(i, j, k) = -cos_wp ( xc ) * sin_wp ( yp ) * cos_wp ( zc )
          ! if(jj<=dtmp%ysz(2)/2 + 1) then
          !   uy(i, j, k) =  yp
          ! else
          !   uy(i, j, k) = -yp + twopi
          ! end if
        end do
      end do
    end do
!---------------------------------------------------------------------------------------------------------------------------------------------
!   uz in x-pencil
!--------------------------------------------------------------------------------------------------------------------------------------------- 
    !uz(:, :, :) =  ZERO
    dtmp = dm%dccp
    do k = 1, dtmp%xsz(3)
      kk = dtmp%xst(3) + k - 1
      zp = dm%h(3) * real(kk - 1, WP) - PI
      do j = 1, dtmp%xsz(2)
        yc = dm%yc(jj)
        do i = 1, dtmp%xsz(1)
          xc = dm%h(1) * (real(ii - 1, WP) + HALF) - PI
          !uz(i, j, k) = -cos_wp ( xc ) *  cos_wp ( yc ) * sin_wp ( zp )
        end do
      end do
    end do
!---------------------------------------------------------------------------------------------------------------------------------------------
!   p in x-pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
    dtmp = dm%dccc
    do k = 1, dtmp%xsz(3)
      kk = dtmp%xst(3) + k - 1
      zc = dm%h(3) * (real(kk - 1, WP) + HALF) - PI
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        yc = dm%yc(jj)
        do i = 1, dtmp%xsz(1)
          ii = dtmp%xst(1) + i - 1
          xc = dm%h(1) * (real(ii - 1, WP) + HALF) - PI
          p(i, j, k)= ONE / SIXTEEN * ( cos(TWO * xc) + cos(TWO * yc) ) * &
                      (cos(TWO * zc) + TWO)
        end do
      end do
    end do

    if(nrank == 0) call Print_debug_end_msg
    
    return
  end subroutine Initialize_vortexgreen_3dflow

end module flow_thermo_initialiasation