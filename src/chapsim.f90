!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>                      CHAPSim version 2.0.0
!>                      --------------------------
!> This file is part of CHAPSim, a general-purpose CFD tool.
!>
!> This program is free software; you can redistribute it and/or modify it under
!> the terms of the GNU General Public License as published by the Free Software
!> Foundation; either version 3 of the License, or (at your option) any later
!> version.
!>
!> This program is distributed in the hope that it will be useful, but WITHOUT
!> ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
!> FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
!> details.
!>
!> You should have received a copy of the GNU General Public License along with
!> this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
!> Street, Fifth Floor, Boston, MA 02110-1301, USA.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================================================================================================
!> \file chapsim.f90
!> \brief the main program.
!> \author Wei Wang, wei.wang@stfc.ac.uk
!> \date 
!==========================================================================================================
program chapsim
  implicit none

  call Initialize_chapsim
  call Solve_eqs_iteration
  call Finalise_chapsim
  
end program
!==========================================================================================================
!> \brief Initialisation and preprocessing of geometry, mesh and tools
!> This subroutine is called at beginning of the main program
!==========================================================================================================
subroutine Initialize_chapsim
  use boundary_conditions_mod
  use code_performance_mod
  use continuity_eq_mod
  use decomp_2d_poisson
  use domain_decomposition_mod
  use files_io_mod
  use flow_thermo_initialiasation
  use geometry_mod
  use input_general_mod
  use io_monitor_mod
  use io_tools_mod
  use mpi_mod
  use operations
  use poisson_interface_mod
  use solver_tools_mod
  use statistics_mod
  use thermo_info_mod
  implicit none
  integer :: i

!----------------------------------------------------------------------------------------------------------
! initialisation of mpi, nrank, nproc
!----------------------------------------------------------------------------------------------------------
  call create_directory
  call call_cpu_time(CPU_TIME_CODE_START, 0, 0)
  call Initialize_mpi
!----------------------------------------------------------------------------------------------------------
! reading input parameters
!----------------------------------------------------------------------------------------------------------
  call Read_input_parameters
!----------------------------------------------------------------------------------------------------------
! build up geometry information
!----------------------------------------------------------------------------------------------------------
  do i = 1, nxdomain
    call configure_bc_type(domain(i)) 
    call Buildup_geometry_mesh_info(domain(i))
  end do
!----------------------------------------------------------------------------------------------------------
! build up domain decomposition
!----------------------------------------------------------------------------------------------------------
  do i = 1, nxdomain
    call Buildup_mpi_domain_decomposition(domain(i)) 
  end do
!----------------------------------------------------------------------------------------------------------
! build up thermo_mapping_relations, independent of any domains
!----------------------------------------------------------------------------------------------------------
  do i = 1, nxdomain
    if(domain(i)%is_thermo) then
      call Buildup_thermo_mapping_relations(thermo(i))
      cycle
    end if
  end do
!----------------------------------------------------------------------------------------------------------
! build up initial boundary values
!----------------------------------------------------------------------------------------------------------
  do i = 1, nxdomain
    if(domain(i)%is_thermo) call buildup_thermo_bc_geo(domain(i), thermo(i)) ! should be before flow.
    call buildup_flow_bc_geo(domain(i), flow(i)) 
    if(domain(i)%is_thermo) call apply_gxgygz_bc_geo(domain(i), flow(i), thermo(i))
  end do
!----------------------------------------------------------------------------------------------------------
! build up output_io
!----------------------------------------------------------------------------------------------------------
  do i = 1, nxdomain
    call initialize_decomp_io(domain(i))
    call write_monitor_ini(domain(i))
    call write_snapshot_ini(domain(i))
    call init_statistics_flow(domain(i), flow(i))
    if(domain(i)%is_thermo) call init_statistics_thermo(domain(i), thermo(i))
    call mesh_output(domain(i))
  end do
!----------------------------------------------------------------------------------------------------------
! build up operation coefficients for all x-subdomains
!----------------------------------------------------------------------------------------------------------
  call Prepare_LHS_coeffs_for_operations

#ifdef DEBUG_TEST
  call Test_algorithms()
#endif
!----------------------------------------------------------------------------------------------------------
! build up fft basic info
!----------------------------------------------------------------------------------------------------------
  do i = 1, nxdomain
    call build_up_poisson_interface(domain(i))
    if(nrank == 0 ) call Print_debug_start_msg("Initializing Poisson solver ...")
    call decomp_2d_poisson_init()
    if(nrank == 0 ) call Print_debug_end_msg
  end do
!----------------------------------------------------------------------------------------------------------
! Initialize flow and thermo fields
!----------------------------------------------------------------------------------------------------------
  do i = 1, nxdomain
    call Initialize_flow_fields(domain(i), flow(i))
    if(domain(i)%is_thermo) then
      call Initialize_thermo_fields(domain(i), flow(i), thermo(i))
      call Find_maximum_velocity(domain(i), flow(i)%gx, flow(i)%gy, flow(i)%gz)
      call Check_mass_conservation(domain(i), flow(i), 'initialization', thermo(i))
    else
      call Find_maximum_velocity(domain(i), flow(i)%qx, flow(i)%qy, flow(i)%qz)
      call Check_mass_conservation(domain(i), flow(i), 'initialization') 
    end if    
    call write_snapshot_flow(domain(i), flow(i))
  end do
!----------------------------------------------------------------------------------------------------------
! update interface values for multiple domain
!----------------------------------------------------------------------------------------------------------
  do i = 1, nxdomain - 1
    call update_flow_bc_2dm_halo(domain(i), flow(i), domain(i+1), flow(i+1))
    if(domain(i)%is_thermo) then
      call update_thermo_bc_2dm_halo(domain(i), thermo(i), domain(i+1), thermo(i+1))
    end if
  end do
  do i = 1, nxdomain
    if(domain(i)%is_thermo) then
      call apply_gxgygz_bc_geo(domain(i), flow(i), thermo(i))
    end if
  end do
  
  return
end subroutine Initialize_chapsim

!==========================================================================================================
!==========================================================================================================
!> \brief solve the governing equations in iteration
!>
!> This subroutine is the main solver. 
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
subroutine Solve_eqs_iteration
  !use iso_fortran_env
  use solver_tools_mod!,   only : Check_cfl_diffusion, Check_cfl_convection
  use continuity_eq_mod
  !use poisson_mod
  use eq_energy_mod
  use eq_momentum_mod
  use flow_thermo_initialiasation 
  use code_performance_mod
  use thermo_info_mod
  use vars_df_mod
  use input_general_mod
  use mpi_mod
  use wtformat_mod
  use io_monitor_mod
  use io_tools_mod
  use io_restart_mod
  use statistics_mod
  use typeconvert_mod
  use boundary_conditions_mod
  implicit none

  logical, allocatable :: is_flow  (:)
  logical, allocatable :: is_thermo(:)
  integer :: i
  integer :: iter, isub
  integer :: iteration
  integer :: niter

  
  !==========================================================================================================
  ! flow advancing/marching iteration/time control
  !==========================================================================================================
  iteration = HUGE(0)
  niter     = 0
  do i = 1, nxdomain
     if( flow(i)%iteration    < iteration) iteration = flow(i)%iteration
     if( flow(i)%nIterFlowEnd > niter)     niter     = flow(i)%nIterFlowEnd
     if( domain(i)%is_thermo) then
       if (thermo(i)%iteration      < iteration) iteration = thermo(i)%iteration
       if (thermo(i)%nIterThermoEnd > niter)     niter     = thermo(i)%nIterThermoEnd
     end if
  end do

  allocate(is_flow  (nxdomain)) 
  allocate(is_thermo(nxdomain)) 
  is_flow   = .false.
  is_thermo = .false.


  call call_cpu_time(CPU_TIME_STEP_START, iteration, niter)

  if(nrank == 0) call Print_debug_start_msg("Solving the governing equations ...")

  do iter = iteration + 1, niter
    if( mod(iter, cpu_nfre) == 0) call call_cpu_time(CPU_TIME_ITER_START, iteration, niter, iter)
    
    !==========================================================================================================
    ! Solver preparation for each domain
    !==========================================================================================================
    do i = 1, nxdomain
      !----------------------------------------------------------------------------------------------------------
      !      setting up 1/re, 1/re/prt, gravity, etc
      !----------------------------------------------------------------------------------------------------------
      call Update_Re(iter, flow(i))
      if(domain(i)%is_thermo) &
      call Update_PrGr(flow(i), thermo(i))
      !----------------------------------------------------------------------------------------------------------
      !      setting up flow solver
      !----------------------------------------------------------------------------------------------------------
      if ( (iter >= flow(i)%nIterFlowStart) .and. (iter <=flow(i)%nIterFlowEnd)) then
        is_flow(i) = .true.
        if (nrank == 0) write(*, wrtfmt1r) "flow field physical time (s) = ", flow(i)%time
        flow(i)%time = flow(i)%time + domain(i)%dt
        flow(i)%iteration = flow(i)%iteration + 1
        call Check_cfl_diffusion (domain(i)%h2r(:), domain(i)%rci, flow(i)%rre, domain(i)%dt)
        call Check_cfl_convection(domain(i), flow(i))
      end if
      !----------------------------------------------------------------------------------------------------------
      !     setting up thermo solver
      !----------------------------------------------------------------------------------------------------------
      if(domain(i)%is_thermo) then
        if ( (iter >= thermo(i)%nIterThermoStart) .and. (iter <= thermo(i)%nIterThermoEnd)) then
          is_thermo(i) = .true.
          if (nrank == 0) write(*, wrtfmt1r) "thermal field physical time = ", thermo(i)%time
          thermo(i)%time = thermo(i)%time  + domain(i)%dt
          thermo(i)%iteration = thermo(i)%iteration + 1
        end if
      end if
    end do
!==========================================================================================================
!  main solver, domain coupling in each sub-iteration (check)
!==========================================================================================================
    do isub = 1, domain(1)%nsubitr
!----------------------------------------------------------------------------------------------------------
!     solve governing equations and update b.c. if necessary
!----------------------------------------------------------------------------------------------------------
      do i = 1, nxdomain
        if(is_thermo(i)) then
          call Solve_energy_eq  (domain(i), flow(i), thermo(i), isub)
          call update_thermo_bc_1dm_halo(domain(i), thermo(i))
        end if
        if(is_flow(i)) then
          if(is_thermo(i)) then
            call Solve_momentum_eq(domain(i), flow(i), isub, thermo(i))
          else
            call Solve_momentum_eq(domain(i), flow(i), isub)
          end if
          call update_flow_bc_1dm_halo(domain(i), flow(i))
        end if
      end do
!----------------------------------------------------------------------------------------------------------
! update interface values for multiple domain
!----------------------------------------------------------------------------------------------------------
      do i = 1, nxdomain - 1
        if(is_thermo(i)) call update_thermo_bc_2dm_halo(domain(i), thermo(i), domain(i+1), thermo(i+1))
        if(is_flow(i)) then
          if(is_thermo(i+1)) then
            call update_flow_bc_2dm_halo(domain(i), flow(i), domain(i+1), flow(i+1), thermo(i+1))
          else
            call update_flow_bc_2dm_halo(domain(i), flow(i), domain(i+1), flow(i+1))
          end if
        end if
      end do
    end do

    !==========================================================================================================
    ! result post-processing for each domain
    !==========================================================================================================
    do i = 1, nxdomain
      !----------------------------------------------------------------------------------------------------------
      !  validation for each time step
      !----------------------------------------------------------------------------------------------------------
      !if(nrank == 0) call Print_debug_mid_msg("For domain id = "//trim(int2str(i)))
      if(is_flow(i)) then
        call Find_maximum_velocity(domain(i), flow(i)%qx, flow(i)%qy, flow(i)%qz)
        call Check_mass_conservation(domain(i), flow(i), 'iteration = '//trim(int2str(iter))) 
      end if

      !----------------------------------------------------------------------------------------------------------
      !  update statistics
      !----------------------------------------------------------------------------------------------------------
      if (iter > domain(i)%stat_istart .and. is_flow(i)) then
        call update_statistics_flow(domain(i), flow(i))
      end if
      if(domain(i)%is_thermo .and. is_thermo(i)) then
        if (iter > domain(i)%stat_istart) then
          call update_statistics_thermo(domain(i), thermo(i))
        end if
      end if
      !----------------------------------------------------------------------------------------------------------
      !  monitoring 
      !----------------------------------------------------------------------------------------------------------
      !if(domain(i)%icase == ICASE_TGV2D) call Validate_TGV2D_error (domain(i), flow(i))
      if(is_flow(i)) then
        call write_monitor_total(domain(i), flow(i))
        call write_monitor_flow(domain(i), flow(i))
      end if
      if(domain(i)%is_thermo .and. is_thermo(i)) then
        call write_monitor_thermo(domain(i), thermo(i))
      end if
      !----------------------------------------------------------------------------------------------------------
      !  write out check point data for restart
      !----------------------------------------------------------------------------------------------------------
      if (mod(iter, domain(i)%ckpt_nfre) == 0) then
        if(is_flow(i)) then
          call write_instantanous_flow_raw_data(domain(i), flow(i))
          if(iter > domain(i)%stat_istart) call write_statistics_flow(domain(i), flow(i))
        end if
        if(domain(i)%is_thermo .and. is_thermo(i)) then
          call write_instantanous_thermo_raw_data(domain(i), thermo(i))
          if(iter > domain(i)%stat_istart) call write_statistics_thermo(domain(i), thermo(i))
        end if
      end if
      !----------------------------------------------------------------------------------------------------------
      ! write data for visualisation
      !----------------------------------------------------------------------------------------------------------
      if(MOD(iter, domain(i)%visu_nfre) == 0) then
        if(is_flow(i)) call write_snapshot_flow(domain(i), flow(i))
        if(domain(i)%is_thermo .and. is_thermo(i)) then
          call write_snapshot_thermo(domain(i), thermo(i))
        end if
        if(iter > domain(i)%stat_istart ) then
          if(is_flow(i)) call write_snapshot_flow_stat(domain(i), flow(i))
          if(domain(i)%is_thermo .and. is_thermo(i)) then
            call write_snapshot_thermo_stat(domain(i), thermo(i))
          end if
        end if
      end if

    end do ! domain

    if( mod(iter, cpu_nfre) == 0) call call_cpu_time(CPU_TIME_ITER_END, iteration, niter, iter)

  end do ! iteration

  call call_cpu_time(CPU_TIME_STEP_END, iteration, niter)
  
  return
end subroutine Solve_eqs_iteration


subroutine Finalise_chapsim
  use mpi_mod
  use code_performance_mod
  implicit none

  call call_cpu_time(CPU_TIME_CODE_END, 0, 0)
  call Finalise_mpi()
  return
end subroutine
