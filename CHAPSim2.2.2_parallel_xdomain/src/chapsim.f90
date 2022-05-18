!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================================
!> \file chapsim.f90
!> \brief the main program.
!> \author Wei Wang wei.wang@stfc.ac.uk
!> \date 
!===============================================================================
program chapsim
  implicit none

  call Initialize_chapsim
  call Solve_eqs_iteration
  
end program
!===============================================================================
!> \brief Initialisation and preprocessing of geometry, mesh and tools
!! This subroutine is called at beginning of the main program
!===============================================================================
subroutine Initialize_chapsim
  use code_performance_mod
  use input_general_mod
  use geometry_mod
  use thermo_info_mod
  use operations
  use domain_decomposition_mod
  use flow_thermo_initialiasation
  use mpi_mod
  use code_performance_mod
  use decomp_2d_poisson
  implicit none
  integer :: i

  !-------------------------------------------------------------------------------
  ! initialisation of mpi, nrank, nproc
  !-------------------------------------------------------------------------------
  call Initialize_mpi
  call Call_cpu_time(CPU_TIME_CODE_START, 0, 0)
  !-------------------------------------------------------------------------------
  ! reading input parameters
  !-------------------------------------------------------------------------------
  call Read_input_parameters
  !-------------------------------------------------------------------------------
  ! build up geometry information
  !-------------------------------------------------------------------------------
  do i = 1, nxdomain
    call Buildup_geometry_mesh_info(domain(i)) 
  end do
  !-------------------------------------------------------------------------------
  ! build up thermo_mapping_relations, independent of any domains
  !-------------------------------------------------------------------------------
  if(is_any_energyeq) then
    i = 1 
    call Buildup_thermo_mapping_relations(thermo(i), domain(i))
  end if
!-------------------------------------------------------------------------------
! build up operation coefficients for all x-subdomains
!-------------------------------------------------------------------------------
  call Prepare_LHS_coeffs_for_operations
!-------------------------------------------------------------------------------
! build up domain decomposition
!-------------------------------------------------------------------------------
  call Buildup_mpi_domain_decomposition
!-------------------------------------------------------------------------------
! build up fft basic info
!-------------------------------------------------------------------------------
  do i = 1, nxdomain
    call build_up_poisson_interface(domain(i))
    call decomp_2d_poisson_init()
  end do
!-------------------------------------------------------------------------------
! Initialize flow and thermo fields
!-------------------------------------------------------------------------------
  call Initialize_flow_thermal_fields

  return
end subroutine Initialize_chapsim

!===============================================================================
!===============================================================================
!> \brief solve the governing equations in iteration
!>
!> This subroutine is the main solver. 
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
subroutine Solve_eqs_iteration
  use iso_fortran_env
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
  implicit none

  logical :: is_flow   = .false.
  logical :: is_thermo = .false.
  integer :: i
  integer :: iter, isub
  integer :: iteration
  integer :: niter

  if(nrank == 0) call Print_debug_start_msg("Solving the governing equations ...")

  !===============================================================================
  ! flow advancing/marching iteration/time control
  !===============================================================================
  iteration = HUGE(0)
  niter     = 0
  do i = 1, nxdomain
     if( flow(i)%iteration    < iteration) iteration = flow(i)%iteration
     if( flow(i)%nIterFlowEnd > niter)     niter     = flow(i)%nIterFlowEnd
     if( is_any_energyeq) then
       if (thermo(i)%iteration      < iteration) iteration = thermo(i)%iteration
       if (thermo(i)%nIterThermoEnd > niter)     niter     = thermo(i)%nIterThermoEnd
     end if
  end do

  do iter = iteration + 1, niter
    call Call_cpu_time(CPU_TIME_ITER_START, iteration, niter, iter)
    !===============================================================================
    ! Solver for each domain
    !===============================================================================
    do i = 1, nxdomain
      !===============================================================================
      !      setting up 1/re, 1/re/prt, gravity, etc
      !===============================================================================
      call Update_Re(iter, flow(i))
      if(domain(i)%ithermo == 1) &
      call Update_PrGr(flow(i), thermo(i))
      !===============================================================================
      !      setting up flow solver
      !===============================================================================
      if ( (iter >= flow(i)%nIterFlowStart) .and. (iter <=flow(i)%nIterFlowEnd)) then
        is_flow = .true.
        flow(i)%time = flow(i)%time + domain(i)%dt
        flow(i)%iteration = flow(i)%iteration + 1
        call Check_cfl_diffusion (domain(i)%h2r(:), flow(i)%rre, domain(i)%dt)
        call Check_cfl_convection(flow(i)%qx, flow(i)%qy, flow(i)%qz, domain(i))
      end if
      !===============================================================================
      !     setting up thermo solver
      !===============================================================================
      if(domain(i)%ithermo == 1) then
        if ( (iter >= thermo(i)%nIterThermoStart) .and. (iter <= thermo(i)%nIterThermoEnd)) then
          is_thermo = .true.
          thermo(i)%time = thermo(i)%time  + domain(i)%dt
          thermo(i)%iteration = thermo(i)%iteration + 1
        end if
      end if
      !===============================================================================
      !     main solver
      !===============================================================================
      do isub = 1, domain(i)%nsubitr
        if(is_thermo) call Solve_energy_eq  (flow(i), thermo(i), domain(i), isub)
        if(is_flow)   call Solve_momentum_eq(flow(i), domain(i), isub)
#ifdef DEBUG
        write (OUTPUT_UNIT, '(A, I1)') "  Sub-iteration in RK = ", isub
        call Check_mass_conservation(flow(i), domain(i)) 
        call Check_maximum_velocity(flow(i)%qx, flow(i)%qy, flow(i)%qz)
#endif
      end do

      !comment this part code for testing 
      ! below is for validation
      ! cpu time will be calculated later today 
      !===============================================================================
      !     validation
      !===============================================================================
      call Check_mass_conservation(flow(i), domain(i)) 
      if(domain(i)%icase == ICASE_TGV2D) call Validate_TGV2D_error (flow(i), domain(i))

      call Call_cpu_time(CPU_TIME_ITER_END, iteration, niter, iter)
      !===============================================================================
      !   visualisation
      !===============================================================================
      if(MOD(iter, domain(i)%nvisu) == 0) then
        !call Display_vtk_slice(domain, 'xy', 'u', 1, flow(i)%qx, iter)
        !call Display_vtk_slice(domain, 'xy', 'v', 2, flow(i)%qy, iter)
        !call Display_vtk_slice(domain, 'xy', 'p', 0, flow(i)%pres, iter)
      end if
    end do ! domain

  end do ! iteration


  call Call_cpu_time(CPU_TIME_CODE_END, iteration, niter)
  call Finalise_mpi()
  return
end subroutine Solve_eqs_iteration

