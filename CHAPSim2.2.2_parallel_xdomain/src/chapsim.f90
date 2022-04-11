!-------------------------------------------------------------------------------
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

!-------------------------------------------------------------------------------
!===============================================================================
!> \file chapsim.f90
!>
!> \brief The main program.
!>
!===============================================================================
program chapsim
  implicit none

  call Initialize_chapsim
  call Solve_eqs_iteration
  
end program
!===============================================================================
!> \brief Initialisation and preprocessing of geometry, mesh and tools
!>
!> This subroutine is called at beginning of the main program
!>
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[in]     none          NA
!> \param[out]    none          NA
!===============================================================================
subroutine Initialize_chapsim
  use code_performance_mod
  use input_general_mod
  use geometry_mod
  use thermo_info_mod
  use operations
  use domain_decomposition_mod
  !use poisson_mod
  use flow_thermo_initialiasation
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
  if(is_any_energyeq) call Buildup_thermo_mapping_relations
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
  !call Prepare_poisson_fft
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
  implicit none

  logical :: is_flow   = .false.
  logical :: is_thermo = .false.
  integer :: i
  integer :: iter, isub
  integer :: nrsttckpt
  integer :: niter

  if(nrank == 0) call Print_debug_start_msg("Solving the governing equations ...")

  nrsttckpt = HUGE(0)
  niter     = 0
  do i = 1, nxdomain
     if( flow(i)%nrsttckpt < nrsttckpt) nrsttckpt = flow(i)%nrsttckpt
     if( flow(i)%nIterFlowEnd > niter)  niter     = flow(i)%nIterFlowEnd
     if( is_any_energyeq) then
       if (thermo(i)%nIterThermoEnd > niter) niter = thermo(i)%nIterThermoEnd
     end if
  end do

  do iter = nrsttckpt + 1, niter
    call Call_cpu_time(CPU_TIME_ITER_START, nrsttckpt, niter, iter)
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
        end if
      end if
!===============================================================================
!     main solver
!===============================================================================
      do isub = 1, domain(i)%nsubitr
        if(is_thermo) call Solve_energy_eq  (flow(i), thermo(i), domain(i), isub)
        if(is_flow)   call Solve_momentum_eq(flow(i), domain(i), isub)
!#ifdef DEBUG
        write (OUTPUT_UNIT, '(A, I1)') "  Sub-iteration in RK = ", isub
        call Check_mass_conservation(flow(i), domain(i)) 
        call Check_maximum_velocity(flow(i)%qx, flow(i)%qy, flow(i)%qz)
!#endif
      end do
!
    !comment this part code for testing 
    ! below is for validation
    ! cpu time will be calculated later today 
!===============================================================================
!     validation
!===============================================================================
      call Check_mass_conservation(flow(i), domain(i)) 
      if(domain(i)%icase == ICASE_TGV2D) call Validate_TGV2D_error (flow(i), domain(i))

      call Call_cpu_time(CPU_TIME_ITER_END, nrsttckpt, niter, iter)
!===============================================================================
!   visualisation
!===============================================================================
      if(MOD(iter, domain(i)%nvisu) == 0) then
        !call Display_vtk_slice(domain, 'xy', 'u', 1, flow(i)%qx, iter)
        !call Display_vtk_slice(domain, 'xy', 'v', 2, flow(i)%qy, iter)
        !call Display_vtk_slice(domain, 'xy', 'p', 0, flow(i)%pres, iter)
      end if
    end do
  end do


  call Call_cpu_time(CPU_TIME_CODE_END, nrsttckpt, niter)
  call Finalise_mpi()
  return
end subroutine Solve_eqs_iteration

