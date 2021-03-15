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

  call Initialize_chapsim ()
  

  call Test_schemes()

  call Finalise_chapsim ()
  
end program

!===============================================================================
!===============================================================================
!> \brief Initialisation and preprocessing of the flow solver
!>
!> This subroutine is called at beginning of the main program
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
subroutine Initialize_chapsim()
  !use mpi_mod
  use input_general_mod
  use input_thermo_mod
  !use domain_decomposition_mod
  use geometry_mod
  use flow_variables_mod
  use operations
  implicit none

  !call Initialize_mpi()
  call Initialize_general_input ()
  call Initialize_thermo_input ()
  call Initialize_geometry_variables ()
  call Prepare_coeffs_for_operations()

  !call Initialize_domain_decompsition ()
  call Initialize_flow_variables ()
  return
end subroutine Initialize_chapsim
!===============================================================================
!===============================================================================
!> \brief Finalising the flow solver
!>
!> This subroutine is called at the end of the main program
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
subroutine Finalise_chapsim()
  use mpi_mod
  implicit none

  !call Deallocate_all_variables
  !call Finalise_mpi()
  return
end subroutine Finalise_chapsim
!===============================================================================
!===============================================================================
!> \brief In-code independent test code for algorithms and schemes
!>
!> This subroutine is only called in the main program for testing.
!> Please select the test options which you are interested in.
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
subroutine Test_schemes()
  use geometry_mod
  use operations
  use tridiagonal_matrix_algorithm
  implicit none

  !call Test_TDMA_cyclic
  !call Test_TDMA_noncyclic
  call Test_interpolation(domain)
  call Test_1st_derivative(domain)
  return 
end subroutine 
!===============================================================================



