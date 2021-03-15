program chapsim
  implicit none

  call Initialize_chapsim ()
  

  call Test_schemes()

  call Finalise_chapsim ()
  
end program

!##############################################################################
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

subroutine Finalise_chapsim()
  use mpi_mod
  implicit none

  !call Deallocate_all_variables
  !call Finalise_mpi()
  return
end subroutine Finalise_chapsim

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




