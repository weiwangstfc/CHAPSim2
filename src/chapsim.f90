program chapsim
  implicit none

  call Initialize_chapsim ()

  call Finalise_chapsim ()
  
end program
!##############################################################################
subroutine Initialize_chapsim()
  use mpi_mod
  use input_general_mod
  use input_thermo_mod
  use domain_decomposition_mod
  use geometry_initialisation_mod
  use flow_variables_mod
  implicit none

  call Initialize_mpi()
  call Initialize_general_input ()
  call Initialize_thermo_input ()
  call Initialize_geometry_variables ()
  call Initialize_domain_decompsition ()
  call Allocate_flow_variables ()

end subroutine Initialize_chapsim

subroutine Finalise_chapsim()
  use mpi_mod
  implicit none

  call Finalise_mpi()

end subroutine Finalise_chapsim




