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
  use domain_decompistion_mod
  use geometry_mod
  use time_stepping_mod
  use flow_variables_mod
  implicit none

  call Initialize_mpi()
  call Initialize_general_input ()
  call Initialize_thermo_input ()
  call Initialize_domain_decompsition ()
  call Initialize_geometry_variables ()
  call Set_timestepping_coefficients ()
  call Allocate_flow_variables ()

end subroutine Initialize_chapsim

subroutine Finalise_chapsim()
  use mpi_mod
  implicit none


  call Finalise_mpi()

end subroutine Finalise_chapsim




