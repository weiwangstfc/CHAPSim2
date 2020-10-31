program chapsim
  use mpi_mod
  implicit none

  call Initialize_mpi()

  call Initialize_chapsim()


  call Finalise_mpi()
end program
!##############################################################################
subroutine Initialize_chapsim()
  use input_mod, only : Initialize_input
  use geometry_mod, only : Initialize_geometry_variables
  use time_stepping_mod, only : Set_timestepping_coefficients
  implicit none

  
  call Initialize_input ()
  call Initialize_domain_decompsition ()
  call Initialize_geometry_variables ()
  call Set_timestepping_coefficients ()
  !call Initialize_flow_variables ()

end subroutine Initialize_chapsim

!##############################################################################

