program chapsim
  implicit none

  call Initialize_chapsim()


  call Finalise_chapsim
end program
!##############################################################################
subroutine Initialize_chapsim()
use  geometry_variables_mod

  call Initialize_mpi
  call Initialize_input
  call Initialize_domain_decompsition
  call Initialize_geometry_variables

end subroutine Initialize_chapsim

!##############################################################################
subroutine Finalise_chapsim()
  use mpi_info_mod
  implicit none

  CALL MPI_FINALIZE(IERROR)
end subroutine Finalise_chapsim
