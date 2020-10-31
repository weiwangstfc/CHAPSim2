!##############################################################################
module mpi_mod
  include "mpif.h"

  integer :: ierror
  integer :: nrank
  integer :: nproc

  public :: Initialize_mpi, Finalise_mpi

contains 

  subroutine Initialize_mpi()

    call MPI_INIT(ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)

  end subroutine Initialize_mpi

  subroutine Finalise_mpi()  

    CALL MPI_FINALIZE(IERROR)

  end subroutine Finalise_mpi

end module mpi_mod