!##############################################################################
module mpi_mod
  include "mpif.h"

  integer :: ierr
  integer :: nrank
  integer :: nproc

  public :: Initialize_mpi, Finalise_mpi

contains 

  subroutine Initialize_mpi()
    implicit none
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, npar, ierr)

  end subroutine Initialize_mpi

  subroutine Finalise_mpi()  

    CALL MPI_FINALIZE(IERROR)

  end subroutine Finalise_mpi

end module mpi_mod