!##############################################################################
module mpi_mod
  include "mpif.h"

  integer :: myid
  integer :: npar
  integer :: nrow ! ==> p_row
  integer :: ncol ! ==> p_col
  integer :: ierror

  public :: Initialize_mpi, Finalise_mpi

contains 

  subroutine Initialize_mpi()
    implicit none
    ! MPI_CART_CREAT is included in the module of decomp_2d_init
    call MPI_INIT(ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, npar, ierror)

  end subroutine Initialize_mpi

  subroutine Finalise_mpi()  

    CALL MPI_FINALIZE(IERROR)

  end subroutine Finalise_mpi

end module mpi_mod
