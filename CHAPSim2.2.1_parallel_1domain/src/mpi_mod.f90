!##############################################################################
module mpi_mod
  include "mpif.h"
  use decomp_2d

  integer :: nrow
  integer :: ncol
  integer :: ierror

  public :: Initialize_mpi, Finalise_mpi

contains 
!===============================================================================
!> \brief mpi initialisation.   
!>
!> this initialisation is a simple one.
!  only used before calling decomp_2d_init, 
!  where there is a complicted one used for 2-d decompoistion.
!  nrank = myid
!  nproc = size of processor 
!  both wil be replaced after calling decomp_2d_init
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d          domain type
!===============================================================================
  subroutine Initialize_mpi()
    implicit none
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)
    return
  end subroutine Initialize_mpi

  subroutine Finalise_mpi()  

    CALL MPI_FINALIZE(IERROR)

  end subroutine Finalise_mpi

end module mpi_mod
