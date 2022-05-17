!##############################################################################
module mpi_mod
  !include "mpif.h"
  use MPI
  use decomp_2d
  use iso_fortran_env
  implicit none
  integer :: ierror
  integer :: nxdomain
  integer :: p_row
  integer :: p_col

  public :: Initialize_mpi
  public :: Finalise_mpi

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
    call MPI_INIT(IERROR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, IERROR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, IERROR)
    return
  end subroutine Initialize_mpi
!===============================================================================
!===============================================================================
  subroutine Finalise_mpi()  
    call MPI_FINALIZE(IERROR)
  end subroutine Finalise_mpi

end module mpi_mod
