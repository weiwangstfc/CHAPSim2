!##############################################################################
module domain_decomposition_mod
  use decomp_2d
  implicit none

  integer :: ierror
  type(DECOMP_INFO) :: phu, phv, phw, php, spp
  
  public :: Initialize_domain_decompsition
  public :: Initialize_mpi
  public :: Finalise_mpi

contains
!===============================================================================
!===============================================================================
  subroutine Initialize_domain_decompsition (d)
    use udf_type_mod,      only : t_domain
    use input_general_mod, only : p_row, p_col
    use decomp_2d
    implicit none
    type(t_domain), intent(in)   :: d

    call decomp_2d_init  (d%np(1), d%np(2), d%np(3), p_row, p_col)
    call decomp_info_init(d%np(1), d%nc(2), d%nc(3), phu)
    call decomp_info_init(d%nc(1), d%np(2), d%nc(3), phv)
    call decomp_info_init(d%nc(1), d%nc(2), d%np(3), phw)
    call decomp_info_init(d%nc(1), d%nc(2), d%nc(3), php)

    return
  end subroutine Initialize_domain_decompsition
!===============================================================================
!===============================================================================
  subroutine Initialize_mpi()
    use MPI
    implicit none
    call MPI_INIT(ierror)
    return
  end subroutine Initialize_mpi
!===============================================================================
!===============================================================================
  subroutine Finalise_mpi()
    use MPI
    implicit none
    CALL MPI_FINALIZE(ierror)
    return
  end subroutine Finalise_mpi

end module domain_decomposition_mod
