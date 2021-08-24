!##############################################################################
module domain_decomposition_mod
  use decomp_2d
  implicit none

  integer :: ierror
  type(DECOMP_INFO) :: phu, phv, phw, php, spp
  
  public :: Initialize_domain_decompsition

contains
!===============================================================================
!===============================================================================
  subroutine Initialize_domain_decompsition (d)
    use udf_type_mod,      only : t_domain
    use input_general_mod, only : p_row, p_col
    use decomp_2d
    implicit none
    type(t_domain), intent(in)   :: d
!_______________________________________________________________________________
! basic 2D decompistion API
! limits: p_row <= min(nx, ny)
!         p_col <= min(ny, nz)
! xsize(i), ysize(i), zsize(i), i = 1,2,3 :
! sizes of the sub-domains held by the current process. 
! The first letter refers to the pencil orientation and the three 1D array elements 
! contain the sub-domain sizes in X, Y and Z directions, respectively. 
! In a 2D pencil decomposition, there is always one dimension which completely 
! resides in local memory. So by definition 
! xsize(1)==nx_global, ysize(2)==ny_global and zsize(3)==nz_global
! xstart(i), ystart(i), zstart(i), xend(i), yend(i), zend(i), i=1,2,3 :
! the starting and ending indices for each sub-domain, as in the global coordinate system. 
! Obviously, it can be seen that xsize(i)=xend(i)-xstart(i)+1. 
! It may be convenient for certain applications to use global coordinate 
! (for example when extracting a 2D plane from a 3D domain, it is easier to know which 
! process owns the plane if global index is used).
!_______________________________________________________________________________
    call decomp_2d_init  (d%np(1), d%np(2), d%np(3), p_row, p_col)
!_______________________________________________________________________________
    call decomp_info_init(d%np(1), d%nc(2), d%nc(3), phu)
    call decomp_info_init(d%nc(1), d%np(2), d%nc(3), phv)
    call decomp_info_init(d%nc(1), d%nc(2), d%np(3), phw)
    call decomp_info_init(d%nc(1), d%nc(2), d%nc(3), php)

    return
  end subroutine Initialize_domain_decompsition

end module domain_decomposition_mod
