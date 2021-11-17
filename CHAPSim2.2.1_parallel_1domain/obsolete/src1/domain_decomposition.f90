!##############################################################################
module domain_decomposition_mod
  use decomp_2d
  use input_general_mod
  implicit none

  integer :: ierror
  type(DECOMP_INFO) :: decomp_pcc ! eg, ux
  type(DECOMP_INFO) :: decomp_cpc ! eg, uy
  type(DECOMP_INFO) :: decomp_ccp ! eg, uz
  type(DECOMP_INFO) :: decomp_ccc ! eg, p

  type(DECOMP_INFO) :: decomp_ppc ! eg, <ux>^y, <uy>^x
  type(DECOMP_INFO) :: decomp_pcp ! eg, <ux>^z, <uz>^x
  type(DECOMP_INFO) :: decomp_cpp ! eg, <uy>^z, <uz>^y
  
  public :: Initialize_domain_decompsition

contains
!===============================================================================
!> \brief domain decompistion.   
!>
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d          domain type
!===============================================================================
  subroutine Initialize_domain_decompsition (d)
    use mpi_mod
    use udf_type_mod,      only : t_domain
    implicit none
    type(t_domain), intent(in)   :: d
!_______________________________________________________________________________
! basic 2D decompistion API
! limits: nrow <= min(nx, ny)
!         ncol <= min(ny, nz)
! xsize(i), ysize(i), zsize(i), i = 1,2,3 :
!   sizes of the sub-domains held by the current process. 
!   The first letter refers to the pencil orientation and the three 1D array elements 
!   contain the sub-domain sizes in X, Y and Z directions, respectively. 
!   example: xsize(1:3) means the subdomain size in x, y, z direction of x-pencil
!   In a 2D pencil decomposition, there is always one dimension which completely 
!   resides in local memory. So by definition, below relations hold 
!   xsize(1)==nx_global, ysize(2)==ny_global and zsize(3)==nz_global
! xstart(i), ystart(i), zstart(i), xend(i), yend(i), zend(i), i=1,2,3 :
!   the starting and ending indices for each sub-domain, as in the global coordinate system. 
!   Obviously, it can be seen that xsize(i)=xend(i)-xstart(i)+1. 
!   It may be convenient for certain applications to use global coordinate 
!   (for example when extracting a 2D plane from a 3D domain, it is easier to know which 
!   process owns the plane if global index is used).
!_______________________________________________________________________________
    call decomp_2d_init  (d%nc(1), d%nc(2), d%nc(3), nrow, ncol)
!_______________________________________________________________________________
    call decomp_info_init(d%np(1), d%nc(2), d%nc(3), decomp_pcc)
    call decomp_info_init(d%nc(1), d%np(2), d%nc(3), decomp_cpc)
    call decomp_info_init(d%nc(1), d%nc(2), d%np(3), decomp_ccp)
    call decomp_info_init(d%nc(1), d%nc(2), d%nc(3), decomp_ccc)

    call decomp_info_init(d%np(1), d%np(2), d%nc(3), decomp_ppc)
    call decomp_info_init(d%nc(1), d%np(2), d%np(3), decomp_cpp)
    call decomp_info_init(d%np(1), d%nc(2), d%np(3), decomp_pcp)

    d%ux_xst(1:3) = decomp_pcc%xst(1:3)
    d%ux_xen(1:3) = decomp_pcc%xen(1:3)
    d%ux_xsz(1:3) = decomp_pcc%xsz(1:3)
    d%uy_xst(1:3) = decomp_cpc%xst(1:3)
    d%uy_xen(1:3) = decomp_cpc%xen(1:3)
    d%uy_xsz(1:3) = decomp_cpc%xsz(1:3)
    d%uz_xst(1:3) = decomp_ccp%xst(1:3)
    d%uz_xen(1:3) = decomp_ccp%xen(1:3)
    d%uz_xsz(1:3) = decomp_ccp%xsz(1:3)
    d%ps_xst(1:3) = decomp_ccc%xst(1:3)
    d%ps_xen(1:3) = decomp_ccc%xen(1:3)
    d%ps_xsz(1:3) = decomp_ccc%xsz(1:3)

    d%ux_yst(1:3) = decomp_pcc%yst(1:3)
    d%ux_yen(1:3) = decomp_pcc%yen(1:3)
    d%ux_ysz(1:3) = decomp_pcc%ysz(1:3)
    d%uy_yst(1:3) = decomp_cpc%yst(1:3)
    d%uy_yen(1:3) = decomp_cpc%yen(1:3)
    d%uy_ysz(1:3) = decomp_cpc%ysz(1:3)
    d%uz_yst(1:3) = decomp_ccp%yst(1:3)
    d%uz_yen(1:3) = decomp_ccp%yen(1:3)
    d%uz_ysz(1:3) = decomp_ccp%ysz(1:3)
    d%ps_yst(1:3) = decomp_ccc%yst(1:3)
    d%ps_yen(1:3) = decomp_ccc%yen(1:3)
    d%ps_ysz(1:3) = decomp_ccc%ysz(1:3)

    d%ux_zst(1:3) = decomp_pcc%zst(1:3)
    d%ux_zen(1:3) = decomp_pcc%zen(1:3)
    d%ux_zsz(1:3) = decomp_pcc%zsz(1:3)
    d%uy_zst(1:3) = decomp_cpc%zst(1:3)
    d%uy_zen(1:3) = decomp_cpc%zen(1:3)
    d%uy_zsz(1:3) = decomp_cpc%zsz(1:3)
    d%uz_zst(1:3) = decomp_ccp%zst(1:3)
    d%uz_zen(1:3) = decomp_ccp%zen(1:3)
    d%uz_zsz(1:3) = decomp_ccp%zsz(1:3)
    d%ps_zst(1:3) = decomp_ccc%zst(1:3)
    d%ps_zen(1:3) = decomp_ccc%zen(1:3)
    d%ps_zsz(1:3) = decomp_ccc%zsz(1:3)

    d%dpcc = decomp_pcc
    d%dcpc = decomp_cpc
    d%dccp = decomp_ccp
    d%dccc = decomp_ccc

    d%dppc = decomp_ppc
    d%dcpp = decomp_cpp
    d%dpcp = decomp_pcp

#ifdef DEBUG
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    if(nrank == 0) write(*, *) 'Basic/Default 2D decompisiton : '
    write(*, *) 'In x-pencil, rank = ', nrank, ' size in x-dir = ', xsize(1), ' x id from ', xstart(1), xend(1)
    write(*, *) 'In x-pencil, rank = ', nrank, ' size in y-dir = ', xsize(2), ' y id from ', xstart(2), xend(2)
    write(*, *) 'In x-pencil, rank = ', nrank, ' size in z-dir = ', xsize(3), ' z id from ', xstart(3), xend(3)
    write(*, *) 'In y-pencil, rank = ', nrank, ' size in x-dir = ', ysize(1), ' x id from ', ystart(1), yend(1)
    write(*, *) 'In y-pencil, rank = ', nrank, ' size in y-dir = ', ysize(2), ' y id from ', ystart(2), yend(2)
    write(*, *) 'In y-pencil, rank = ', nrank, ' size in z-dir = ', ysize(3), ' z id from ', ystart(3), yend(3)
    write(*, *) 'In z-pencil, rank = ', nrank, ' size in x-dir = ', zsize(1), ' x id from ', zstart(1), zend(1)
    write(*, *) 'In z-pencil, rank = ', nrank, ' size in y-dir = ', zsize(2), ' y id from ', zstart(2), zend(2)
    write(*, *) 'In z-pencil, rank = ', nrank, ' size in z-dir = ', zsize(3), ' z id from ', zstart(3), zend(3)
#endif

#ifdef DEBUG
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    write(*, *) 'Ux in x-pencil, rank = ', nrank, ' size in x-dir = ', d%ux_xsz(1), ' x id = ', d%ux_xst(1), ' to ', d%ux_xen(1)
    write(*, *) 'Ux in x-pencil, rank = ', nrank, ' size in y-dir = ', d%ux_xsz(2), ' y id = ', d%ux_xst(2), ' to ', d%ux_xen(2)
    write(*, *) 'Ux in x-pencil, rank = ', nrank, ' size in z-dir = ', d%ux_xsz(3), ' z id = ', d%ux_xst(3), ' to ', d%ux_xen(3)
    write(*, *) 'Uy in x-pencil, rank = ', nrank, ' size in x-dir = ', d%uy_xsz(1), ' x id = ', d%uy_xst(1), ' to ', d%uy_xen(1)
    write(*, *) 'Uy in x-pencil, rank = ', nrank, ' size in y-dir = ', d%uy_xsz(2), ' y id = ', d%uy_xst(2), ' to ', d%uy_xen(2)
    write(*, *) 'Uy in x-pencil, rank = ', nrank, ' size in z-dir = ', d%uy_xsz(3), ' z id = ', d%uy_xst(3), ' to ', d%uy_xen(3)
    write(*, *) 'Uz in x-pencil, rank = ', nrank, ' size in x-dir = ', d%uz_xsz(1), ' x id = ', d%uz_xst(1), ' to ', d%uz_xen(1)
    write(*, *) 'Uz in x-pencil, rank = ', nrank, ' size in y-dir = ', d%uz_xsz(2), ' y id = ', d%uz_xst(2), ' to ', d%uz_xen(2)
    write(*, *) 'Uz in x-pencil, rank = ', nrank, ' size in z-dir = ', d%uz_xsz(3), ' z id = ', d%uz_xst(3), ' to ', d%uz_xen(3)
    write(*, *) 'p  in x-pencil, rank = ', nrank, ' size in x-dir = ', d%ps_xsz(1), ' x id = ', d%ps_xst(1), ' to ', d%ps_xen(1)
    write(*, *) 'p  in x-pencil, rank = ', nrank, ' size in y-dir = ', d%ps_xsz(2), ' y id = ', d%ps_xst(2), ' to ', d%ps_xen(2)
    write(*, *) 'p  in x-pencil, rank = ', nrank, ' size in z-dir = ', d%ps_xsz(3), ' z id = ', d%ps_xst(3), ' to ', d%ps_xen(3)
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    write(*, *) 'Ux in y-pencil, rank = ', nrank, ' size in x-dir = ', d%ux_ysz(1), ' x id = ', d%ux_yst(1), ' to ', d%ux_yen(1)
    write(*, *) 'Ux in y-pencil, rank = ', nrank, ' size in y-dir = ', d%ux_ysz(2), ' y id = ', d%ux_yst(2), ' to ', d%ux_yen(2)
    write(*, *) 'Ux in y-pencil, rank = ', nrank, ' size in z-dir = ', d%ux_ysz(3), ' z id = ', d%ux_yst(3), ' to ', d%ux_yen(3)
    write(*, *) 'Uy in y-pencil, rank = ', nrank, ' size in x-dir = ', d%uy_ysz(1), ' x id = ', d%uy_yst(1), ' to ', d%uy_yen(1)
    write(*, *) 'Uy in y-pencil, rank = ', nrank, ' size in y-dir = ', d%uy_ysz(2), ' y id = ', d%uy_yst(2), ' to ', d%uy_yen(2)
    write(*, *) 'Uy in y-pencil, rank = ', nrank, ' size in z-dir = ', d%uy_ysz(3), ' z id = ', d%uy_yst(3), ' to ', d%uy_yen(3)
    write(*, *) 'Uz in y-pencil, rank = ', nrank, ' size in x-dir = ', d%uz_ysz(1), ' x id = ', d%uz_yst(1), ' to ', d%uz_yen(1)
    write(*, *) 'Uz in y-pencil, rank = ', nrank, ' size in y-dir = ', d%uz_ysz(2), ' y id = ', d%uz_yst(2), ' to ', d%uz_yen(2)
    write(*, *) 'Uz in y-pencil, rank = ', nrank, ' size in z-dir = ', d%uz_ysz(3), ' z id = ', d%uz_yst(3), ' to ', d%uz_yen(3)
    write(*, *) 'p  in y-pencil, rank = ', nrank, ' size in x-dir = ', d%ps_ysz(1), ' x id = ', d%ps_yst(1), ' to ', d%ps_yen(1)
    write(*, *) 'p  in y-pencil, rank = ', nrank, ' size in y-dir = ', d%ps_ysz(2), ' y id = ', d%ps_yst(2), ' to ', d%ps_yen(2)
    write(*, *) 'p  in y-pencil, rank = ', nrank, ' size in z-dir = ', d%ps_ysz(3), ' z id = ', d%ps_yst(3), ' to ', d%ps_yen(3)
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    write(*, *) 'Ux in z-pencil, rank = ', nrank, ' size in x-dir = ', d%ux_zsz(1), ' x id = ', d%ux_zst(1), ' to ', d%ux_zen(1)
    write(*, *) 'Ux in z-pencil, rank = ', nrank, ' size in y-dir = ', d%ux_zsz(2), ' y id = ', d%ux_zst(2), ' to ', d%ux_zen(2)
    write(*, *) 'Ux in z-pencil, rank = ', nrank, ' size in z-dir = ', d%ux_zsz(3), ' z id = ', d%ux_zst(3), ' to ', d%ux_zen(3)
    write(*, *) 'Uy in z-pencil, rank = ', nrank, ' size in x-dir = ', d%uy_zsz(1), ' x id = ', d%uy_zst(1), ' to ', d%uy_zen(1)
    write(*, *) 'Uy in z-pencil, rank = ', nrank, ' size in y-dir = ', d%uy_zsz(2), ' y id = ', d%uy_zst(2), ' to ', d%uy_zen(2)
    write(*, *) 'Uy in z-pencil, rank = ', nrank, ' size in z-dir = ', d%uy_zsz(3), ' z id = ', d%uy_zst(3), ' to ', d%uy_zen(3)
    write(*, *) 'Uz in z-pencil, rank = ', nrank, ' size in x-dir = ', d%uz_zsz(1), ' x id = ', d%uz_zst(1), ' to ', d%uz_zen(1)
    write(*, *) 'Uz in z-pencil, rank = ', nrank, ' size in y-dir = ', d%uz_zsz(2), ' y id = ', d%uz_zst(2), ' to ', d%uz_zen(2)
    write(*, *) 'Uz in z-pencil, rank = ', nrank, ' size in z-dir = ', d%uz_zsz(3), ' z id = ', d%uz_zst(3), ' to ', d%uz_zen(3)
    write(*, *) 'p  in z-pencil, rank = ', nrank, ' size in x-dir = ', d%ps_zsz(1), ' x id = ', d%ps_zst(1), ' to ', d%ps_zen(1)
    write(*, *) 'p  in z-pencil, rank = ', nrank, ' size in y-dir = ', d%ps_zsz(2), ' y id = ', d%ps_zst(2), ' to ', d%ps_zen(2)
    write(*, *) 'p  in z-pencil, rank = ', nrank, ' size in z-dir = ', d%ps_zsz(3), ' z id = ', d%ps_zst(3), ' to ', d%ps_zen(3)
#endif

    return
  end subroutine Initialize_domain_decompsition

end module domain_decomposition_mod
