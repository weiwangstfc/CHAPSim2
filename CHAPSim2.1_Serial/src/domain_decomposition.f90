!##############################################################################
module domain_decomposition_mod
  use decomp_2d
  implicit none

  integer :: ierror
  type(DECOMP_INFO) :: physu, physv, physw, physs, specp
  
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
    call decomp_2d_init  (d%np(1), d%np(2), d%np(3), p_row, p_col)
!_______________________________________________________________________________
    call decomp_info_init(d%np(1), d%nc(2), d%nc(3), physu)
    call decomp_info_init(d%nc(1), d%np(2), d%nc(3), physv)
    call decomp_info_init(d%nc(1), d%nc(2), d%np(3), physw)
    call decomp_info_init(d%nc(1), d%nc(2), d%nc(3), physs)

#ifdef DEBUG
    !if(myid == 0) write(*, *) 'Basic/Default 2D decompisiton : '
    write(*, *) 'In x-pencil size in x-dir = ', xsize(1), ' x id from ', xstart(1), xend(1)
    write(*, *) 'In x-pencil size in y-dir = ', xsize(2), ' y id from ', xstart(2), xend(2)
    write(*, *) 'In x-pencil size in z-dir = ', xsize(3), ' z id from ', xstart(3), xend(3)
    write(*, *) 'In y-pencil size in x-dir = ', ysize(1), ' x id from ', ystart(1), yend(1)
    write(*, *) 'In y-pencil size in y-dir = ', ysize(2), ' y id from ', ystart(2), yend(2)
    write(*, *) 'In y-pencil size in z-dir = ', ysize(3), ' z id from ', ystart(3), yend(3)
    write(*, *) 'In z-pencil size in x-dir = ', zsize(1), ' x id from ', zstart(1), zend(1)
    write(*, *) 'In z-pencil size in y-dir = ', zsize(2), ' y id from ', zstart(2), zend(2)
    write(*, *) 'In z-pencil size in z-dir = ', zsize(3), ' z id from ', zstart(3), zend(3)
#endif

#ifdef DEBUG
    !if(myid == 0) write(*, *) 'Physical Field of U velocity, 2D decompisiton : '
    write(*, *) 'In x-pencil size in x-dir = ', physu%xsz(1), ' x id from ', physu%xst(1), physu%xen(1)
    write(*, *) 'In x-pencil size in y-dir = ', physu%xsz(2), ' y id from ', physu%xst(2), physu%xen(2)
    write(*, *) 'In x-pencil size in z-dir = ', physu%xsz(3), ' z id from ', physu%xst(3), physu%xen(3)
    write(*, *) 'In y-pencil size in x-dir = ', physu%ysz(1), ' x id from ', physu%yst(1), physu%yen(1)
    write(*, *) 'In y-pencil size in y-dir = ', physu%ysz(2), ' y id from ', physu%yst(2), physu%yen(2)
    write(*, *) 'In y-pencil size in z-dir = ', physu%ysz(3), ' z id from ', physu%yst(3), physu%yen(3)
    write(*, *) 'In z-pencil size in x-dir = ', physu%zsz(1), ' x id from ', physu%zst(1), physu%zen(1)
    write(*, *) 'In z-pencil size in y-dir = ', physu%zsz(2), ' y id from ', physu%zst(2), physu%zen(2)
    write(*, *) 'In z-pencil size in z-dir = ', physu%zsz(3), ' z id from ', physu%zst(3), physu%zen(3)
#endif

#ifdef DEBUG
    !if(myid == 0) write(*, *) 'Physical Field of V velocity, 2D decompisiton : '
    write(*, *) 'In x-pencil size in x-dir = ', physv%xsz(1), ' x id from ', physv%xst(1), physv%xen(1)
    write(*, *) 'In x-pencil size in y-dir = ', physv%xsz(2), ' y id from ', physv%xst(2), physv%xen(2)
    write(*, *) 'In x-pencil size in z-dir = ', physv%xsz(3), ' z id from ', physv%xst(3), physv%xen(3)
    write(*, *) 'In y-pencil size in x-dir = ', physv%ysz(1), ' x id from ', physv%yst(1), physv%yen(1)
    write(*, *) 'In y-pencil size in y-dir = ', physv%ysz(2), ' y id from ', physv%yst(2), physv%yen(2)
    write(*, *) 'In y-pencil size in z-dir = ', physv%ysz(3), ' z id from ', physv%yst(3), physv%yen(3)
    write(*, *) 'In z-pencil size in x-dir = ', physv%zsz(1), ' x id from ', physv%zst(1), physv%zen(1)
    write(*, *) 'In z-pencil size in y-dir = ', physv%zsz(2), ' y id from ', physv%zst(2), physv%zen(2)
    write(*, *) 'In z-pencil size in z-dir = ', physv%zsz(3), ' z id from ', physv%zst(3), physv%zen(3)

    !if(myid == 0) write(*, *) 'Physical Field of W velocity, 2D decompisiton : '
    write(*, *) 'In x-pencil size in x-dir = ', physw%xsz(1), ' x id from ', physw%xst(1), physw%xen(1)
    write(*, *) 'In x-pencil size in y-dir = ', physw%xsz(2), ' y id from ', physw%xst(2), physw%xen(2)
    write(*, *) 'In x-pencil size in z-dir = ', physw%xsz(3), ' z id from ', physw%xst(3), physw%xen(3)
    write(*, *) 'In y-pencil size in x-dir = ', physw%ysz(1), ' x id from ', physw%yst(1), physw%yen(1)
    write(*, *) 'In y-pencil size in y-dir = ', physw%ysz(2), ' y id from ', physw%yst(2), physw%yen(2)
    write(*, *) 'In y-pencil size in z-dir = ', physw%ysz(3), ' z id from ', physw%yst(3), physw%yen(3)
    write(*, *) 'In z-pencil size in x-dir = ', physw%zsz(1), ' x id from ', physw%zst(1), physw%zen(1)
    write(*, *) 'In z-pencil size in y-dir = ', physw%zsz(2), ' y id from ', physw%zst(2), physw%zen(2)
    write(*, *) 'In z-pencil size in z-dir = ', physw%zsz(3), ' z id from ', physw%zst(3), physw%zen(3)
#endif

#ifdef DEBUG
    !if(myid == 0) write(*, *) 'Physical Field of pressure, 2D decompisiton : '
    write(*, *) 'In x-pencil size in x-dir = ', physs%xsz(1), ' x id from ', physs%xst(1), physs%xen(1)
    write(*, *) 'In x-pencil size in y-dir = ', physs%xsz(2), ' y id from ', physs%xst(2), physs%xen(2)
    write(*, *) 'In x-pencil size in z-dir = ', physs%xsz(3), ' z id from ', physs%xst(3), physs%xen(3)
    write(*, *) 'In y-pencil size in x-dir = ', physs%ysz(1), ' x id from ', physs%yst(1), physs%yen(1)
    write(*, *) 'In y-pencil size in y-dir = ', physs%ysz(2), ' y id from ', physs%yst(2), physs%yen(2)
    write(*, *) 'In y-pencil size in z-dir = ', physs%ysz(3), ' z id from ', physs%yst(3), physs%yen(3)
    write(*, *) 'In z-pencil size in x-dir = ', physs%zsz(1), ' x id from ', physs%zst(1), physs%zen(1)
    write(*, *) 'In z-pencil size in y-dir = ', physs%zsz(2), ' y id from ', physs%zst(2), physs%zen(2)
    write(*, *) 'In z-pencil size in z-dir = ', physs%zsz(3), ' z id from ', physs%zst(3), physs%zen(3)
#endif

#ifdef DEBUG
    !if(myid == 0) write(*, *) 'Spectral domain, 2D decompisiton : '
    write(*, *) 'In x-pencil size in x-dir = ', specp%xsz(1), ' x id from ', specp%xst(1), specp%xen(1)
    write(*, *) 'In x-pencil size in y-dir = ', specp%xsz(2), ' y id from ', specp%xst(2), specp%xen(2)
    write(*, *) 'In x-pencil size in z-dir = ', specp%xsz(3), ' z id from ', specp%xst(3), specp%xen(3)
    write(*, *) 'In y-pencil size in x-dir = ', specp%ysz(1), ' x id from ', specp%yst(1), specp%yen(1)
    write(*, *) 'In y-pencil size in y-dir = ', specp%ysz(2), ' y id from ', specp%yst(2), specp%yen(2)
    write(*, *) 'In y-pencil size in z-dir = ', specp%ysz(3), ' z id from ', specp%yst(3), specp%yen(3)
    write(*, *) 'In z-pencil size in x-dir = ', specp%zsz(1), ' x id from ', specp%zst(1), specp%zen(1)
    write(*, *) 'In z-pencil size in y-dir = ', specp%zsz(2), ' y id from ', specp%zst(2), specp%zen(2)
    write(*, *) 'In z-pencil size in z-dir = ', specp%zsz(3), ' z id from ', specp%zst(3), specp%zen(3)
#endif

    return
  end subroutine Initialize_domain_decompsition

end module domain_decomposition_mod
