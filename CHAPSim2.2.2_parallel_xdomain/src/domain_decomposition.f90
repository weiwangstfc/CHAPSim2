!##############################################################################
module domain_decomposition_mod
  use mpi_mod
  use decomp_2d
  implicit none

  private :: Initialize_domain_decomposition
  public  :: Buildup_mpi_domain_decomposition

contains
!===============================================================================
!> \brief domain decompistion.   
!------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   priviate
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[in]     d          domain type
!===============================================================================
  subroutine Initialize_domain_decomposition (dm)
    use udf_type_mod
    use iso_fortran_env
    implicit none
    type(t_domain), intent(inout)   :: dm
    
#ifdef DEBUG
    type(DECOMP_INFO) :: dtmp
    integer :: i
#endif 
!-------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! initialize decomp 
!-------------------------------------------------------------------------------
    call decomp_info_init(dm%np(1), dm%nc(2), dm%nc(3), dm%dpcc)
    call decomp_info_init(dm%nc(1), dm%np(2), dm%nc(3), dm%dcpc)
    call decomp_info_init(dm%nc(1), dm%nc(2), dm%np(3), dm%dccp)
    call decomp_info_init(dm%nc(1), dm%nc(2), dm%nc(3), dm%dccc)
    call decomp_info_init(dm%np(1), dm%np(2), dm%nc(3), dm%dppc)
    call decomp_info_init(dm%nc(1), dm%np(2), dm%np(3), dm%dcpp)
    call decomp_info_init(dm%np(1), dm%nc(2), dm%np(3), dm%dpcp)

    call decomp_info_init(dm%np(1), dm%np(2), dm%np(3), dm%dppp) ! this is only used in test.

#ifdef DEBUG
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    do i = 1, 7
      select case(i)
      case(1)
        dtmp = dm%dpcc
        write (OUTPUT_UNIT, *) 'In the decomp - pcc grids (for ux, gx) :'
      case(2)
        dtmp = dm%dcpc
        write (OUTPUT_UNIT, *) 'In the decomp - pcc grids (for uy, gy) :'
      case(3)
        dtmp = dm%dccp
        write (OUTPUT_UNIT, *) 'In the decomp - cpp grids (for uz, gz) :'
      case(4)
        dtmp = dm%dccc
        write (OUTPUT_UNIT, *) 'In the decomp - ccc grids (for rho, p) :'
      case(5)
        dtmp = dm%dppc
        write (OUTPUT_UNIT, *) 'In the decomp - ppc grids (for dux/dy, duy/dx) :'
      case(6)
        dtmp = dm%dcpp
        write (OUTPUT_UNIT, *) 'In the decomp - ppc grids (for duy/dz, duz/dy) :'
      case(7)
        dtmp = dm%dpcp
        write (OUTPUT_UNIT, *) 'In the decomp - ppc grids (for dux/dz, duz/dx) :'
      case default
      end select
      call mpi_barrier(MPI_COMM_WORLD, ierror)
      write (OUTPUT_UNIT, *) 'x-pencil, x id in rank ', nrank, ' : ', dtmp%xst(1), dtmp%xen(1), dtmp%xsz(1)
      write (OUTPUT_UNIT, *) 'x-pencil, y id in rank ', nrank, ' : ', dtmp%xst(2), dtmp%xen(2), dtmp%xsz(2)
      write (OUTPUT_UNIT, *) 'x-pencil, z id in rank ', nrank, ' : ', dtmp%xst(3), dtmp%xen(3), dtmp%xsz(3)
      call mpi_barrier(MPI_COMM_WORLD, ierror)
      write (OUTPUT_UNIT, *) 'y-pencil, x id in rank ', nrank, ' : ', dtmp%yst(1), dtmp%yen(1), dtmp%ysz(1)
      write (OUTPUT_UNIT, *) 'y-pencil, y id in rank ', nrank, ' : ', dtmp%yst(2), dtmp%yen(2), dtmp%ysz(2)
      write (OUTPUT_UNIT, *) 'y-pencil, z id in rank ', nrank, ' : ', dtmp%yst(3), dtmp%yen(3), dtmp%ysz(3)
      call mpi_barrier(MPI_COMM_WORLD, ierror)
      write (OUTPUT_UNIT, *) 'z-pencil, x id in rank ', nrank, ' : ', dtmp%zst(1), dtmp%zen(1), dtmp%zsz(1)
      write (OUTPUT_UNIT, *) 'z-pencil, y id in rank ', nrank, ' : ', dtmp%zst(2), dtmp%zen(2), dtmp%zsz(2)
      write (OUTPUT_UNIT, *) 'z-pencil, z id in rank ', nrank, ' : ', dtmp%zst(3), dtmp%zen(3), dtmp%zsz(3)
    end do
#endif

    return
  end subroutine Initialize_domain_decomposition
!===============================================================================
!> \brief domain decompistion.  
!------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain   module
!>         all    once           all       public
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[in]     none          NA
!===============================================================================
  subroutine Buildup_mpi_domain_decomposition
    use vars_df_mod, only : domain
    implicit none
    integer :: i

    call decomp_2d_init(domain(1)%nc(1), domain(1)%nc(2), domain(1)%nc(3), nrow, ncol)

    do i = 1, nxdomain
      call Initialize_domain_decomposition(domain(i))
    end do

    return
  end subroutine Buildup_mpi_domain_decomposition

end module domain_decomposition_mod
