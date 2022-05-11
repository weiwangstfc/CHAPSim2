!##############################################################################
module mpi_mod
  !include "mpif.h"
  use MPI
  use decomp_2d
  use iso_fortran_env
  implicit none
  integer :: ierror
  integer :: nxdomain

  public :: Initialize_mpi
  public :: Finalise_mpi

contains 
!===============================================================================
  subroutine Initialize_mpi()
    call MPI_INIT(IERROR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, IERROR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, IERROR)
    return
  end subroutine Initialize_mpi
!===============================================================================
  subroutine Finalise_mpi()  
    call MPI_FINALIZE(IERROR)
  end subroutine Finalise_mpi
end module mpi_mod
!##############################################################################
program convert_decomp_data_to_vtk
  use decomp_2d
  use decomp_2d_io
  use mpi_mod
  
  implicit none 
  integer, parameter :: p_row = 1
  integer, parameter :: p_col = 1
  real(mytype), allocatable, dimension(:,:,:) :: input_var

  character( len = 128) :: filename
  logical :: file_exists = .FALSE.
  integer :: outputunit
  integer :: i, j, k
  real(mytype) :: x, y, z

  write(*, *) 'Data postprocessing begins ... '

  ! initialisation 
  call Initialize_mpi

  ! prepare files
  filename = 'display_test.vtk'
  call SYSTEM('cp display_mesh_domain1.vtk '//trim(filename))

  call decomp_2d_init(nx,ny,nz,p_row,p_col)

  ! prepare files
  filename = 'display_test.vtk'
  call SYSTEM('cp display_mesh_domain1.vtk '//trim(filename))

  ! reading data based in xpencil
  write(*, *) 'Reading data ... '
  allocate(u1b(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
! read back to different arrays
  call decomp_2d_read_one(1,u1b,'display_ux_0.0000.dat')

  
  INQUIRE(FILE = trim(filename), exist = file_exists)

  if(.not.file_exists) then

    open(newunit = outputunit, file = trim(filename), action = "write", status = "new")
    write(outputunit, '(A)') '# vtk DataFile Version 2.0'
    write(outputunit, '(A)') 'xy_mesh'
    write(outputunit, '(A)') 'ASCII'
    write(outputunit, '(A)') 'DATASET STRUCTURED_GRID'
    write(outputunit, '(A, 3I10.1)') 'DIMENSIONS', 8, 8, 8
    write(outputunit, '(A, I10.1, 1X, A)') 'POINTS', 8 * 8 * 8, 'float'
    do k = 1, 8
      z = 0.1d0 * real(k - 1, mytype)
      do j = 1, 8
        y = 0.1d0 * real(j - 1, mytype)
        do i = 1, 8
          x = 0.1d0 * real(i - 1, mytype)
          write(outputunit, *) x, y, z
        end do
      end do
    end do

  else
    open(newunit = outputunit, file = trim(filename), action = "write", status = "old", position="append")
  end if

  ! write out points data...
  write(outputunit, '(A, 1I10.1)') 'CELL_DATA', 8 * 8 * 8
  write(outputunit, '(A, 1I10.1)') 'SCALARS ux float',  1
  write(outputunit, '(A)') 'LOOKUP_TABLE default'
  write(*, *) 'Writing data ... '
  do k = 1, 8
    do j = 1, 8
      do i = 1, 8
        write(outputunit, *) u1b(i, j, k)
      end do
    end do
  end do
  close(outputunit)


  open(newunit = outputunit, file = 'display_test2.vtk', action = "write", status = "new")
  write(outputunit, *)'<VTKFile type="RectilinearGrid" version="0.1"',' byte_order="LittleEndian">'
  write(outputunit, *)'  <RectilinearGrid WholeExtent=','"1 ',8,' 1 ',8,' 1 ',8,'">'
  write(outputunit, *)'    <Piece Extent=','"1 ',8,' 1 ',8,' 1 ',8,'">'
  write(outputunit, *)'      <Coordinates>'
  write(outputunit, *)'        <DataArray type="Float32"',' Name="X_COORDINATES"',' NumberOfComponents="1">'
  write(outputunit, *) (0.1d0 * real(i - 1, mytype), i = 1, 8)
  write(outputunit, *)'        </DataArray>'
  write(outputunit, *)'        <DataArray type="Float32"',' Name="Y_COORDINATES"',' NumberOfComponents="1">'
  write(outputunit, *) (0.1d0 * real(j - 1, mytype), j = 1, 8)
  write(outputunit, *)'        </DataArray>'
  write(outputunit, *)'        <DataArray type="Float32"',' Name="Z_COORDINATES"',' NumberOfComponents="1">'
  write(outputunit, *) (0.1d0 * real(k - 1, mytype), k = 1, 8)
  write(outputunit, *)'        </DataArray>'
  write(outputunit, *)'      </Coordinates>'
  write(outputunit, *)'      <PointData Scalars="scalar">'

  write(outputunit, *)'        <DataArray Name="ux" type="Float32"', ' NumberOfComponents="1"',' format="ascii">'
  write(outputunit, *) (((u1b(i,j,k),i=1,8),j=1,8),k=1,8)
  write(outputunit, *)'        </DataArray>'
 ! End of writing data
  write(outputunit, *)'      </PointData>'
  ! Write the end of the file
  write(outputunit, *)'    </Piece>'
  write(outputunit, *)'  </RectilinearGrid>'
  write(outputunit, *)'</VTKFile>'
  close(outputunit)


  deallocate (u1b)
  write(*, *) 'Code finished sucessfully.'
  call Finalise_mpi()

end program


