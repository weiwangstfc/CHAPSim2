module VTK_mod
  use parameters_constant_mod
  implicit none
  
  private 
  public :: Generate_vtk_mesh_slice

contains

  subroutine Generate_vtk_mesh_slice(n1, n2, x1, x2, str)
    integer, intent( in ) :: n1, n2
    real(WP), intent( in ) :: x1(:), x2(:)
    character(*), intent( in ) :: str

    integer :: output_unit
    integer :: i, j

    open(newunit = output_unit, file = 'mesh_'//str//'.vtk', action = "write", status = "replace")
    write(output_unit, '(A)') '# vtk DataFile Version 2.0'
    write(output_unit, '(A)') str//'_mesh'
    write(output_unit, '(A)') 'ASCII'
    write(output_unit, '(A)') 'DATASET STRUCTURED_GRID'
    write(output_unit, '(A, 3I10.1)') 'DIMENSIONS', n1, n2, 1
    write(output_unit, '(A, I10.1, X, A)') 'POINTS', n1 * n2, 'float'
    do j = 1, n2
      do i = 1, n1
          write(output_unit, *) x1(i), x2(j), ZERO
      end do
    end do
    close(output_unit)

  end subroutine
end module VTK_mod






