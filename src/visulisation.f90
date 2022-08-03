module visulisation_mod


  public :: view_data_in_rank

contains

  subroutine view_data_in_rank(var, dtmp, dm, str, iter)
    use decomp_2d
    use udf_type_mod
    use typeconvert_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    type(t_domain), intent(in) :: dm
    real(WP), dimension(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3)), intent(in) :: var
    character(2), intent(in) :: str
    integer, intent(in) :: iter

    logical :: file_exists = .FALSE.
    integer :: outputunit

    character(128) :: filename
    integer :: i, j, k
    real(WP) :: lx, ly, lz


    filename = 'check_'//str//'_rank'//trim(int2str(nrank))//'_iter'//trim(int2str(iter))//'.vtk'

    INQUIRE(FILE = trim(filename), exist = file_exists)

    if(.not.file_exists) then
  
      open(newunit = outputunit, file = trim(filename), action = "write", status = "new")
      write(outputunit, '(A)') '# vtk DataFile Version 2.0'
      write(outputunit, '(A)') 'check data'
      write(outputunit, '(A)') 'ASCII'
      write(outputunit, '(A)') 'DATASET RECTILINEAR_GRID'
      write(outputunit, '(A, 3I10.1)') 'DIMENSIONS', dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3)

      write(outputunit, '(A, I10.1, 1X, A)') 'X_COORDINATES', dtmp%xsz(1), 'float'
      do i = 1, dtmp%xsz(1)
          lx = dm%h(1) * real( i + dtmp%xst(1) - 1, WP)
          write(outputunit, '(1ES17.7E3)') lx
      end do

      write(outputunit, '(A, I10.1, 1X, A)') 'Y_COORDINATES', dtmp%xsz(2), 'float'
      do j = 1, dtmp%xsz(2)
        ly = dm%yp( j + dtmp%xst(2) - 1 )
        write(outputunit, '(1ES17.7E3)') ly
      end do

      write(outputunit, '(A, I10.1, 1X, A)') 'Z_COORDINATES', dtmp%xsz(3), 'float'
      do k = 1, dtmp%xsz(3)
          lz = dm%h(3) * real( k + dtmp%xst(3) - 1, WP)
          write(outputunit, '(1ES17.7E3)') lz
      end do
  
    else
      open(newunit = outputunit, file = trim(filename), action = "write", status = "old", position="append")
    end if
  
    write(outputunit, '(A, 1I10.1)') 'POINT_DATA', dtmp%xsz(1) * dtmp%xsz(2) * dtmp%xsz(3)
    write(outputunit, '(A)') 'SCALARS '//str//' float'
    write(outputunit, '(A)') 'LOOKUP_TABLE default'
    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
          do i = 1, dtmp%xsz(1)
              write(outputunit, '(1ES17.7E3)') var(i, j, k)
          end do
      end do
    end do

    close(outputunit)

    return
  end subroutine




end module