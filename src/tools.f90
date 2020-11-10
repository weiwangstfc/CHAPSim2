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


module table_index_locating_mod
  implicit none
  private 
  public :: Map_variables_from_list

contains 

  pure function Map_variables_from_list(a, alist, blist) result(b)
    use precision_mod
    use parameters_constant_mod, only : MINP, ONE
    implicit none
    real(WP), intent(in) :: a
    real(WP), intent(in) :: alist(:)
    real(WP), intent(in) :: blist(:)
    real(WP) :: b

    integer :: i1, i2, im
    real(WP) :: d1, dm
    real(WP) :: w1, w2

    i1 = 1
    i2 = size(alist)

    do while ( (i2 - i1) > 1)
      im = i1 + (i2 - i1) / 2
      d1 = alist(i1) - a
      dm = alist(im) - a
      if ( (d1 * dm) > MINP ) then
        i1 = im
      else
        i2 = im
      end if
    end do

    w1 = (alist(i2) - a) / (alist(i2) - alist(i1)) 
    w2 = ONE - w1

    b = w1 * blist(i1) + w2 * blist(i2)

  end function Map_variables_from_list
end module table_index_locating_mod




module sort_array_mod
  implicit none
  private 
  public :: Sort_array_small2big

contains 

  subroutine Sort_array_small2big(t, d, h)
    use precision_mod
    real(WP),intent(inout) :: t(:)
    real(WP),intent(inout) :: d(:)
    real(WP),intent(inout),optional :: h(:)

    integer :: i, k
    integer :: n, n1, n2
    real(WP) :: buf_t, buf_d, buf_h

    n = size(t)
    n1 = size(d)

    if (n1 /= n) return
    if(present(h)) then
      n2 = size(h)
      if(n2 /= n) return
    end if

    do i = 1, n
      k = minloc( t(i:n), dim = 1) + i - 1

      buf_t = t(i)
      buf_d = d(i)

      t(i) = t(k)
      d(i) = d(k)

      t(k) = buf_t
      d(k) = buf_d

      if( present(h) ) then
        buf_h = h(i)
        h(i) = h(k)
        h(k) = buf_h
      end if
    end do 

  end subroutine Sort_array_small2big
end module sort_array_mod









