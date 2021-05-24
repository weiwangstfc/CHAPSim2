
subroutine Display_vtk_slice(d, str, varnm, vartp, var0)
  use udf_type_mod
  use parameters_constant_mod, only: ZERO
  use operations, only: Get_midp_interpolation_1D
  implicit none
  type(t_domain), intent( in ) :: d
  integer(4) :: vartp
  character( len = *), intent( in ) :: str
  character( len = *), intent( in ) :: varnm
  real(WP), intent( in ) :: var0(:, :, :)

  real(WP), allocatable :: var1(:, :, :)
  real(WP), allocatable :: fi(:), fo(:)

  integer :: output_unit
  integer :: i, j, k
  integer :: nd1, nd2, nd3
  integer :: nc1, nc2, nc3
  real(WP) :: x, y, z
  character( len = 128) :: filename
  logical :: file_exists = .FALSE.

  nd1 = d%np_geo(1)
  nd2 = d%np_geo(2)
  nd3 = d%np_geo(3)

  nc1 = d%nc(1)
  nc2 = d%nc(2)
  nc3 = d%nc(3)

  if ( trim( str ) == 'xy' ) then
    nd3 = 1
    nc3 = 1
  end if
  if ( trim( str ) == 'yz' ) then
    nd1 = 1
    nc1 = 1
  end if
  if ( trim( str ) == 'zx' ) then
    nd2 = 1
    nc2 = 1
  end if

  filename = 'display_'//str//'.vtk'

  INQUIRE(FILE = trim(filename), exist = file_exists)

  if(.not.file_exists) then

    open(newunit = output_unit, file = trim(filename), action = "write", status = "new")
    write(output_unit, '(A)') '# vtk DataFile Version 2.0'
    write(output_unit, '(A)') str//'_mesh'
    write(output_unit, '(A)') 'ASCII'
    write(output_unit, '(A)') 'DATASET STRUCTURED_GRID'
    write(output_unit, '(A, 3I10.1)') 'DIMENSIONS', nd1, nd2, nd3
    write(output_unit, '(A, I10.1, X, A)') 'POINTS', nd1 * nd2 * nd3, 'float'
    do k = 1, nd3
      z = d%h(3) * real(k - 1, WP)
      do j = 1, nd2
        y = d%yp(j)
        do i = 1, nd1
          x = d%h(1) * real(i - 1, WP)
          write(output_unit, *) x, y, z
        end do
      end do
    end do

  else
    open(newunit = output_unit, file = trim(filename), action = "write", status = "old", position="append")
  end if
    ! data convert to cell centre data...

  allocate ( var1(d%nc(1), d%nc(2), d%nc(3)) ); var1 = ZERO

  if (vartp == 0) then
    ! no convert. for example, p, T, etc.
    var1(:, :, :) = var0(:, :, :)

  else if (vartp == 1 ) then
    ! convert np data to nc data
    
    allocate ( fi(size(var0, vartp)) ); fi = ZERO
    allocate ( fo(d%nc(1)) ); fo = ZERO

    do k = 1, nc3
      do j = 1, nc2
        fi(:) = var0(:, j, k)
        call Get_midp_interpolation_1D('x', 'P2C', d, fi(:), fo(:))
        var1(:, j, k) = fo(:)
      end do
    end do
    deallocate (fi)
    deallocate (fo)

  else if (vartp == 2) then

    ! convert np data to nc data
    
    allocate ( fi(size(var0, vartp)) ); fi = ZERO
    allocate ( fo(d%nc(2)) ); fo = ZERO

    do k = 1, nc3
      do i = 1, nc1
        fi(:) = var0(i, :, k)
        call Get_midp_interpolation_1D('y', 'P2C', d, fi(:), fo(:))
        var1(i, :, k) = fo(:)
      end do
    end do
    deallocate (fi)
    deallocate (fo)

  else if (vartp == 3) then
    ! convert np data to nc data
    
    allocate ( fi(size(var0, vartp)) ); fi = ZERO
    allocate ( fo(d%nc(3)) ); fo = ZERO

    do j = 1, nc2
      do i = 1, nc1
        fi(:) = var0(i, j, :)
        call Get_midp_interpolation_1D('z', 'P2C', d, fi(:), fo(:))
        var1(i, j, :) = fo(:)
      end do
    end do
    deallocate (fi)
    deallocate (fo)

  else 
    call Print_error_msg(" No such direction in Subroutine: "//"Display_vtk_slice")
  end if

  ! write out points data...
  write(output_unit, '(A, 1I10.1)') 'CELL_DATA', nc1 * nc2 * nc3
  write(output_unit, '(A, 1I10.1)') 'SCALARS '//varnm(:)//' float',  nc1 * nc2 * nc3
  write(output_unit, '(A)') 'LOOKUP_TABLE default'

  do k = 1, nc3
    do j = 1, nc2
      do i = 1, nc1
        write(output_unit, *) var1(i, j, k)
      end do
    end do
  end do

  close(output_unit)
  deallocate (var1)

  return

end subroutine Display_vtk_slice
