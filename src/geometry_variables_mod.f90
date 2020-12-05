module geometry_initialisation_mod
  implicit none

  private
  public :: Initialize_geometry_variables
  public :: Display_vtk_mesh

contains

  subroutine Display_vtk_mesh(nd, str)
    use geometry_mod
    type(node_t), intent( in ) :: nd(:, :, :)
    character( len = *), intent( in ) :: str

    integer :: output_unit
    integer :: i, j, k
    integer :: n1, n2, n3

    n1 = size( nd(:, 1, 1) )
    n2 = size( nd(1, :, 1) )
    n3 = size( nd(1, 1, :) )

    if ( trim( str ) == 'xy' ) then
      n3 = 1
    end if
    if ( trim( str ) == 'yz' ) then
      n1 = 1
    end if
    if ( trim( str ) == 'zx' ) then
      n2 = 1
    end if

    open(newunit = output_unit, file = 'mesh_'//str//'.vtk', action = "write", status = "replace")
    write(output_unit, '(A)') '# vtk DataFile Version 2.0'
    write(output_unit, '(A)') str//'_mesh'
    write(output_unit, '(A)') 'ASCII'
    write(output_unit, '(A)') 'DATASET STRUCTURED_GRID'
    write(output_unit, '(A, 3I10.1)') 'DIMENSIONS', n1, n2, n3
    write(output_unit, '(A, I10.1, X, A)') 'POINTS', n1 * n2 * n3, 'float'
    do k = 1, n3
      do j = 1, n2
        do i = 1, n1
          write(output_unit, *) nd(i, j, k)%x, nd(i, j, k)%y, nd(i, j, k)%z 
        end do
      end do
    end do
    close(output_unit)

  end subroutine Display_vtk_mesh

  subroutine Initialize_geometry_variables ()
    use mpi_mod
    use input_general_mod
    use geometry_mod
    use math_mod
    use chapsim_abort_mod
    use parameters_constant_mod, only : ONE, HALF, ZERO
    integer :: i, j, k
    real(WP) :: s, yy, c1, c2, c3, c4

    ! Build up domain info
    domain%bcx(:) = ifbcx(:)
    domain%bcy(:) = ifbcy(:)
    domain%bcz(:) = ifbcz(:)

    domain%is_periodic(:) = is_periodic(:)
    
    domain%nc(1) = ncx
    domain%nc(2) = ncy
    domain%nc(3) = ncz

    ! to note: geometric node number is always cell numbe + 1,
    ! different from the flow node number
    domain%np(1) = ncx + 1
    domain%np(2) = ncy + 1
    domain%np(3) = ncz + 1

    domain%dx = lxx / real(domain%nc(1), WP)
    domain%dz = lzz / real(domain%nc(3), WP)

    domain%dx2 = domain%dx * domain%dx
    domain%dz2 = domain%dz * domain%dz

    domain%dxi = ONE / domain%dx
    domain%dzi = ONE / domain%dz

    ! to build up node and cell info
    allocate ( node (domain%np(1), domain%np(2), domain%np(3)) )
    allocate ( cell (domain%nc(1), domain%nc(2), domain%nc(3)) )
    ! to initialize 
    do k = 1, domain%np(3)
      do j = 1, domain%np(2)
        do i = 1, domain%np(1)
          node(i, j, k)%x = ZERO
          node(i, j, k)%y = ZERO
          node(i, j, k)%z = ZERO
        end do
      end do
    end do
    
    do k = 1, domain%nc(3)
      do j = 1, domain%nc(2)
        do i = 1, domain%nc(1)
          cell(i, j, k)%x = ZERO
          cell(i, j, k)%y = ZERO
          cell(i, j, k)%z = ZERO
        end do
      end do
    end do
    
    ! x 
    block_xcoordinate: do i = 1, domain%np(1)
      node(i, :, :)%x = real( (i - 1), WP ) * domain%dx
    end do block_xcoordinate

    block_xcoordinate: do i = 1, domain%nc(1)
      cell(i, :, :)%x = node(i, :, :)%x + domain%dx * HALF
    end do block_xcoordinate

    ! z
    block_zcoordinate: do k = 1, domain%np(3)
      node(:, :, k)%z = real( (k - 1), WP ) * domain%dz
    end do block_zcoordinate

    block_zcoordinate: do k = 1, domain%nc(3)
      cell(:, :, k)%z = node(:, :, k)%z + domain%dz * HALF
    end do block_zcoordinate

    ! y
    block_ycnst: if (istret == ISTRET_SIDES) then
      c1 = rstret * HALF
      c2 = tanh_wp (c1)
      c3 = ONE
      c4 = HALF
    else if (istret == ISTRET_BOTTOM) then
      c1 = rstret * ONE
      c2 = tanh_wp (c1)
      c3 = ONE
      c4 = ONE
    else if (istret == ISTRET_TOP) then
      c1 = rstret * ZERO
      c2 = tanh_wp (rstret)
      c3 = ZERO
      c4 = ONE
    else
      c1 = rstret * HALF
      c2 = tanh_wp (c1)
      c3 = ONE
      c4 = HALF
    end if block_ycnst

    block_ynd: do j = 1, domain%np(2)
      yy = real ((j - 1), WP) / real ( (npy - 1), WP)
      if(istret == ISTRET_NO) then
        s = yy
      else 
        s = (tanh_wp( (rstret * yy) - c1 ) / c2 + c3) * c4
      end if
      node(:, j, :)%y = s * (lyt - lyb) + lyb
    end do block_ynd 

    block_ycl: do j = 1, domain%nc(2)
      cell(:, j, :)%y = ( node(1, j, 1)%y + node(1, j + 1, 1)%y ) * HALF
    end do block_ycl

    block_dy: do j = 1, domain%nc(2)
      cell(:, j, :)%dy = node(1, j + 1, 1)%y - node(1, j, 1)%y 
    end do block_dy

    ! to validate the y distance
    yy = 0
    do j = 1, domain%nc(2)
      yy = yy + cell(:, j, :)%dy
    end do
    if ( abs_wp(yy - (lyt - lyb)) > TRUNCERR ) &
    call Print_error_msg(1, "Error in meshing.")

    ! ri
    cell(:, :, :)%ri = ONE
    node(:, :, :)%ri = ONE
    if(icoordinate == ICYLINDRICAL) then

      do j = 1, domain%np(2)
        if ( abs_wp( node(:, j, :)%y ) < MINP ) then
          node(:, j, :)%ri = MAXP
        else
          node(:, j, :)%ri = ONE / node(:, j, :)%y
        end if
      end do 

      do j = 1, domain%nc(2)
        cell(:, j, :)%ri = ONE / cell(:, j, :)%y
      end do

    end if

    ! to print geometry/mesh domain
    if(myid == 0) then ! test
      call Display_vtk_mesh ( node(:, :, :), 'xy' )
      call Display_vtk_mesh ( node(:, :, :), 'yz' )
      call Display_vtk_mesh ( node(:, :, :), 'zx' )
      call Display_vtk_mesh ( node(:, :, :), '3D' )
    end if

  end subroutine  Initialize_geometry_variables

end module geometry_initialisation_mod


