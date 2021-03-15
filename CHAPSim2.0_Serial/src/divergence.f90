module divergence

contains

  subroutine check_divergence(d)
    use flow_variables_mod
    use parameters_constant_mod
    implicit none
    type(t_domain), intent(in) :: d
    real(WP), allocatable :: fi(:), fo(:)
    real(WP), allocatable :: div(:, :, :)
    real(WP) :: divmax


    allocate ( div(d%nc(1), d%nc(2), d%nc(3)) ); div = ZERO

    ! to get du / dx at cell centre, P2C, unknow = nc
    nsz = size(qx, 1)
    allocate ( fi(nsz)     ); fi = ZERO
    allocate ( fo(d%nc(1)) ); fo = ZERO

    do k = 1, d%nc(3)
      do j = 1, d%nc(2)
        fix(:) = qx(:, j, k)
        call Get_1st_derivative('x', 'P2C', d, fi(:), fo(:) )
        div(:, j, k) = fox(:)
      end do
    end do
    deallocate (fi)
    deallocate (fo)

    ! to get dv / dy at cell centre, P2C, unknow = nc
    nsz = size(qx, 2)
    allocate ( fi(nsz)     ); fi = ZERO
    allocate ( fo(d%nc(2)) ); fo = ZERO
    do k = 1, d%nc(3)
      do i = 1, d%nc(1)
        fi(:) = qx(i, :, k)
        call Get_1st_derivative('y', 'P2C', d, fi(:), fo(:) )
        div(i, :, k) = div(i, :, k) + fo(:)
      end do
    end do
    deallocate (fi)
    deallocate (fo)

    ! to get dw / dz at cell centre, P2C, unknow = nc
    nsz = size(qx, 3)
    allocate ( fi(nsz)     ); fi = ZERO
    allocate ( fo(d%nc(3)) ); fo = ZERO
    do j = 1, d%nc(2)
      do i = 1, d%nc(1)
        fi(:) = qx(i, j, :)
        call Get_1st_derivative('z', 'P2C', d, fi(:), fo(:) )
        div(i, j, :) = div(i, j, :) + fo(:)
      end do
    end do
    deallocate (fi)
    deallocate (fo)

    ! to get the max. divergence
    divmax = MAXVAL(div(:, :, :))

    deallocate (div)
    
    return
  end subroutine


end module 