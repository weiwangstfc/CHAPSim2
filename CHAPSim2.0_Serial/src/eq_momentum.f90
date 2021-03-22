module eq_momentum_mod

contains

  subroutine Calculate_rhs_convection(f, d)

    
    real(WP), dimension( d%np(1), d%nc(2), d%nc(3) ) :: rhs_m1
    real(WP), dimension( d%nc(1), d%np(2), d%nc(3) ) :: rhs_m2
    real(WP), dimension( d%nc(1), d%nc(2), d%np(3) ) :: rhs_m3

    integer(4), parameter :: II = 1, JJ = 2, KK = 3
    real(WP), allocatable :: fi(:), fo(:)

    real(WP), dimension( d%nc(1), d%nc(2), d%nc(3) ) :: gxxa, qxxa
    real(WP), dimension( d%np(1), d%np(2), d%nc(3) ) :: gyxa, qyxa
    real(WP), dimension( d%np(1), d%nc(2), d%np(3) ) :: gzxa, qzxa

    real(WP), dimension( d%np(1), d%np(2), d%nc(3) ) :: gxya, qxya
    real(WP), dimension( d%nc(1), d%nc(2), d%nc(3) ) :: gyya, qyya
    real(WP), dimension( d%nc(1), d%np(2), d%np(3) ) :: gzya, qzya

    real(WP), dimension( d%np(1), d%nc(2), d%np(3) ) :: gxza, qxza
    real(WP), dimension( d%nc(1), d%np(2), d%np(3) ) :: gyza, qyza
    real(WP), dimension( d%nc(1), d%nc(2), d%nc(3) ) :: gzza, qzza

!-------------------------------------------------------------------------------
! interpolation
!_______________________________________________________________________________
    ! interpolation operation in x direction
    allocate ( fi( d%np(II) ) ); fi = ZERO
    allocate ( fo( d%np(II) ) ); fo = ZERO

    do k = 1, d%np(KK)
      do j = 1, d%np(JJ)

        ! qx, gx, x average
        if( k <= d%nc(KK) .and. j <= d%nc(JJ) ) then
          fi(:) = f%gx(:, j, k)
          call Get_midp_interpolation( 'x', 'P2C', d, fi(:), fo( 1 : d%nc(II) ) ) ! gx(i, j, k)
          gxxa(:, j, k) = fo( 1 : d%nc(II) )
          
          fi( 1 : d%np(II) ) = f%qx(:, j, k)
          call Get_midp_interpolation( 'x', 'P2C', d, fi(:), fo( 1 : d%nc(II) ) ) ! qx(i, j, k)
          qxxa(:, j, k) = fo( 1 : d%nc(II) )
        end if

        ! qy, gy, x average
        if( k <= d%nc(KK) ) then
          fi( 1 : d%nc(II) ) = f%gy(:, j, k)
          call Get_midp_interpolation('x', 'C2P', d, fi( 1 : d%nc(II) ), fo(:) ) ! gy(i', j', k)
          gyxa(:, j, k) = fo(:)

          fi( 1 : d%nc(II) ) = f%qy(:, j, k)
          call Get_midp_interpolation('x', 'C2P', d, fi( 1 : d%nc(II) ), fo(:) ) ! qy(i', j', k)
          qyxa(:, j, k) = fo(:)
        end if

        ! qz, gz, x average
        if( j <= d%nc(JJ) ) then
          fi( 1 : d%nc(II) ) = f%gz(:, j, k)
          call Get_midp_interpolation('x', 'C2P', d, fi( 1 : d%nc(II) ), fo(:) ) ! gz(i', j, k')
          gzxa(:, j, k) = fo(:)

          fi( 1 : d%nc(II) ) = f%qz(:, j, k)
          call Get_midp_interpolation('x', 'C2P', d, fi( 1 : d%nc(II) ), fo(:) ) ! qz(i', j, k')
          qzxa(:, j, k) = fo(:)
        end if


      end do
    end do
    deallocate (fi)
    deallocate (fo)


    ! interpolation  operation in y direction
    allocate ( fi( d%np(JJ) ) ); fi = ZERO
    allocate ( fo( d%np(JJ) ) ); fo = ZERO
    
    do k = 1, d%np(KK)
      do i = 1, d%np(II)

        ! qx, gx, y average
        if ( k<=d%nc(KK) ) then
          fi( 1 : d%nc(JJ) ) = f%gx(i, :, k)
          call Get_midp_interpolation('y', 'C2P', d, fi( 1 : d%nc(JJ) ), fo(:) ) ! gx(i', j', k)
          gxya(i, :, k) = fo(:)

          fi( 1 : d%nc(JJ) ) = f%qx(i, :, k)
          call Get_midp_interpolation('y', 'C2P', d, fi( 1 : d%nc(JJ) ), fo(:) ) ! qx(i', j', k)
          qxya(i, :, k) = fo(:)
        end if
        
        ! qy, gy, y average
        if ( k<=d%nc(KK) .and.  i<=d%nc(II) ) then
          fi(:) = f%gy(i, :, k)
          call Get_midp_interpolation('y', 'P2C', d, fi(:), fo( 1 : d%nc(JJ) ) ) ! gy(i, j, k)
          gyya(i, :, k) = fo( 1 : d%nc(JJ) )

          fi(:) = f%qy(i, :, k)
          call Get_midp_interpolation('y', 'P2C', d, fi(:), fo( 1 : d%nc(JJ) ) ) ! qy(i, j, k)
          qyya(i, :, k) = fo( 1 : d%nc(JJ) )
        end if

        !qz, gz, y average
        if ( i<=d%nc(II) ) then
          fi( 1 : d%nc(JJ) ) = f%gz(i, :, k)
          call Get_midp_interpolation('y', 'C2P', d, fi( 1 : d%nc(JJ) ), fo(:)) ! gz(i, j', k')
          gzya(i, :, k) = fo(:)
  
          fi( 1 : d%nc(JJ) ) = f%qz(i, :, k)
          call Get_midp_interpolation('y', 'C2P', d, fi( 1 : d%nc(JJ) ), fo(:)) ! qz(i, j', k')
          qzya(i, :, k) = fo(:)
        end if

      end do
    end do
    deallocate (fi)
    deallocate (fo)


    ! interpolation  operation in z direction
    allocate ( fi( d%np(KK) ) ); fi = ZERO
    allocate ( fo( d%np(KK) ) ); fo = ZERO

    do j = 1, d%np(JJ)
      do i = 1, d%np(II)

        ! qx, gx, z average
        if( j <= d%nc(JJ) ) then
          fi( 1 : d%nc(KK) ) = f%gx(i, j, :)
          call Get_midp_interpolation('z', 'C2P', d, fi( 1 : d%nc(KK) ), fo(:) ) ! gx(i', j, k')
          gxza(i, j, :) = fo(:)

          fi( 1 : d%nc(KK) ) = f%qx(i, j, :)
          call Get_midp_interpolation('z', 'C2P', d, fi( 1 : d%nc(KK) ), fo(:) ) ! qx(i', j, k')
          qxza(i, j, :) = fo(:)
        end if

        ! qy, gy, z average
        if( i <= d%nc(II) ) then
          fi( 1 : d%nc(KK) ) = f%gy(i, j, :)
          call Get_midp_interpolation('z', 'C2P', d, fi( 1 : d%nc(KK) ), fo(:) ) ! gy(i, j', k')
          gyza(i, j, :) = fo(:)
  
          fi(:) = f%qy(i, j, :)
          call Get_midp_interpolation('z', 'C2P', d, fi( 1 : d%nc(KK) ), fo(:) ) ! qy(i, j', k')
          qyza(i, j, :) = fo(:)
        end if

        ! qz, gz, z average
        if( i <= d%nc(II) .and. j <= d%nc(JJ) ) then
          fi(:) = f%gz(i, j, :)
          call Get_midp_interpolation('z', 'P2C', d, fi(:), fo( 1 : d%nc(KK) ) ) ! gz(i, j, k)
          gzza(i, j, :) = fo( 1 : d%nc(KK) )
  
          fi(:) = f%qz(i, j, :)
          call Get_midp_interpolation('z', 'P2C', d, fi(:), fo( 1 : d%nc(KK) ) ) ! qz(i, j, k)
          qzza(i, j, :) = fo( 1 : d%nc(KK) )
        end if

      end do
    end do
    deallocate (fi)
    deallocate (fo)

!-------------------------------------------------------------------------------
! 1st deriviative
!_______________________________________________________________________________
    rhs_m1(:, :, :) = ZERO
    rhs_m2(:, :, :) = ZERO
    rhs_m3(:, :, :) = ZERO

    ! 1st deriviate operation in x direction
    allocate ( fi( d%np(II) ) ); fi = ZERO
    allocate ( fo( d%np(II) ) ); fo = ZERO
    do k = 1, d%np(KK)
      do j = 1, d%np(JJ)

        ! for x convection
        if( k <= d%nc(KK) .and. j <= d%nc(JJ) ) then
          fi( 1 : d%nc(II) ) = gxxa(:, j, k) * qxxa(:, j, k)                !(i, j, k)
          call Get_1st_derivative('x', 'C2P', d, fi( 1 : d%nc(II) ), fo(:)) ! dx(gx * qx) at (i', j, k)
          rhs_m1(:, j, k) = fo(:)
        end if

        ! for y convection
        if( k <= d%nc(KK)) then
          fi(:) = gxya(:, j, k) * qyxa(:, j, k)                !(i', j', k)
          call Get_1st_derivative('x', 'P2C', d, fi(:), fo( 1 : d%nc(II) )) ! dx(gy * qy) at (i, j', k)
          rhs_m2(:, j, k) = fo( 1 : d%nc(II) )
        end if

        ! for z convection
        if( j <= d%nc(JJ)) then
          fi(:) = gxza(:, j, k) * qzxa(:, j, k)                !(i', j, k')
          call Get_1st_derivative('x', 'P2C', d, fi(:), fo( 1 : d%nc(II) )) ! dx(gx * qz) at (i, j, k')
          rhs_m3(:, j, k) = fo( 1 : d%nc(II) )
        end if

      end do
    end do
    deallocate (fi)
    deallocate (fo)


    ! 1st deriviate operation in y direction
    allocate ( fi( d%np(JJ) ) ); fi = ZERO
    allocate ( fo( d%np(JJ) ) ); fo = ZERO
    do k = 1, d%np(KK)
      do i = 1, d%np(II)

        ! for x convection
        if ( k <= d%nc(KK) ) then
          fi(:) = gyxa(i, :, k) * qxya(i, :, k)                              ! (i', j', k)
          call Get_1st_derivative('y', 'P2C', d, fi(:), fo( 1 : d%nc(JJ) ) ) ! dy(gy * qx) at (i', j, k)
          rhs_m1(:, j, k) = rhs_m1(:, j, k) + fo( 1 : d%nc(JJ) )
        end if

        ! for y convection
        if ( k <= d%nc(KK) .and. i <= d%nc(II)) then
          fi( 1 : d%nc(JJ) ) = gyya(i, :, k) * qyya(i, :, k)                ! (i, j, k)
          call Get_1st_derivative('y', 'C2P', d, fi( 1 : d%nc(JJ) ), fo(:)) ! dy(gy * qy) at (i, j', k)
         rhs_m2(i, :, k) = rhs_m2(i, :, k) + fo(:)
        end if

        ! for z convection
        if ( i <= d%nc(II)) then
          fi(:) = gyza(i, :, k) * qzya(i, :, k)                !(i, j', k')
          call Get_1st_derivative('y', 'P2C', d, fi(:), fo( 1 : d%nc(JJ) )) ! dy(gy * qz) at (i, j, k')
          rhs_m3(i, :, k) = rhs_m3(i, :, k) + fo( 1 : d%nc(JJ) )
        end if

      end do
    end do
    deallocate (fi)
    deallocate (fo)

    ! 1st deriviate operation in z direction
    allocate ( fi( d%np(KK) ) ); fi = ZERO
    allocate ( fo( d%np(KK) ) ); fo = ZERO
    do j = 1, d%np(JJ)
      do i = 1, d%np(II)

        ! for x convection
        if( j <= d%nc(JJ) ) then
          fi(:) = gzxa(i, :, k) * qxza(i, :, k)                              ! (i', j, k')
          call Get_1st_derivative('z', 'P2C', d, fi(:), fo( 1 : d%nc(KK) ) ) ! dz(gz * qx) at (i', j, k)
          rhs_m1(:, j, k) = rhs_m1(:, j, k) + fo( 1 : d%nc(KK) )
        end if

        ! for y convection
        if( i <= d%nc(II) ) then
          fi(:) = gzya(i, j, :) * qyza(i, j, :)                              ! (i, j', k')
          call Get_1st_derivative('z', 'P2C', d, fi(:), fo( 1 : d%nc(KK) ) ) ! dz(gz * qy) at (i, j', k)
          rhs_m2(i, j, :) = rhs_m2(i, j, :) + fo( 1 : d%nc(KK) )
        end if

        ! for z convection
        if( j <= d%nc(JJ) .and. i <= d%nc(II) ) then
          fi( 1 : d%nc(KK) ) = gzza(i, j, :) * qzza(i, j, :)                !(i, j, k)
          call Get_1st_derivative('z', 'C2P', d, fi( 1 : d%nc(KK) ), fo(:)) ! dz(gz * qz) at (i, j, k')
          rhs_m3(i, j, :) = rhs_m3(i, j, :) + fo(:)
        end if

      end do
    end do
    deallocate (fi)
    deallocate (fo)
 
    return
  end subroutine Calculate_rhs_convection







  subroutine Calculate_momentum_rhs()

    call Calculate_rhs_convection(, m1_rhs)
    call Calculate_rhs_diffusion(, m1_rhs)
    call Calculate_rhs_pressure_gradient(, m1_rhs)
    call Calculate_rhs_buoyancy_force(, m1_rhs)
    call Calculate_rhs_periodic_driven(, m1_rhs)

    return
  end subroutine



  subroutine Solve_momentum_eq

    call Calculate_momentum_rhs


  end subroutine

end module eq_momentum_mod