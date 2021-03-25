module solver_tools_mod

  ! procedure
  private
  public  :: Compute_CFL_diffusion

  public  :: Calculate_massflux_from_velocity

contains
!===============================================================================
!===============================================================================
!> \brief Calculate the conservative variables from primary variable.     
!>
!> This subroutine is called to update $\rho u_i$ from $u_i$.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!> \param[in]     f             flow
!_______________________________________________________________________________
  subroutine Calculate_massflux_from_velocity(f, d)
    use parameters_constant_mod, only: ZERO
    use udf_type_mod
    use operations
    implicit none
!===============================================================================
! Arguments
!===============================================================================
    type(t_domain), intent(in )   :: d
    type(t_flow  ), intent(inout) :: f
!===============================================================================
! Local arguments
!===============================================================================
    real(WP), allocatable :: fi(:), fo(:)
    integer(4) :: i, j, k
!===============================================================================
! Code
!===============================================================================
!-------------------------------------------------------------------------------
! u1 -> g1
!_______________________________________________________________________________
    allocate ( fi( d%nc(1) ) ); fi = ZERO
    allocate ( fo( d%np(1) ) ); fo = ZERO
    do k = 1, d%nc(3)
      do j = 1, d%nc(2)
        fi(:) = f%dDens(:, j, k)
        call Get_midp_interpolation( 'x', 'C2P', d, fi(:), fo(:) )
        f%gx(:, j, k) = fo(:) * f%qx(:, j, k)
      end do
    end do
    deallocate (fi)
    deallocate (fo)
!-------------------------------------------------------------------------------
! u2 -> g2
!_______________________________________________________________________________
    allocate ( fi( d%nc(2) ) ); fi = ZERO
    allocate ( fo( d%np(2) ) ); fo = ZERO
    do k = 1, d%nc(3)
      do i = 1, d%nc(1)
        fi(:) = f%dDens(i, :, k)
        call Get_midp_interpolation( 'y', 'C2P', d, fi(:), fo(:) )
        f%gy(i, :, k) = fo(:) * f%qy(i, :, k)
      end do
    end do
    deallocate (fi)
    deallocate (fo)
!-------------------------------------------------------------------------------
! u3 -> g3
!_______________________________________________________________________________
    allocate ( fi( d%nc(3) ) ); fi = ZERO
    allocate ( fo( d%np(3) ) ); fo = ZERO
    do j = 1, d%nc(2)
      do i = 1, d%nc(1)
        fi(:) = f%dDens(i, j, :)
        call Get_midp_interpolation( 'z', 'C2P', d, fi(:), fo(:) )
        f%gz(i, j, :) = fo(:) * f%qz(i, j, :)
      end do
    end do
    deallocate (fi)
    deallocate (fo)

    return
  end subroutine Calculate_massflux_from_velocity

  subroutine Check_cfl_diffusion(x2r, rre)
    use input_general_mod, only: dt
    use parameters_constant_mod, only: TWO, ONE
    use precision_mod
    implicit none
    real(WP), intent(in) :: x2r(3)
    real(WP), intent(in) :: rre
    real(WP) :: cfl_diff

    ! check, ours is two times of the one in xcompact3d.
    cfl_diff = sum(x2r) * TWO * dt * rre

    write(*,*) "-------------------------------------------------------------------------------"
    if(cfl_diff > ONE) call Print_warning_msg("Warning: Diffusion number is larger than 1.")
    write(*,*) "Diffusion number :"
    write(*,"(12X, F13.8)") cfl_diff
    write(*,*) "-------------------------------------------------------------------------------"
    
    return
  end subroutine

  subroutine Check_cfl_convection(u, v, w, d)
    use input_general_mod, only: dt
    use parameters_constant_mod, only: ZERO, ONE
    use precision_mod
    use udf_type_mod, only: t_domain
    implicit none

    type(t_domain),               intent(in) :: d
    real(WP), dimension(:, :, :), intent(in) :: u, v, w

    real(WP), allocatable :: fi(:), fo(:)
    real(WP), allocatable :: udx(:, :, :)
    real(WP)              :: cfl_convection
    integer(4)            :: i, j, k

    allocate ( udx( d%nc(1), d%nc(2), d%nc(3) ) ); udx = ZERO

!-------------------------------------------------------------------------------
! \overline{u}^x/dx at cell centre
!_______________________________________________________________________________
    allocate ( fi( d%np(1) ) ); fi = ZERO
    allocate ( fo( d%nc(1) ) ); fo = ZERO
    do k = 1, d%nc(3)
      do j = 1, d%nc(2)
        fi(:) = u(:, j, k)
        call Get_midp_interpolation('x', 'P2C', d, fi(:), fo(:))
        udx(:, j, k) = fo(:) * d%h1r(3)
      end do
    end do
    deallocate (fi)
    deallocate (fo)
!-------------------------------------------------------------------------------
! \overline{v}^y/dy at cell centre
!_______________________________________________________________________________
    allocate ( fi( d%np(2) ) ); fi = ZERO
    allocate ( fo( d%nc(2) ) ); fo = ZERO
    do k = 1, d%nc(3)
      do i = 1, d%nc(1)
        fi(:) = v(i, :, k)
        call Get_midp_interpolation('y', 'P2C', d, fi(:), fo(:))
        udx(i, :, k) = udx(i, :, k) + fo(:) / d%yc(:)
      end do
    end do
    deallocate (fi)
    deallocate (fo)
!-------------------------------------------------------------------------------
! \overline{w}^z/dz at cell centre
!_______________________________________________________________________________
    allocate ( fi( d%np(3) ) ); fi = ZERO
    allocate ( fo( d%nc(3) ) ); fo = ZERO
    do j = 1, d%nc(2)
      do i = 1, d%nc(1)
        fi(:) = w(i, j, :)
        call Get_midp_interpolation('z', 'P2C', d, fi(:), fo(:))
        udx(i, :, k) = udx(i, :, k) + fo(:) * d%h1r(3)
      end do
    end do
    deallocate (fi)
    deallocate (fo)

    write(*,*) "-------------------------------------------------------------------------------"
    if(cfl_convection > ONE) call Print_warning_msg("Warning: CFL is larger than 1.")
    write(*,*) "CFL (convection) :"
    write(*,"(12X, F13.8)") cfl_convection
    write(*,*) "-------------------------------------------------------------------------------"
    
    deallocate (udx)

    return
  end subroutine

end module