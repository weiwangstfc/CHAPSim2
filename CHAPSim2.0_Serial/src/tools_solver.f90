module solver_tools_mod

  ! procedure
  private
  public  :: Compute_CFL_diffusion
  private :: Get_max_timestep


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
!> \param[in]     u1            velocity q1
!> \param[in]     u2            velocity q2
!> \param[in]     u3            velocity q3
!> \param[out]    g1            rho * u1
!> \param[out]    g2            rho * u2
!> \param[out]    g3            rho * u3
!> \param[in]     deny          density
!> \param[in]     d             domain
!_______________________________________________________________________________
  subroutine Calculate_massflux_from_velocity(u1, u2, u3, g1, g2, g3, den, d)
    use parameters_constant_mod, only: ZERO
    implicit none
!===============================================================================
! Arguments
!===============================================================================
    real(WP), intent( in )  :: u1(:, :, :)
    real(WP), intent( in )  :: u2(:, :, :)
    real(WP), intent( in )  :: u3(:, :, :)
    real(WP), intent( out ) :: g1(:, :, :)
    real(WP), intent( out ) :: g2(:, :, :)
    real(WP), intent( out ) :: g3(:, :, :)
    real(WP), intent( in )  :: den(:, :, :)
    type(t_domain), intent(in) :: d
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
        fi(:) = den(:, j, k)
        call Get_midp_interpolation( 'x', 'C2P', d, fi(:), fo(:) )
        g1(:, j, k) = fo(:) * u1(:, j, k)
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
        fi(:) = den(i, :, k)
        call Get_midp_interpolation( 'y', 'C2P', d, fi(:), fo(:) )
        g2(i, :, k) = fo(:) * u2(i, :, k)
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
        fi(:) = den(i, j, :)
        call Get_midp_interpolation( 'z', 'C2P', d, fi(:), fo(:) )
        g1(i, j, :) = fo(:) * u3(i, j, :)
      end do
    end do
    deallocate (fi)
    deallocate (fo)

    return
  end subroutine Calculate_massflux_from_velocity


  subroutine Compute_CFL_diffusion(d, dtvis)
    implicit none
    type(t_domain), intent( in  ) :: d
    real(WP),       intent( out ) :: dtvis


    cfl



    dtvis = ZERO
    do j = 1, d%nc(2)
      dtvis = dtvis 

    end do

    return
  end subroutine


  subroutine Get_max_timestep



end module