module solver_tools_mod

  ! procedure
  private
  public  :: Compute_CFL_diffusion
  public  :: Calculate_parameters_in_eqs
  public  :: Calculate_massflux_from_velocity
  public  :: Check_cfl_convection
  public  :: Check_cfl_diffusion

contains

subroutine Calculate_parameters_in_eqs(f, t, iter)
  use input_general_mod, only: ithermo, nIterFlow0, ren, renIni, lenRef
  use input_thermo_mod, only: tpRef0
  use parameters_constant_mod, only: GRAVITY
  implicit none
  type(t_flow),   intent(inout)         :: f
  type(t_thermo), intent(inout)         :: t
  integer(4),     intent(in ), optional :: iter  

  real(WP) :: u0

  if ( present(dim) ) then
    if(iter < nIterIniRen) then
      f%rre = ONE / renIni
    else
      f%rre = ONE / ren
    end if
  else 
    f%rre = ONE / ren
  end if

  if(ithermo == 1) then

    t%rPrRen = f%rre * tpRef0%k / tpRef0%m / tpRef0%cp

    u0 = ONE / f%rre * tpRef0%m / tpRef0%d / lenRef
    if (igravity == 0) then
      ! no gravity
      f%fgravity = ZERO
    else if (igravity == 1 .or. igravity == 2 .or. igravity == 3 ) then 
      ! flow/gravity same dirction
      f%fgravity =  lenRef / u0 / u0 * GRAVITY
    else if (igravity == -1 .or. igravity == -2 .or. igravity == -3 ) then 
      ! flow/gravity opposite dirction
      f%fgravity = -lenRef / u0 / u0 * GRAVITY
    else
      ! no gravity
      f%fgravity = ZERO
    end if

  end if

  return
end subroutine 

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
    real(WP), dimension( d%nc(1) ) :: fix
    real(WP), dimension( d%np(1) ) :: fox
    real(WP), dimension( d%nc(2) ) :: fiy
    real(WP), dimension( d%np(2) ) :: foy
    real(WP), dimension( d%nc(3) ) :: fiz
    real(WP), dimension( d%np(3) ) :: foz
    integer(4) :: i, j, k
!===============================================================================
! Code
!===============================================================================
!-------------------------------------------------------------------------------
! u1 -> g1
!_______________________________________________________________________________
    do k = 1, d%nc(3)
      do j = 1, d%nc(2)
        fix(:) = f%dDens(:, j, k)
        call Get_midp_interpolation_1D( 'x', 'C2P', d, fix(:), fox(:) )
        f%gx(:, j, k) = fox(:) * f%qx(:, j, k)
      end do
    end do
!-------------------------------------------------------------------------------
! u2 -> g2
!_______________________________________________________________________________
    do k = 1, d%nc(3)
      do i = 1, d%nc(1)
        fiy(:) = f%dDens(i, :, k)
        call Get_midp_interpolation_1D( 'y', 'C2P', d, fiy(:), foy(:) )
        f%gy(i, :, k) = foy(:) * f%qy(i, :, k)
      end do
    end do
!-------------------------------------------------------------------------------
! u3 -> g3
!_______________________________________________________________________________
    do j = 1, d%nc(2)
      do i = 1, d%nc(1)
        fiz(:) = f%dDens(i, j, :)
        call Get_midp_interpolation_1D( 'z', 'C2P', d, fiz(:), foz(:) )
        f%gz(i, j, :) = foz(:) * f%qz(i, j, :)
      end do
    end do

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
    use parameters_constant_mod, only: ZERO, ONE
    use precision_mod
    use udf_type_mod, only: t_domain
    use operations, only: Get_midp_interpolation_1D
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
        call Get_midp_interpolation_1D('x', 'P2C', d, fi(:), fo(:))
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
        call Get_midp_interpolation_1D('y', 'P2C', d, fi(:), fo(:))
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
        call Get_midp_interpolation_1D('z', 'P2C', d, fi(:), fo(:))
        udx(i, j, :) = udx(i, j, :) + fo(:) * d%h1r(3)
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
