module solver_tools_mod

  ! procedure
  private

  public  :: Calculate_massflux_from_velocity
  public  :: Check_cfl_convection
  public  :: Check_cfl_diffusion

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
    use parameters_constant_mod, only : ZERO
    use udf_type_mod
    use operations
    implicit none
    type(t_domain), intent(in )   :: d
    type(t_flow  ), intent(inout) :: f
    real(WP), dimension( d%nc(1) ) :: fix
    real(WP), dimension( d%np(1) ) :: fox
    real(WP), dimension( d%nc(2) ) :: fiy
    real(WP), dimension( d%np(2) ) :: foy
    real(WP), dimension( d%nc(3) ) :: fiz
    real(WP), dimension( d%np(3) ) :: foz
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3)) ::  uy_ypencil
    real(WP), dimension( d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3)) :: duy_ypencil
    real(WP), dimension( d%uz_ysz(1), d%uz_ysz(2), d%uz_ysz(3)) ::  uz_ypencil
    real(WP), dimension( d%uz_ysz(1), d%uz_ysz(2), d%uz_ysz(3)) :: duz_ypencil
    real(WP), dimension( d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3)) ::   d_ypencil
    
    real(WP), dimension( d%uz_zsz(1), d%uz_zsz(2), d%uz_zsz(3)) ::  uz_zpencil
    real(WP), dimension( d%uz_zsz(1), d%uz_zsz(2), d%uz_zsz(3)) :: duz_zpencil
    real(WP), dimension( d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3)) ::   d_zpencil
    
    integer(4) :: i, j, k
!-------------------------------------------------------------------------------
! Default x-pencil
!-------------------------------------------------------------------------------
    if(d%ux_xsz(1) /= d%np(1)) call Print_error_msg("Error, not X-pencil")
!-------------------------------------------------------------------------------
! x-pencil : u1 -> g1
!-------------------------------------------------------------------------------
    do k = 1, d%ux_xsz(3)
      do j = 1, d%ux_xsz(2)
        fix(:) = f%dDens(:, j, k)
        call Get_midp_interpolation_1D( 'x', 'C2P', d, fix(:), fox(:) )
        f%gx(:, j, k) = fox(:) * f%qx(:, j, k)
      end do
    end do
!-------------------------------------------------------------------------------
! x-pencil --> y-pencil
!-------------------------------------------------------------------------------
    call transpose_x_to_y(f%qy,    uy_ypencil, d%dcpc)
    call transpose_x_to_y(f%dDens,  d_ypencil, d%dccc)
    call transpose_x_to_y(f%qz,    uz_ypencil, d%dccp)
!-------------------------------------------------------------------------------
! y-pencil : u2 -> g2
!-------------------------------------------------------------------------------
    do k = 1, d%uy_ysz(3)
      do i = 1, d%uy_ysz(1)
        fiy(:) = d_ypencil(i, :, k)
        call Get_midp_interpolation_1D( 'y', 'C2P', d, fiy(:), foy(:) )
        duy_ypencil(i, :, k) = foy(:) * uy_ypencil
      end do
    end do
!-------------------------------------------------------------------------------
! y-pencil --> z-pencil
!-------------------------------------------------------------------------------
    call transpose_y_to_z( d_ypencil,  d_zpencil, d%dccc)
    call transpose_y_to_z(uz_ypencil, uz_zpencil, d%dccp)
!-------------------------------------------------------------------------------
! Z-pencil : u3 -> g3
!-------------------------------------------------------------------------------
    do j = 1, d%uz_zsz(2)
      do i = 1, d%uz_zsz(1)
        fiz(:) = d_zpencil(i, j, :)
        call Get_midp_interpolation_1D( 'z', 'C2P', d, fiz(:), foz(:) )
        duz_zpencil(i, j, :) = foz(:) * uz_zpencil(i, j, :)
      end do
    end do
!-------------------------------------------------------------------------------
! z-pencil --> y-pencil
!-------------------------------------------------------------------------------
    call transpose_z_to_y(duz_zpencil, duz_ypencil, d%dccp)
!-------------------------------------------------------------------------------
! y-pencil --> x-pencil
!-------------------------------------------------------------------------------
    call transpose_y_to_x(duz_ypencil, f%gz, d%dccp)
    call transpose_y_to_x(duy_ypencil, f%gy, d%dcpc)

    return
  end subroutine Calculate_massflux_from_velocity
!===============================================================================
!===============================================================================
  subroutine Check_cfl_diffusion(x2r, rre)
    use input_general_mod, only : dt
    use parameters_constant_mod, only : TWO, ONE
    use precision_mod
    implicit none
    real(WP), intent(in) :: x2r(3)
    real(WP), intent(in) :: rre
    real(WP) :: cfl_diff

    ! check, ours is two times of the one in xcompact3d.
    cfl_diff = sum(x2r) * TWO * dt * rre

    if(cfl_diff > ONE) call Print_warning_msg("Warning: Diffusion number is larger than 1.")
    write (OUTPUT_UNIT,*) "  Diffusion number :"
    write (OUTPUT_UNIT,"(12X, F13.8)") cfl_diff
    
    return
  end subroutine
!===============================================================================
!===============================================================================
  subroutine Check_cfl_convection(u, v, w, d)
    use parameters_constant_mod, only : ZERO, ONE
    use precision_mod
    use input_general_mod, only : dt
    use udf_type_mod, only : t_domain
    use operations, only : Get_midp_interpolation_1D
    use domain_decomposition_mod
    implicit none

    type(t_domain),               intent(in) :: d
    real(WP), dimension(:, :, :), intent(in) :: u, v, w

    real(WP), allocatable :: fi(:), fo(:)
    real(WP) :: udx_xpencil (d%ps_xsz(1), d%ps_xsz(2), d%ps_xsz(3))
    real(WP) :: udx_ypencil (d%ps_ysz(1), d%ps_ysz(2), d%ps_ysz(3))
    real(WP) :: udx_ypencil (d%ps_zsz(1), d%ps_zsz(2), d%ps_zsz(3))
    real(WP) ::   v_ypencil (d%uy_ysz(1), d%uy_ysz(2), d%uy_ysz(3))
    real(WP) ::   w_ypencil (d%uz_ysz(1), d%uz_ysz(2), d%uz_ysz(3))
    real(WP) ::   w_zpencil (d%uz_zsz(1), d%uz_zsz(2), d%uz_zsz(3))
    real(WP)   :: cfl_convection, cfl_convection_work
    integer(4) :: i, j, k


!-------------------------------------------------------------------------------
!   Ensure it is in x-pencil
!-------------------------------------------------------------------------------
    if(d%ux_xsz(1) /= d%np(1)) call Print_error_msg("Error, not X-pencil")
!-------------------------------------------------------------------------------
! X-pencil
!-------------------------------------------------------------------------------
    allocate ( udx_xpencil( d%ps_xsz(1), d%ps_xsz(2), d%ps_xsz(3) ) )
    udx_xpencil = ZERO
!-------------------------------------------------------------------------------
! X-pencil : \overline{u}^x/dx at cell centre
!-------------------------------------------------------------------------------
    allocate ( fi( d%np(1) ) ); fi = ZERO
    allocate ( fo( d%nc(1) ) ); fo = ZERO
    udx_pencil(:, :, :) = ZERO
    do k = 1, d%ux_xsz(3)
      do j = 1, d%ux_xsz(2)
        fi(:) = u(:, j, k)
        call Get_midp_interpolation_1D('x', 'P2C', d, fi(:), fo(:))
        udx_xpencil(:, j, k) = fo(:) * d%h1r(3) * dt
      end do
    end do
    deallocate (fi)
    deallocate (fo)
!-------------------------------------------------------------------------------
! Convert X-pencil to Y-Pencil
!-------------------------------------------------------------------------------
    call transpose_x_to_y(udx_xpencil, udx_ypencil, d%dccc)
    call transpose_x_to_y(v,             v_ypencil, d%dcpc)
    call transpose_x_to_y(w,             w_ypencil, d%dccp)
!-------------------------------------------------------------------------------
! Y-pencil : \overline{v}^y/dy at cell centre
!-------------------------------------------------------------------------------
    allocate ( fi( d%np(2) ) ); fi = ZERO
    allocate ( fo( d%nc(2) ) ); fo = ZERO
    do k = 1, d%uy_ysz(3)
      do i = 1, d%uy_ysz(1)
        fi(:) = v(i, :, k)
        call Get_midp_interpolation_1D('y', 'P2C', d, fi(:), fo(:))
        udx_ypencil(i, :, k) = udx_ypencil(i, :, k) + fo(:) * d%h1r(2) * dt
      end do
    end do
    deallocate (fi)
    deallocate (fo)
!-------------------------------------------------------------------------------
! Convert Y-pencil to Z-Pencil
!-------------------------------------------------------------------------------
    call transpose_y_to_z(udx_ypencil, udx_zpencil, d%dccc)
    call transpose_y_to_z(  w_ypencil,   w_zpencil, d%dccp)
!-------------------------------------------------------------------------------
! Z-pencil : \overline{w}^z/dz at cell centre
!-------------------------------------------------------------------------------
    allocate ( fi( d%np(3) ) ); fi = ZERO
    allocate ( fo( d%nc(3) ) ); fo = ZERO
    do j = 1, d%uz_zsz(2)
      do i = 1, d%uz_zsz(1)
        fi(:) = w_zpencil(i, j, :)
        call Get_midp_interpolation_1D('z', 'P2C', d, fi(:), fo(:))
        udx_zpencil(i, j, :) = udx_zpencil(i, j, :) + fo(:) * d%h1r(3) * dt
      end do
    end do
    deallocate (fi)
    deallocate (fo)
!-------------------------------------------------------------------------------
! Z-pencil : Find the maximum 
!-------------------------------------------------------------------------------
    cfl_convection = MAXVAL(udx_zpencil(:, :, :))
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(cfl_convection, cfl_convection_work, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)

    if(nrank == 0) then
      if(cfl_convection_work > ONE) call Print_warning_msg("Warning: CFL is larger than 1.")
      write (OUTPUT_UNIT,*) "  CFL (convection) :"
      write (OUTPUT_UNIT,"(12X, F13.8)") cfl_convection_work
    end if
    
    return
  end subroutine

end module
