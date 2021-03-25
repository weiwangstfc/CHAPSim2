module boundary_conditions_mod


  public :: Apply_BC_velocity

contains

  subroutine Apply_BC_velocity (f, d)
    use parameters_constant_mod, only: ZERO
    use udf_type_mod, only: t_domain, t_flow
    use input_general_mod, only: IBC_UDIRICHLET
    implicit none
    type(t_domain), intent(in )   :: d
    type(t_flow  ), intent(inout) :: f

    integer :: m, n
!-------------------------------------------------------------------------------
! ux at i = 1 
!-------------------------------------------------------------------------------
    m = 1
    n = 1
    if(d%bc(m, n) == IBC_UDIRICHLET) then
      f%qx(1, :, :) = d%ubc(m, n)
    end if
!-------------------------------------------------------------------------------
! ux at i = np 
!-------------------------------------------------------------------------------
    m = 2
    n = 1
    if(d%bc(m, n) == IBC_UDIRICHLET) then
      f%qx(d%np(1), :, :) = d%ubc(m, n)
    end if    
!-------------------------------------------------------------------------------
! uy at j = 1 
!-------------------------------------------------------------------------------
    m = 1
    n = 2
    if(d%bc(m, n) == IBC_UDIRICHLET) then
      f%qy(:, 1, :) = d%ubc(m, n)
    end if
!-------------------------------------------------------------------------------
! uy at j = np 
!-------------------------------------------------------------------------------
    m = 2
    n = 2
    if(d%bc(m, n) == IBC_UDIRICHLET) then
      f%qy(:, d%np(2), :) = d%ubc(m, n)
    end if
!-------------------------------------------------------------------------------
! uz at k = 1 
!-------------------------------------------------------------------------------
    m = 1
    n = 3
    if(d%bc(m, n) == IBC_UDIRICHLET) then
      f%qz(:, :, 1) = d%ubc(m, n)
    end if
!-------------------------------------------------------------------------------
! uz at k = np 
!-------------------------------------------------------------------------------
    m = 2
    n = 3
    if(d%bc(m, n) == IBC_UDIRICHLET) then
      f%qz(:, :, d%np(3)) = d%ubc(m, n)
    end if

    return
  end subroutine

end module