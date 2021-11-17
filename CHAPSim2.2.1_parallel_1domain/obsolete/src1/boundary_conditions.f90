module boundary_conditions_mod

  public :: Apply_BC_thermo
  public :: Apply_BC_velocity

contains

  subroutine Apply_BC_thermo(t, d)
    use input
!-------------------------------------------------------------------------------
!   Build up B.C. info, undimensional, constant temperature
!-------------------------------------------------------------------------------
    do i = 1, 2
      if( d%ibcx(5, i) = IBC_DIRICHLET ) then
        tpbcx(i)%t = d%fbcx(5, i) / t%t0Ref
        call tpbcx(i)%Refresh_thermal_properties_from_T_undim
      end if
      if( d%ibcy(5, i) = IBC_DIRICHLET ) then
        tpbcy(i)%t = d%fbcy(5, i) / t%t0Ref
        call tpbcy(i)%Refresh_thermal_properties_from_T_undim
      end if
      if( d%ibcz(5, i) = IBC_DIRICHLET ) then
        tpbcz(i)%t = d%fbcz(5, i) / t%t0Ref
        call tpbcz(i)%Refresh_thermal_properties_from_T_undim
      end if
    end do

    return
  end subroutine

  subroutine Apply_BC_velocity (ux, uy, uz, d)
    use precision_mod
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain, t_flow
    use input_general_mod, only : IBC_UDIRICHLET
    implicit none
    type(t_domain), intent(in )   :: d
    real(WP), intent(inout)       :: ux(:, :, :), &
                                     uy(:, :, :), &
                                     uz(:, :, :)
    integer :: m, n

!-------------------------------------------------------------------------------
! x-pencil : ux at i = 1 
!-------------------------------------------------------------------------------
    m = 1
    n = 1
    if(d%bc(m, n) == IBC_UDIRICHLET) then
      if(d%ux_xst(1) == 1) then
        ux(1, :, :) = d%ubc(m, n)
      end if
    end if
!-------------------------------------------------------------------------------
! x-pencil : ux at i = np 
!-------------------------------------------------------------------------------
    m = 2
    n = 1
    if(d%bc(m, n) == IBC_UDIRICHLET) then
      if(d%ux_xen(1) == d%np(1)) then
        ux(d%ux_xsz(1), :, :) = d%ubc(m, n)
      end if
    end if    
!-------------------------------------------------------------------------------
! x-pencil : uy at j = 1 
!-------------------------------------------------------------------------------
    m = 1
    n = 2
    if(d%bc(m, n) == IBC_UDIRICHLET) then
      if(d%uy_xst(2) == 1) then
        uy(:, 1, :) = d%ubc(m, n)
      end if
    end if
!-------------------------------------------------------------------------------
! x-pencil : uy at j = np 
!-------------------------------------------------------------------------------
    m = 2
    n = 2
    if(d%bc(m, n) == IBC_UDIRICHLET) then
      if(d%uy_xen(2) == d%np(2)) then
        uy(:, d%uy_xsz(2), :) = d%ubc(m, n)
      end if
    end if
!-------------------------------------------------------------------------------
! x-pencil : uz at k = 1 
!-------------------------------------------------------------------------------
    m = 1
    n = 3
    if(d%bc(m, n) == IBC_UDIRICHLET) then
      if(d%uz_xst(3) == 1) then
        uz(:, :, 1) = d%ubc(m, n)
      end if
    end if
!-------------------------------------------------------------------------------
! x-pencil : uz at k = np 
!-------------------------------------------------------------------------------
    m = 2
    n = 3
    if(d%bc(m, n) == IBC_UDIRICHLET) then
      if(d%uz_xen(3) == d%np(3)) then
        uz(:, :, d%uz_xsz(3)) = d%ubc(m, n)
      end if
    end if

    return
  end subroutine

end module