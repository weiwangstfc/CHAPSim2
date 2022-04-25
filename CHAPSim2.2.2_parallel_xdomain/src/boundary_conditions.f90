module boundary_conditions_mod

  public :: Convert_thermo_BC_from_dim_to_undim
  public :: Apply_BC_velocity

contains

  subroutine Convert_thermo_BC_from_dim_to_undim(dm, th)
    use parameters_constant_mod
    use udf_type_mod
    use thermo_info_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    type(t_thermo), intent(inout) :: th

    integer :: i 

!-------------------------------------------------------------------------------
!   Build up B.C. info, undimensional, constant temperature
!-------------------------------------------------------------------------------
    do i = 1, 2

      if( dm%ibcx(i, 5) == IBC_DIRICHLET ) then
        dm%fbcx(i, 5) = dm%fbcx(i, 5) / th%t0Ref ! dimensional T --> undimensional T
        th%tpbcx(i)%t = dm%fbcx(i, 5)
        call th%tpbcx(i)%Refresh_thermal_properties_from_T_undim
      else if (dm%ibcx(i, 5) == IBC_NEUMANN) then
        ! dimensional heat flux (k*dT/dx) --> undimensional heat flux (k*dT/dx)
        dm%fbcx(i, 5) = dm%fbcx(i, 5) * th%lenRef / tpRef0%k / tpRef0%t 
      else
      end if

      if( dm%ibcy(i, 5) == IBC_DIRICHLET ) then
        dm%fbcy(i, 5) = dm%fbcy(i, 5) / th%t0Ref ! dimensional T --> undimensional T
        th%tpbcy(i)%t = dm%fbcy(i, 5)
        call th%tpbcy(i)%Refresh_thermal_properties_from_T_undim
      else if (dm%ibcy(i, 5) == IBC_NEUMANN) then
        ! dimensional heat flux (k*dT/dy) --> undimensional heat flux (k*dT/dy)
        dm%fbcy(i, 5) = dm%fbcy(i, 5) * th%lenRef / tpRef0%k / tpRef0%t 
      else
      end if

      if( dm%ibcz(i, 5) == IBC_DIRICHLET ) then
        dm%fbcz(i, 5) = dm%fbcz(i, 5) / th%t0Ref ! dimensional T --> undimensional T
        th%tpbcz(i)%t = dm%fbcz(i, 5)
        call th%tpbcz(i)%Refresh_thermal_properties_from_T_undim
      else if (dm%ibcz(i, 5) == IBC_NEUMANN) then
        ! dimensional heat flux (k*dT/dz) --> undimensional heat flux (k*dT/dz)
        dm%fbcz(i, 5) = dm%fbcz(i, 5) * th%lenRef / tpRef0%k / tpRef0%t 
      else
      end if

    end do

    return
  end subroutine
!===============================================================================
!> \brief Apply b.c. conditions 
!------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   public
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[in]     d             domain
!> \param[out]    f             flow
!===============================================================================
  subroutine Apply_BC_velocity (dm, ux, uy, uz)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in )   :: dm
    real(WP), intent(inout)       :: ux(:, :, :), &
                                     uy(:, :, :), &
                                     uz(:, :, :)
    integer :: m, n
    type(DECOMP_INFO) :: dtmp

!-------------------------------------------------------------------------------
!   ux at x-pencil , x-id = 1
!-------------------------------------------------------------------------------
    dtmp = dm%dpcc
    m = 1
    n = 1
    if(dm%ibcx(n, m) == IBC_DIRICHLET) then
      if(dtmp%xst(m) == 1) then
        ux(dtmp%xst(m), :, :) = dm%fbcx(n, m)
      end if
    end if
!-------------------------------------------------------------------------------
!   ux at x-pencil , x-id = np
!-------------------------------------------------------------------------------
    n = 2
    if(dm%ibcx(n, m) == IBC_DIRICHLET) then
      if(dtmp%xen(m) == dm%np(m)) then
        ux(dtmp%xen(m), :, :) = dm%fbcx(n, m)
      end if
    end if    
!-------------------------------------------------------------------------------
!   uy at x-pencil , y-id = 1
!-------------------------------------------------------------------------------
    dtmp = dm%dcpc
    m = 2
    n = 1
    if(dm%ibcy(n, m) == IBC_DIRICHLET) then
      if(dtmp%xst(m) == 1) then
        uy(:, dtmp%xst(m), :) = dm%fbcy(n, m)
      end if
    end if
!-------------------------------------------------------------------------------
!   uy at x-pencil , y-id = np
!-------------------------------------------------------------------------------
    n = 2
    if(dm%ibcy(n, m) == IBC_DIRICHLET) then
      if(dtmp%xen(m) == dm%np(m)) then
        uy(:, dtmp%xsz(m), :) = dm%fbcy(n, m)
      end if
    end if
!-------------------------------------------------------------------------------
!   uz at x-pencil , y-id = 1
!-------------------------------------------------------------------------------
    dtmp = dm%dccp
    m = 3
    n = 1
    if(dm%ibcz(n, m) == IBC_DIRICHLET) then
      if(dtmp%xst(m) == 1) then
        uz(:, :, dtmp%xst(m)) = dm%fbcz(n, m)
      end if
    end if
!-------------------------------------------------------------------------------
!   uz at x-pencil , y-id = np
!-------------------------------------------------------------------------------
    n = 2
    if(dm%ibcz(n, m) == IBC_DIRICHLET) then
      if(dtmp%xen(m) == dm%np(m)) then
        uz(:, :, dtmp%xsz(m)) = dm%fbcz(n, m)
      end if
    end if

    return
  end subroutine

end module
