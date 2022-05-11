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
    !   for x-pencil
    !   scale the given thermo b.c. in dimensional to undimensional
    !-------------------------------------------------------------------------------
    do i = 1, 2

      if( dm%ibcx(i, 5) == IBC_DIRICHLET ) then
        ! dimensional T --> undimensional T
        dm%fbcx(i, 5) = dm%fbcx(i, 5) / th%t0Ref 
        th%tpbcx(i)%t = dm%fbcx(i, 5)
        call th%tpbcx(i)%Refresh_thermal_properties_from_T_undim
      else if (dm%ibcx(i, 5) == IBC_NEUMANN) then
        ! dimensional heat flux (k*dT/dx) --> undimensional heat flux (k*dT/dx)
        dm%fbcx(i, 5) = dm%fbcx(i, 5) * th%lenRef / tpRef0%k / tpRef0%t 
      else
      end if

      if( dm%ibcy(i, 5) == IBC_DIRICHLET ) then
        ! dimensional T --> undimensional T
        dm%fbcy(i, 5) = dm%fbcy(i, 5) / th%t0Ref 
        th%tpbcy(i)%t = dm%fbcy(i, 5)
        call th%tpbcy(i)%Refresh_thermal_properties_from_T_undim
      else if (dm%ibcy(i, 5) == IBC_NEUMANN) then
        ! dimensional heat flux (k*dT/dy) --> undimensional heat flux (k*dT/dy)
        dm%fbcy(i, 5) = dm%fbcy(i, 5) * th%lenRef / tpRef0%k / tpRef0%t 
      else
      end if

      if( dm%ibcz(i, 5) == IBC_DIRICHLET ) then
        ! dimensional T --> undimensional T
        dm%fbcz(i, 5) = dm%fbcz(i, 5) / th%t0Ref 
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
!   ux at x-pencil
!-------------------------------------------------------------------------------
    m = 1
    dtmp = dm%dpcc
    do s = 1, 2
      if(dm%ibcx(s, m) == IBC_DIRICHLET) then
        if(dtmp%xst(m) == 1) then
          ux(1, :, :) = dm%fbcx(s, m)
        end if
        if(dtmp%xen(m) == dm%np(m)) then
          ux(dtmp%xsz(m), :, :) = dm%fbcx(s, m)
        end if
      end if
    end do
!-------------------------------------------------------------------------------
!   uy at x-pencil
!-------------------------------------------------------------------------------
    m = 2
    dtmp = dm%dcpc
    do s = 1, 2
      if(dm%ibcy(s, m) == IBC_DIRICHLET) then
        if(dtmp%xst(m) == 1) then
          uy(:, 1, :) = dm%fbcy(s, m)
        end if
        if(dtmp%xen(m) == dm%np(m)) then
          uy(:, dtmp%xsz(m), :) = dm%fbcy(s, m)
        end if
      end if
    end do
!-------------------------------------------------------------------------------
!   uz at x-pencil
!-------------------------------------------------------------------------------
    m = 3
    dtmp = dm%dccp
    do s = 1, 2
      if(dm%ibcz(s, m) == IBC_DIRICHLET) then
        if(dtmp%xst(m) == 1) then
          uz(:, :, 1) = dm%fbcz(s, m)
        end if
        if(dtmp%xen(m) == dm%np(m)) then
          uz(:, :, dtmp%xsz(m)) = dm%fbcz(n, m)
        end if
      end if
    end do

    return
  end subroutine

end module
