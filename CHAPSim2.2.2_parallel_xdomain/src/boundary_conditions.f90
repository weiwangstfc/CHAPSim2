module boundary_conditions_mod

  public :: Apply_BC_thermo
  public :: Apply_BC_velocity

contains

  subroutine Apply_BC_thermo(dm, th)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in )   :: dm
    type(t_thermo), intent(inout) :: th

    integer :: i 

!-------------------------------------------------------------------------------
!   Build up B.C. info, undimensional, constant temperature
!-------------------------------------------------------------------------------
    do i = 1, 2
      if( dm%ibcx(5, i) = IBC_DIRICHLET ) then
        tpbcx(i)%t = dm%fbcx(5, i) / th%t0Ref
        call tpbcx(i)%Refresh_thermal_properties_from_T_undim
      end if
      if( dm%ibcy(5, i) = IBC_DIRICHLET ) then
        tpbcy(i)%t = dm%fbcy(5, i) / tm%t0Ref
        call tpbcy(i)%Refresh_thermal_properties_from_T_undim
      end if
      if( dm%ibcz(5, i) = IBC_DIRICHLET ) then
        tpbcz(i)%t = dm%fbcz(5, i) / tm%t0Ref
        call tpbcz(i)%Refresh_thermal_properties_from_T_undim
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
    if(dm%ibcx(m, n) == IBC_DIRICHLET) then
      if(dtmp%xst(m) == 1) then
        ux(dtmp%xst(m), :, :) = dm%fbcx(m, n)
      end if
    end if
!-------------------------------------------------------------------------------
!   ux at x-pencil , x-id = np
!-------------------------------------------------------------------------------
    n = 2
    if(dm%ibcx(m, n) == IBC_DIRICHLET) then
      if(dtmp%xen(m) == dm%np(m)) then
        ux(dtmp%xen(m), :, :) = dm%fbcx(m, n)
      end if
    end if    
!-------------------------------------------------------------------------------
!   uy at x-pencil , y-id = 1
!-------------------------------------------------------------------------------
    dtmp = dm%dcpc
    m = 2
    n = 1
    if(dm%ibcy(m, n) == IBC_DIRICHLET) then
      if(dtmp%xst(m) == 1) then
        uy(:, dtmp%xst(m), :) = dm%fbcy(m, n)
      end if
    end if
!-------------------------------------------------------------------------------
!   uy at x-pencil , y-id = np
!-------------------------------------------------------------------------------
    n = 2
    if(dm%ibcy(m, n) == IBC_DIRICHLET) then
      if(dtmp%xen(m) == dm%np(m)) then
        uy(:, dtmp%xsz(m), :) = dm%fbcy(m, n)
      end if
    end if
!-------------------------------------------------------------------------------
!   uz at x-pencil , y-id = 1
!-------------------------------------------------------------------------------
    dtmp = dm%dccp
    m = 3
    n = 1
    if(dm%ibcz(m, n) == IBC_DIRICHLET) then
      if(dtmp%xst(m) == 1) then
        uz(:, :, dtmp%xst(m)) = dm%fbcz(m, n)
      end if
    end if
!-------------------------------------------------------------------------------
!   uz at x-pencil , y-id = np
!-------------------------------------------------------------------------------
    n = 2
    if(dm%ibcz(m, n) == IBC_DIRICHLET) then
      if(dtmp%xen(m) == dm%np(m)) then
        uz(:, :, dtmp%xsz(m)) = dm%fbcz(m, n)
      end if
    end if

    return
  end subroutine

end module