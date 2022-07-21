module boundary_conditions_mod


  public :: Apply_BC_velocity

contains

  
!=============================================================================================================================================
!> \brief Apply b.c. conditions 
!--------------------------------------------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   public
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!---------------------------------------------------------------------------------------------------------------------------------------------
!> \param[in]     d             domain
!> \param[out]    f             flow
!=============================================================================================================================================
  subroutine Apply_BC_velocity (dm, ux, uy, uz)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in )   :: dm
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent(inout) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent(inout) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent(inout) :: uz
    integer :: m, s
    type(DECOMP_INFO) :: dtmp


    if(dm%ibcx(1, 1) /= IBC_DIRICHLET .and. &
       dm%ibcx(2, 1) /= IBC_DIRICHLET .and. &
       dm%ibcy(1, 2) /= IBC_DIRICHLET .and. &
       dm%ibcy(2, 2) /= IBC_DIRICHLET .and. &
       dm%ibcz(1, 3) /= IBC_DIRICHLET .and. &
       dm%ibcz(2, 3) /= IBC_DIRICHLET ) return

!---------------------------------------------------------------------------------------------------------------------------------------------
!   ux at x-pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------------------------
!   uy at x-pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------------------------
!   uz at x-pencil
!---------------------------------------------------------------------------------------------------------------------------------------------
    m = 3
    dtmp = dm%dccp
    do s = 1, 2
      if(dm%ibcz(s, m) == IBC_DIRICHLET) then
        if(dtmp%xst(m) == 1) then
          uz(:, :, 1) = dm%fbcz(s, m)
        end if
        if(dtmp%xen(m) == dm%np(m)) then
          uz(:, :, dtmp%xsz(m)) = dm%fbcz(s, m)
        end if
      end if
    end do

    return
  end subroutine

end module
