module restart
  implicit none 

  public :: write_instantanous_raw_data
  public :: read_instantanous_raw_data

contains 
!===============================================================================
!===============================================================================
  subroutine write_instantanous_flow_data(fl, dm)
    use udf_type_mod
    use mpi_mod
    implicit none
    type(t_domain), intent(in)    :: dm
    type(t_flow),   intent(in) :: fl

    integer :: ipencil
    

    if(nrank == 0) call Print_debug_start_msg("write instantanous flow data ...")

    ipencil = 1
    call decomp_2d_write_one(ipencil, fl%qx,   'instantanous_ux_'//trim(int2str(fl%iteration))//'.dat', dm%dpcc)
    call decomp_2d_write_one(ipencil, fl%qy,   'instantanous_uy_'//trim(int2str(fl%iteration))//'.dat', dm%dcpc)
    call decomp_2d_write_one(ipencil, fl%qz,   'instantanous_uz_'//trim(int2str(fl%iteration))//'.dat', dm%dccp)
    call decomp_2d_write_one(ipencil, fl%pres, 'instantanous_pr_'//trim(int2str(fl%iteration))//'.dat', dm%dccc)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
!===============================================================================
!===============================================================================
  subroutine write_instantanous_thermo_data(tm, dm)
    use udf_type_mod
    use thermo_info_mod
    use mpi_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(in) :: tm

    integer :: ipencil
    

    if(nrank == 0) call Print_debug_start_msg("write instantanous thermal data ...")

    ipencil = 1
    call decomp_2d_write_one(ipencil, tm%dh,    'instantanous_dh_'//trim(int2str(tm%iteration))//'.dat', dm%dccc)
    call decomp_2d_write_one(ipencil, tm%tTemp, 'instantanous_te_'//trim(int2str(tm%iteration))//'.dat', dm%dccc)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
!===============================================================================
!===============================================================================
  subroutine read_instantanous_flow_raw_data(fl, dm)
    use udf_type_mod
    use mpi_mod
    implicit none
    type(t_domain), intent(in)    :: dm
    type(t_flow),   intent(in) :: fl

    integer :: ipencil

    if(nrank == 0) call Print_debug_start_msg("read instantanous flow data ... ...")
    ipencil = 1
    call decomp_2d_read_one(ipencil, fl%qx,   'instantanous_ux_'//trim(int2str(fl%nrsttckpt))//'.dat', dm%dpcc)
    call decomp_2d_read_one(ipencil, fl%qy,   'instantanous_uy_'//trim(int2str(fl%nrsttckpt))//'.dat', dm%dcpc)
    call decomp_2d_read_one(ipencil, fl%qz,   'instantanous_uz_'//trim(int2str(fl%nrsttckpt))//'.dat', dm%dccp)
    call decomp_2d_read_one(ipencil, fl%pres, 'instantanous_pr_'//trim(int2str(fl%nrsttckpt))//'.dat', dm%dccc)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
!===============================================================================
!===============================================================================
  subroutine restore_flow_variables_from_restart(fl, dm)
    use udf_type_mod
    use mpi_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(in) :: fl
    real(WP) :: ubulk
    
    call Apply_BC_velocity(dm, fl%qx, fl%qx, fl%qx)
    call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), dm%fbcy(:, 1), dm, dm%dpcc, fl%qx, ubulk)
    if(nrank == 0) then
        Call Print_debug_mid_msg("  The restarted mass flux is:")
        write (OUTPUT_UNIT, '(5X, A, 1ES13.5)') ' average[u(x,y,z)]_[x,y,z]: ', ubulk
    end if
    !-------------------------------------------------------------------------------
    ! to check maximum velocity
    !-------------------------------------------------------------------------------
    call Check_maximum_velocity(fl%qx, fl%qy, fl%qz)
    !-------------------------------------------------------------------------------
    ! to set up other parameters for flow only, which will be updated in thermo flow.
    !-------------------------------------------------------------------------------
    fl%pcor(:, :, :) = ZERO
    fl%dDens  (:, :, :) = ONE
    fl%mVisc  (:, :, :) = ONE
    fl%dDensm1(:, :, :) = ONE
    fl%dDensm2(:, :, :) = ONE

    return
  end subroutine
!===============================================================================
!===============================================================================
  subroutine read_instantanous_thermo_raw_data(tm, dm)
    use udf_type_mod
    use thermo_info_mod
    use mpi_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(in) :: tm

    if (dm%ithermo /= 1) return

    if(nrank == 0) call Print_debug_start_msg("read instantanous thermo data ... ...")
    ipencil = 1
    call decomp_2d_read_one(ipencil, tm%dh,    'instantanous_dh_'//trim(int2str(tm%nrsttckpt))//'.dat', dm%dccc)
    call decomp_2d_read_one(ipencil, tm%tTemp, 'instantanous_te_'//trim(int2str(tm%nrsttckpt))//'.dat', dm%dccc)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine

!===============================================================================
  subroutine restore_thermo_variables_from_restart(fl, tm, dm)
    use udf_type_mod
    use thermo_info_mod
    use mpi_mod
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(in) :: fl
    type(t_thermo), intent(in) :: tm

    if (dm%ithermo /= 1) return

    call Update_thermal_properties(fl, tm, dm)

    return
  end subroutine

end module 