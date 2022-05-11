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
  subroutine read_instantanous_thermo_raw_data(tm, dm)
    use udf_type_mod
    use thermo_info_mod
    use mpi_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(in) :: tm

    if(nrank == 0) call Print_debug_start_msg("read instantanous thermo data ... ...")
    ipencil = 1
    call decomp_2d_read_one(ipencil, tm%dh,    'instantanous_dh_'//trim(int2str(tm%nrsttckpt))//'.dat', dm%dccc)
    call decomp_2d_read_one(ipencil, tm%tTemp, 'instantanous_te_'//trim(int2str(tm%nrsttckpt))//'.dat', dm%dccc)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine

end module 