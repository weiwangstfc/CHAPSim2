module io_restart_mod
  implicit none 

  character(len=10), parameter :: io_name = "restart-io"

  private :: read_instantanous_array
  private :: write_instantanous_array
  public :: write_instantanous_flow
  public :: write_instantanous_thermo
  public :: read_instantanous_flow
  public :: read_instantanous_thermo
  public :: restore_flow_variables_from_restart
  public :: restore_thermo_variables_from_restart

contains 
!==========================================================================================================
!==========================================================================================================
  subroutine read_instantanous_array(var, keyword, idom, iter, dtmp)
    use precision_mod
    use files_io_mod
    use io_tools_mod
    use decomp_2d_io
    implicit none 
    integer, intent(in) :: idom
    character(*), intent(in) :: keyword
    integer, intent(in) :: iter
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), dimension(:, :, :), intent(out) :: var( dtmp%xsz(1), &
                                                      dtmp%xsz(2), &
                                                      dtmp%xsz(3))
    character(120):: data_flname

    call generate_file_name(data_flname, idom, trim(keyword), 'bin', iter)
    if(nrank == 0) call Print_debug_mid_msg("Reading "//trim(dir_data)//"/"//trim(data_flname))

    call decomp_2d_read_one(X_PENCIL, var, trim(dir_data), trim(data_flname), io_name, dtmp, reduce_prec=.false.)

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine write_instantanous_array(var, keyword, idom, iter, dtmp)
    use udf_type_mod
    use files_io_mod
    use io_tools_mod
    use decomp_2d_io
    implicit none 
    real(WP), intent(in) :: var( :, :, :)
    type(DECOMP_INFO), intent(in) :: dtmp
    character(*), intent(in) :: keyword
    integer, intent(in) :: idom
    integer, intent(in) :: iter
    
    character(120):: data_flname_path
    logical :: file_exists

    call generate_pathfile_name(data_flname_path, idom, trim(keyword), dir_data, 'bin', iter)
    call decomp_2d_write_one(X_PENCIL, var, trim(data_flname_path), dtmp)

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine write_instantanous_flow(fl, dm)
    use udf_type_mod
    use decomp_2d_io
    use io_tools_mod
    use files_io_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(in) :: fl

    character(120):: data_flname_path
    character(120):: keyword
    logical :: file_exists

    if(nrank == 0) call Print_debug_start_msg("writing out instantanous 3d flow data ...")

    call write_instantanous_array(fl%qx, 'ux', dm%idom, fl%iteration, dm%dpcc)
    call write_instantanous_array(fl%qy, 'qy', dm%idom, fl%iteration, dm%dcpc)
    call write_instantanous_array(fl%qz, 'qz', dm%idom, fl%iteration, dm%dccp)
    call write_instantanous_array(fl%pres, 'pr', dm%idom, fl%iteration, dm%dccc)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine write_instantanous_thermo(tm, dm)
    use udf_type_mod
    use thermo_info_mod
    use decomp_2d_io
    use io_tools_mod
    use files_io_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(in) :: tm

    character(120):: data_flname_path
    character(120):: keyword
    

    if(nrank == 0) call Print_debug_start_msg("writing out instantanous 3d thermo data ...")

    call write_instantanous_array(tm%dh,    'rhoh', dm%idom, tm%iteration, dm%dccc)
    call write_instantanous_array(tm%tTemp, 'temp', dm%idom, tm%iteration, dm%dccc)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine read_instantanous_flow(fl, dm)
    use udf_type_mod
    use decomp_2d_io
    use io_tools_mod
    use precision_mod
    use files_io_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    type(t_flow),   intent(inout) :: fl

    character(120):: data_flname
    character(120):: keyword


    if(nrank == 0) call Print_debug_start_msg("read instantanous flow data ... ...")

    call read_instantanous_array(fl%qx, 'ux', dm%idom, fl%iterfrom, dm%dpcc)
    call read_instantanous_array(fl%qy, 'uy', dm%idom, fl%iterfrom, dm%dcpc)
    call read_instantanous_array(fl%qz, 'uz', dm%idom, fl%iterfrom, dm%dccp)
    call read_instantanous_array(fl%pres, 'pr', dm%idom, fl%iterfrom, dm%dccc)
    
    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine restore_flow_variables_from_restart(fl, dm)
    use udf_type_mod
    use mpi_mod
    use parameters_constant_mod
    use boundary_conditions_mod
    use solver_tools_mod
    use wtformat_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(inout) :: fl
    real(WP) :: ubulk
    

    call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), dm%fbcy_qx(:, :, :), dm, dm%dpcc, fl%qx, ubulk, "ux")
    if(nrank == 0) then
        Call Print_debug_mid_msg("  The restarted mass flux is:")
        write (*, wrtfmt1e) ' average[u(x,y,z)]_[x,y,z]: ', ubulk
    end if
    !----------------------------------------------------------------------------------------------------------
    ! to check maximum velocity
    !----------------------------------------------------------------------------------------------------------
    call Find_maximum_absvar3d(fl%qx, fl%umax(1), "maximum ux:", wrtfmt1e)
    call Find_maximum_absvar3d(fl%qy, fl%umax(2), "maximum uy:", wrtfmt1e)
    call Find_maximum_absvar3d(fl%qz, fl%umax(3), "maximum uz:", wrtfmt1e)
    !----------------------------------------------------------------------------------------------------------
    ! to set up other parameters for flow only, which will be updated in thermo flow.
    !----------------------------------------------------------------------------------------------------------
    fl%pcor(:, :, :) = ZERO
    fl%pcor_zpencil_ggg(:, :, :) = ZERO
    if(dm%is_thermo) then
      fl%dDens  (:, :, :) = ONE
      fl%mVisc  (:, :, :) = ONE
      fl%dDensm1(:, :, :) = ONE
      fl%dDensm2(:, :, :) = ONE
    end if

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine read_instantanous_thermo(tm, dm)
    use udf_type_mod
    use thermo_info_mod
    use decomp_2d_io
    use files_io_mod
    use io_tools_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    type(t_thermo), intent(inout) :: tm

    character(120):: data_flname
    character(120):: keyword

    if (.not. dm%is_thermo) return
    if(nrank == 0) call Print_debug_start_msg("read instantanous thermo data ... ...")

    tm%iteration = tm%iterfrom

    keyword = 'rhoh'
    call generate_file_name(data_flname, dm%idom, keyword, 'bin', tm%iteration)
    call decomp_2d_read_one(X_PENCIL, tm%dh, trim(dir_data), trim(data_flname), io_name, dm%dccc, reduce_prec=.false.)

    keyword = 'temp'
    call generate_file_name(data_flname, dm%idom, keyword, 'bin', tm%iteration)
    call decomp_2d_read_one(X_PENCIL, tm%tTemp, trim(dir_data), trim(data_flname), io_name, dm%dccc, reduce_prec=.false.)
    
    tm%time = real(tm%iterfrom, WP) * dm%dt 

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
!==========================================================================================================
  subroutine restore_thermo_variables_from_restart(fl, tm, dm)
    use udf_type_mod
    use thermo_info_mod
    use eq_energy_mod
    use solver_tools_mod
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm

    if (.not. dm%is_thermo) return

    call Update_thermal_properties(fl, tm, dm)
    call Calculate_massflux_from_velocity (fl, dm)

    fl%dDensm1(:, :, :) = fl%dDens(:, :, :)
    fl%dDensm2(:, :, :) = fl%dDens(:, :, :)

    return
  end subroutine

end module 