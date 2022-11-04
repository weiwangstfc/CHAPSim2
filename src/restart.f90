module restart_mod
  implicit none 

  public :: write_instantanous_flow_data
  public :: write_instantanous_thermo_data
  public :: read_instantanous_flow_raw_data
  public :: restore_flow_variables_from_restart

contains 
!==========================================================================================================
!==========================================================================================================
  subroutine write_instantanous_flow_data(fl, dm)
    use udf_type_mod
    use mpi_mod
    use parameters_constant_mod
    use typeconvert_mod
    use decomp_2d_io
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(in) :: fl

    integer :: ipencil
    character( len = 128) :: filename
    logical :: file_exists = .FALSE.
    integer :: outputunit
    

    if(nrank == 0) then
      call Print_debug_start_msg("write instantanous flow data ...")

      filename = 'results_flow.case'
      INQUIRE(FILE = trim(filename), exist = file_exists)
      if(.not.file_exists) then
        open(newunit = outputunit, file = trim(filename), action = "write", status = "new")
        write(outputunit, '(A, 1I4.1)')  '# Case ID:                  ', dm%icase
        write(outputunit, '(A, 2F10.4)') '# Computational Domain - x: ', zero, dm%lxx
        write(outputunit, '(A, 2F10.4)') '# Computational Domain - y: ', dm%lyb, dm%lyt
        write(outputunit, '(A, 2F10.4)') '# Computational Domain - z: ', zero, dm%lzz
        write(outputunit, '(A, 3I10.1)') '# Computational cell size:  ', dm%nc(1:3)
        write(outputunit, '(A, 3I10.1)') '# Write flow raw data at iteration (1st col) and time(2nd col)'
      else
        open(newunit = outputunit, file = trim(filename), action = "write", status = "old", position="append")
      end if
      write(outputunit, '(1I10.1, F11.5)') fl%iteration, fl%time
      close(outputunit)

    end if

    ipencil = 1
    call decomp_2d_write_one(ipencil, fl%qx,   'instantanous_ux_'//trim(int2str(fl%iteration))//'.dat', dm%dpcc)
    call decomp_2d_write_one(ipencil, fl%qy,   'instantanous_uy_'//trim(int2str(fl%iteration))//'.dat', dm%dcpc)
    call decomp_2d_write_one(ipencil, fl%qz,   'instantanous_uz_'//trim(int2str(fl%iteration))//'.dat', dm%dccp)
    call decomp_2d_write_one(ipencil, fl%pres, 'instantanous_pr_'//trim(int2str(fl%iteration))//'.dat', dm%dccc)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine write_instantanous_thermo_data(tm, dm)
    use udf_type_mod
    use thermo_info_mod
    use mpi_mod
    use typeconvert_mod
    use decomp_2d_io
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(in) :: tm

    integer :: ipencil
    character( len = 128) :: filename
    logical :: file_exists = .FALSE.
    integer :: outputunit
    

    if(nrank == 0) then
      call Print_debug_start_msg("write instantanous thermal data ...")

      filename = 'results_thermo.case'
      INQUIRE(FILE = trim(filename), exist = file_exists)
      if(.not.file_exists) then
        open(newunit = outputunit, file = trim(filename), action = "write", status = "new")
        write(outputunit, '(A, 1I4.1)')  '# Case ID:                  ', dm%icase
        write(outputunit, '(A, 2F10.4)') '# Computational Domain - x: ', zero, dm%lxx
        write(outputunit, '(A, 2F10.4)') '# Computational Domain - y: ', dm%lyb, dm%lyt
        write(outputunit, '(A, 2F10.4)') '# Computational Domain - z: ', zero, dm%lzz
        write(outputunit, '(A, 3I10.1)') '# Computational cell size:  ', dm%nc(1:3)
        write(outputunit, '(A, 3I10.1)') '# Write thermo raw data at iteration (1st col) and time(2nd col)'
      else
        open(newunit = outputunit, file = trim(filename), action = "write", status = "old", position="append")
      end if
      write(outputunit, '(A, 1I10.1, F11.5)') tm%iteration, tm%time
      close(outputunit)

    end if

    ipencil = 1
    call decomp_2d_write_one(ipencil, tm%dh,    'instantanous_dh_'//trim(int2str(tm%iteration))//'.dat', dm%dccc)
    call decomp_2d_write_one(ipencil, tm%tTemp, 'instantanous_te_'//trim(int2str(tm%iteration))//'.dat', dm%dccc)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine read_instantanous_flow_raw_data(fl, dm)
    use udf_type_mod
    use mpi_mod
    use typeconvert_mod
    use decomp_2d_io
    implicit none
    type(t_domain), intent(in)    :: dm
    type(t_flow),   intent(inout) :: fl

    integer :: ipencil
    character( len = 128) :: filename
    logical :: file_exists = .FALSE.
    integer :: inputunit
    integer :: ioerr
    character( len = 128 ) :: strinfo
    integer :: itmp
    real(WP) :: rtmp


    if(nrank == 0) then
      call Print_debug_start_msg("read instantanous flow data ... ...")

      filename = 'results_flow.case'
      INQUIRE(FILE = trim(filename), exist = file_exists)
      if(.not.file_exists) then
        call Print_error_msg ("Cannot finding "//trim(filename) )
      else
        open(newunit = inputunit, file = trim(filename), action = "read", status = "old")
      end if

      read(inputunit, *, iostat = ioerr) strinfo, rtmp, rtmp
      read(inputunit, *, iostat = ioerr) strinfo, rtmp, rtmp
      read(inputunit, *, iostat = ioerr) strinfo, rtmp, rtmp
      read(inputunit, *, iostat = ioerr) strinfo, rtmp, rtmp
      read(inputunit, *, iostat = ioerr) strinfo, itmp, itmp, itmp
      read(inputunit, *, iostat = ioerr) strinfo

      do while (ioerr == 0)
        read(inputunit, *, iostat = ioerr) fl%iteration, fl%time
      end do

      close(inputunit)

      if(fl%iteration /= fl%nrsttckpt) &
      call Print_error_msg ("Flow restart file does not give correct iteration number or physical time.")

    end if

    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_bcast(fl%iteration, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror)
    call mpi_bcast(fl%time,      1, MPI_REAL_WP, 0, MPI_COMM_WORLD, ierror)
    call mpi_barrier(MPI_COMM_WORLD, ierror)

    ipencil = 1
    call decomp_2d_read_one(ipencil, fl%qx,   'instantanous_ux_'//trim(int2str(fl%nrsttckpt))//'.dat', dm%dpcc)
    call decomp_2d_read_one(ipencil, fl%qy,   'instantanous_uy_'//trim(int2str(fl%nrsttckpt))//'.dat', dm%dcpc)
    call decomp_2d_read_one(ipencil, fl%qz,   'instantanous_uz_'//trim(int2str(fl%nrsttckpt))//'.dat', dm%dccp)
    call decomp_2d_read_one(ipencil, fl%pres, 'instantanous_pr_'//trim(int2str(fl%nrsttckpt))//'.dat', dm%dccc)

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
    
    call Apply_BC_velocity(dm, fl%qx, fl%qx, fl%qx)
    call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), dm%fbcy(:, 1), dm, dm%dpcc, fl%qx, ubulk)
    if(nrank == 0) then
        Call Print_debug_mid_msg("  The restarted mass flux is:")
        write (*, wrtfmt1r) ' average[u(x,y,z)]_[x,y,z]: ', ubulk
    end if
    !----------------------------------------------------------------------------------------------------------
    ! to check maximum velocity
    !----------------------------------------------------------------------------------------------------------
    call Find_maximum_absvar3d(fl%qx, "maximum ux:", wrtfmt1e)
    call Find_maximum_absvar3d(fl%qy, "maximum uy:", wrtfmt1e)
    call Find_maximum_absvar3d(fl%qz, "maximum uz:", wrtfmt1e)
    !----------------------------------------------------------------------------------------------------------
    ! to set up other parameters for flow only, which will be updated in thermo flow.
    !----------------------------------------------------------------------------------------------------------
    fl%pcor(:, :, :) = ZERO
    fl%pcor_zpencil_ggg(:, :, :) = ZERO
    fl%dDens  (:, :, :) = ONE
    fl%mVisc  (:, :, :) = ONE
    fl%dDensm1(:, :, :) = ONE
    fl%dDensm2(:, :, :) = ONE

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine read_instantanous_thermo_raw_data(tm, dm)
    use udf_type_mod
    use thermo_info_mod
    use mpi_mod
    use typeconvert_mod
    use decomp_2d_io
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(inout) :: tm

    integer :: ipencil
    character( len = 128) :: filename
    logical :: file_exists = .FALSE.
    integer :: inputunit
    integer :: ioerr
    character( len = 128 ) :: strinfo
    integer :: itmp
    real(WP) :: rtmp

    if (dm%ithermo /= 1) return

    if(nrank == 0) then
      call Print_debug_start_msg("read instantanous thermo data ... ...")

      filename = 'results_thermo.case'
      INQUIRE(FILE = trim(filename), exist = file_exists)
      if(.not.file_exists) then
        call Print_error_msg ("Cannot finding "//trim(filename) )
      else
        open(newunit = inputunit, file = trim(filename), action = "read", status = "old")
      end if

      read(inputunit, *, iostat = ioerr) strinfo, rtmp, rtmp
      read(inputunit, *, iostat = ioerr) strinfo, rtmp, rtmp
      read(inputunit, *, iostat = ioerr) strinfo, rtmp, rtmp
      read(inputunit, *, iostat = ioerr) strinfo, rtmp, rtmp
      read(inputunit, *, iostat = ioerr) strinfo, itmp, itmp, itmp
      read(inputunit, *, iostat = ioerr) strinfo

      do while (ioerr == 0)
        read(inputunit, *, iostat = ioerr) tm%iteration, tm%time
      end do

      close(inputunit)

      if(tm%iteration /= tm%nrsttckpt) &
      call Print_error_msg ("Thermo restart file does not give correct iteration number or physical time.")

    end if

    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_bcast(tm%iteration, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror)
    call mpi_bcast(tm%time,      1, MPI_REAL_WP, 0, MPI_COMM_WORLD, ierror)
    call mpi_barrier(MPI_COMM_WORLD, ierror)


    ipencil = 1
    call decomp_2d_read_one(ipencil, tm%dh,    'instantanous_dh_'//trim(int2str(tm%nrsttckpt))//'.dat', dm%dccc)
    call decomp_2d_read_one(ipencil, tm%tTemp, 'instantanous_te_'//trim(int2str(tm%nrsttckpt))//'.dat', dm%dccc)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
!==========================================================================================================
  subroutine restore_thermo_variables_from_restart(fl, tm, dm)
    use udf_type_mod
    use thermo_info_mod
    use mpi_mod
    use solver_tools_mod
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm

    if (dm%ithermo /= 1) return

    call Update_thermal_properties(fl, tm, dm)
    call Calculate_massflux_from_velocity (fl, dm)

    fl%dDensm1(:, :, :) = fl%dDens(:, :, :)
    fl%dDensm2(:, :, :) = fl%dDens(:, :, :)

    return
  end subroutine

end module 