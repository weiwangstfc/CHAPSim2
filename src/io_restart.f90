module io_restart_mod
  use print_msg_mod
  use parameters_constant_mod
  use decomp_2d_io
  use udf_type_mod
  use io_files_mod
  use io_tools_mod
  implicit none 

  character(len=10), parameter :: io_name = "restart-io"

  private :: write_instantanous_array
  public  :: write_instantanous_flow
  public  :: write_instantanous_thermo
  private :: read_instantanous_array
  public  :: read_instantanous_flow
  public  :: read_instantanous_thermo
  public  :: restore_flow_variables_from_restart
  public  :: restore_thermo_variables_from_restart

  private :: append_instantanous_xoutlet
  private :: write_instantanous_plane
  public  :: write_instantanous_xoutlet

  private :: assign_instantanous_xinlet
  private :: read_instantanous_plane
  public  :: read_instantanous_xinlet

contains 
!==========================================================================================================
!==========================================================================================================
  subroutine read_instantanous_array(var, keyword, idom, iter, dtmp)
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
    implicit none 
    real(WP), intent(in) :: var( :, :, :)
    type(DECOMP_INFO), intent(in) :: dtmp
    character(*), intent(in) :: keyword
    integer, intent(in) :: idom
    integer, intent(in) :: iter
    
    character(120):: data_flname_path

    call generate_pathfile_name(data_flname_path, idom, trim(keyword), dir_data, 'bin', iter)
    call decomp_2d_write_one(X_PENCIL, var, trim(data_flname_path), dtmp)

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine write_instantanous_flow(fl, dm)
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(in) :: fl

    character(120):: data_flname_path
    character(120):: keyword

    if(nrank == 0) call Print_debug_start_msg("writing out instantanous 3d flow data ...")

    call write_instantanous_array(fl%qx, 'qx', dm%idom, fl%iteration, dm%dpcc)
    call write_instantanous_array(fl%qy, 'qy', dm%idom, fl%iteration, dm%dcpc)
    call write_instantanous_array(fl%qz, 'qz', dm%idom, fl%iteration, dm%dccp)
    call write_instantanous_array(fl%pres, 'pr', dm%idom, fl%iteration, dm%dccc)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine write_instantanous_thermo(tm, dm)
    use thermo_info_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(in) :: tm

    character(120):: data_flname_path
    character(120):: keyword
    

    if(nrank == 0) call Print_debug_start_msg("writing out instantanous 3d thermo data ...")

    call write_instantanous_array(tm%rhoh,    'rhoh', dm%idom, tm%iteration, dm%dccc)
    call write_instantanous_array(tm%tTemp, 'temp', dm%idom, tm%iteration, dm%dccc)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine read_instantanous_flow(fl, dm)
    implicit none
    type(t_domain), intent(inout) :: dm
    type(t_flow),   intent(inout) :: fl

    character(120):: data_flname
    character(120):: keyword


    if(nrank == 0) call Print_debug_start_msg("read instantanous flow data ... ...")

    call read_instantanous_array(fl%qx, 'qx', dm%idom, fl%iterfrom, dm%dpcc)
    call read_instantanous_array(fl%qy, 'qy', dm%idom, fl%iterfrom, dm%dcpc)
    call read_instantanous_array(fl%qz, 'qz', dm%idom, fl%iterfrom, dm%dccp)
    call read_instantanous_array(fl%pres, 'pr', dm%idom, fl%iterfrom, dm%dccc)
    
    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine restore_flow_variables_from_restart(fl, dm)
    use mpi_mod
    use boundary_conditions_mod
    use solver_tools_mod
    use wtformat_mod
    use find_max_min_ave_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(inout) :: fl
    real(WP) :: ubulk
    

    !call Get_volumetric_average_3d(.false., dm%ibcy_qx(:), dm%fbcy_qx(:, :, :), dm, dm%dpcc, fl%qx, ubulk, "ux")
    call Get_volumetric_average_3d_for_var_xcx(dm, dm%dpcc, fl%qx, ubulk, LF3D_VOL_AVE, "ux")
    if(nrank == 0) then
        Call Print_debug_mid_msg("  The restarted mass flux is:")
        write (*, wrtfmt1e) ' average[u(x,y,z)]_[x,y,z]: ', ubulk
    end if
    !----------------------------------------------------------------------------------------------------------
    ! to check maximum velocity
    !----------------------------------------------------------------------------------------------------------
    call Find_max_min_3d(fl%qx, "qx: ", wrtfmt2e)
    call Find_max_min_3d(fl%qy, "qy: ", wrtfmt2e)
    call Find_max_min_3d(fl%qz, "qz: ", wrtfmt2e)
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
    use thermo_info_mod
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
    call decomp_2d_read_one(X_PENCIL, tm%rhoh, trim(dir_data), trim(data_flname), io_name, dm%dccc, reduce_prec=.false.)

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
    use convert_primary_conservative_mod
    type(t_domain), intent(inout) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm

    if (.not. dm%is_thermo) return

    call Update_thermal_properties(fl, tm, dm)
    call convert_primary_conservative (fl, dm, IQ2G)

    fl%dDensm1(:, :, :) = fl%dDens(:, :, :)
    fl%dDensm2(:, :, :) = fl%dDens(:, :, :)

    return
  end subroutine


!==========================================================================================================
  subroutine append_instantanous_xoutlet(fl, dm)
    implicit none 
    type(t_flow), intent(in) :: fl
    type(t_domain), intent(inout) :: dm

    integer :: niter, j, k
    type(DECOMP_INFO) :: dtmp

    ! based on x pencil
    if(.not. dm%is_record_xoutlet) return

    ! if dm%ndbfre = 10
    ! store : file_10, store  1,  2, ..., 10
    !         file_20, store 11, 12, ..., 20
    niter = mod(fl%iteration, dm%ndbfre) !
    if(niter == 0) niter =  dm%ndbfre

    dtmp = dm%dpcc
    do j = 1, dtmp%xsz(2)
      do k = 1, dtmp%xsz(3)
        dm%fbcx_qx_outl1(niter, j, k) = fl%qx(dtmp%xsz(1),   j, k)
        dm%fbcx_qx_outl2(niter, j, k) = fl%qx(dtmp%xsz(1)-1, j, k)
      end do
    end do

    !write(*, *) 'j, fl%qx(1, j, 1), dm%fbcx_qx_outl1(niter, j, 1)'
    ! do j = 1, dm%dpcc%xsz(2)
    !   write(*, *) j, fl%qx(dtmp%xsz(1), j, 1), dm%fbcx_qx_outl1(niter, j, 1)
    ! end do

    dtmp = dm%dcpc
    do j = 1, dtmp%xsz(2)
      do k = 1, dtmp%xsz(3)
        dm%fbcx_qy_outl1(niter, j, k) = fl%qy(dtmp%xsz(1),   j, k)
        dm%fbcx_qy_outl2(niter, j, k) = fl%qy(dtmp%xsz(1)-1, j, k)
      end do
    end do

    dtmp = dm%dccp
    do j = 1, dtmp%xsz(2)
      do k = 1, dtmp%xsz(3)
        dm%fbcx_qz_outl1(niter, j, k) = fl%qz(dtmp%xsz(1),   j, k)
        dm%fbcx_qz_outl2(niter, j, k) = fl%qz(dtmp%xsz(1)-1, j, k)
      end do
    end do

    dtmp = dm%dccc
    do j = 1, dtmp%xsz(2)
      do k = 1, dtmp%xsz(3)
        dm%fbcx_pr_outl1(niter, j, k) = fl%pres(dtmp%xsz(1),   j, k)
        dm%fbcx_pr_outl2(niter, j, k) = fl%pres(dtmp%xsz(1)-1, j, k)
      end do
    end do

    return
  end subroutine
!==========================================================================================================
  subroutine write_instantanous_plane(var, keyword, idom, iter, niter, dtmp)
    implicit none 
    real(WP), intent(in) :: var( :, :, :)
    type(DECOMP_INFO), intent(in) :: dtmp
    character(*), intent(in) :: keyword
    integer, intent(in) :: idom
    integer, intent(in) :: iter, niter

    character(120):: data_flname_path

    call generate_pathfile_name(data_flname_path, idom, trim(keyword), dir_data, 'bin', iter)

    if(nrank==0) write(*, *) 'Write outflow data to '//trim(data_flname_path)
 
    call decomp_2d_open_io (io_in2outlet, trim(data_flname_path), decomp_2d_write_mode)
    call decomp_2d_start_io(io_in2outlet, trim(data_flname_path))!

    call decomp_2d_write_outflow(trim(data_flname_path), trim(keyword), niter, var, io_in2outlet, dtmp)
    !call decomp_2d_write_plane(X_PENCIL, var, 1, dtmp%xsz(1), trim(data_flname_path), dtmp)
    !write(*,*)var

    call decomp_2d_end_io(io_in2outlet, trim(data_flname_path))
    call decomp_2d_close_io(io_in2outlet, trim(data_flname_path))

    return
  end subroutine
!==========================================================================================================
  subroutine write_instantanous_xoutlet(fl, dm)
    implicit none 
    type(t_flow), intent(in) :: fl
    type(t_domain), intent(inout) :: dm
    
    character(120):: data_flname_path
    integer :: idom, iter


    if(.not. dm%is_record_xoutlet) return

    call append_instantanous_xoutlet(fl, dm)

    if(mod(fl%iteration, dm%ndbfre) /= 0) return
    call write_instantanous_plane(dm%fbcx_qx_outl1, 'outlet1_qx', dm%idom, fl%iteration, dm%ndbfre, dm%dxcc)
    call write_instantanous_plane(dm%fbcx_qx_outl2, 'outlet2_qx', dm%idom, fl%iteration, dm%ndbfre, dm%dxcc)
    call write_instantanous_plane(dm%fbcx_qy_outl1, 'outlet1_qy', dm%idom, fl%iteration, dm%ndbfre, dm%dxpc)
    call write_instantanous_plane(dm%fbcx_qy_outl2, 'outlet2_qy', dm%idom, fl%iteration, dm%ndbfre, dm%dxpc)
    call write_instantanous_plane(dm%fbcx_qz_outl1, 'outlet1_qz', dm%idom, fl%iteration, dm%ndbfre, dm%dxcp)
    call write_instantanous_plane(dm%fbcx_qz_outl2, 'outlet2_qz', dm%idom, fl%iteration, dm%ndbfre, dm%dxcp)
    call write_instantanous_plane(dm%fbcx_pr_outl1, 'outlet1_pr', dm%idom, fl%iteration, dm%ndbfre, dm%dxcc)
    call write_instantanous_plane(dm%fbcx_pr_outl2, 'outlet2_pr', dm%idom, fl%iteration, dm%ndbfre, dm%dxcc)

    return
  end subroutine
!==========================================================================================================
  subroutine assign_instantanous_xinlet(fl, dm)
    implicit none 
    type(t_flow), intent(in) :: fl
    type(t_domain), intent(inout) :: dm

    integer :: iter, j, k
    type(DECOMP_INFO) :: dtmp

    ! based on x pencil
    if(.not. dm%is_read_xinlet) return

    if(fl%iteration > dm%ndbend) then
      iter = mod(fl%iteration, dm%ndbend) ! database recycle
    else if (fl%iteration == 0) then
      iter = 1
    else
      iter = fl%iteration
    end if

    iter = mod(iter, dm%ndbfre)
    if(iter == 0) iter = dm%ndbfre

    if(dm%ibcx_nominal(1, 1) == IBC_DATABASE) then
      dtmp = dm%dpcc
      do j = 1, dtmp%xsz(2)
        do k = 1, dtmp%xsz(3)
          dm%fbcx_qx(1, j, k) = dm%fbcx_qx_inl1(iter, j, k)
          dm%fbcx_qx(3, j, k) = dm%fbcx_qx_inl2(iter, j, k)
          ! check, below 
          !fl%qx(1, j, k) = dm%fbcx_qx(1, j, k)
        end do
      end do
      !if(nrank == 0) write(*,*) 'fbcx_qx = ', iter, dm%fbcx_qx(1, :, :)
    end if

        ! test
    !write(*,*) 'j, fl%qx(1, j, 1), dm%fbcx_qx(1, j, 1)'
    !do j = 1, dm%dpcc%xsz(2)
      !write(*,*) j, fl%qx(1, j, 1), dm%fbcx_qx(1, j, 1)
    !end do

    if(dm%ibcx_nominal(1, 2) == IBC_DATABASE) then
      dtmp = dm%dcpc
      do j = 1, dtmp%xsz(2)
        do k = 1, dtmp%xsz(3)
          dm%fbcx_qy(1, j, k) = dm%fbcx_qy_inl1(iter, j, k)
          dm%fbcx_qy(3, j, k) = dm%fbcx_qy_inl2(iter, j, k)
        end do
      end do
      !if(nrank == 0) write(*,*) 'fbcx_qy = ', iter, dm%fbcx_qy(1, :, :)
    end if

    if(dm%ibcx_nominal(1, 3) == IBC_DATABASE) then
      dtmp = dm%dccp
      do j = 1, dtmp%xsz(2)
        do k = 1, dtmp%xsz(3)
          dm%fbcx_qz(1, j, k) = dm%fbcx_qz_inl1(iter, j, k)
          dm%fbcx_qz(3, j, k) = dm%fbcx_qz_inl2(iter, j, k)
        end do
      end do
      !if(nrank == 0) write(*,*) 'fbcx_qz = ', iter, dm%fbcx_qz(1, :, :)
    end if

    if(dm%ibcx_nominal(1, 4) == IBC_DATABASE) then
      dtmp = dm%dccc
      do j = 1, dtmp%xsz(2)
        do k = 1, dtmp%xsz(3)
          dm%fbcx_pr(1, j, k) = dm%fbcx_pr_inl1(iter, j, k)
          dm%fbcx_pr(3, j, k) = dm%fbcx_pr_inl2(iter, j, k)
        end do
      end do
      !if(nrank == 0) write(*,*) 'fbcx_pr = ', iter, dm%fbcx_pr(1, :, :)
    end if

    return
  end subroutine
!==========================================================================================================
  subroutine read_instantanous_plane(var, keyword, idom, iter, nfre, dtmp)
    use decomp_2d_io
    implicit none 
    real(WP), intent(inout) :: var( :, :, :)
    type(DECOMP_INFO), intent(in) :: dtmp
    character(*), intent(in) :: keyword
    integer, intent(in) :: idom
    integer, intent(in) :: iter
    integer, intent(in) :: nfre

    character(120):: data_flname_path

    call generate_pathfile_name(data_flname_path, idom, trim(keyword), dir_data, 'bin', iter)

    call decomp_2d_open_io (io_in2outlet, trim(data_flname_path), decomp_2d_read_mode)
    if(nrank == 0) call Print_debug_mid_msg("Read data on a plane from file: "//trim(data_flname_path))
    call decomp_2d_read_inflow(trim(data_flname_path), trim(keyword), nfre, var, io_in2outlet, dtmp)
    !write(*,*) var
    call decomp_2d_close_io(io_in2outlet, trim(data_flname_path))

    return
  end subroutine
!==========================================================================================================
  subroutine read_instantanous_xinlet(fl, dm)
    use typeconvert_mod
    implicit none 
    type(t_flow), intent(in) :: fl
    type(t_domain), intent(inout) :: dm
    
    character(120):: data_flname_path
    integer :: idom, iter, niter


    if(.not. dm%is_read_xinlet) return

    if(fl%iteration > dm%ndbend) then
      iter = mod(fl%iteration, dm%ndbend) ! database recycle
    else
      iter = fl%iteration
    end if

    niter = (iter + dm%ndbfre - 1) / dm%ndbfre ! integer operation
    niter = niter * dm%ndbfre

    if(fl%iteration == 0) niter = dm%ndbfre

    if(mod(iter - 1, dm%ndbfre) == 0 .or. &
       fl%iteration == 0) then
      if(nrank == 0) call Print_debug_mid_msg("Read inlet database at "&
        //trim(int2str(iter))//' in '//trim(int2str(niter)))

      call read_instantanous_plane(dm%fbcx_qx_inl1, 'outlet1_qx', dm%idom, niter, dm%ndbfre, dm%dxcc)
      call read_instantanous_plane(dm%fbcx_qx_inl2, 'outlet2_qx', dm%idom, niter, dm%ndbfre, dm%dxcc)
      call read_instantanous_plane(dm%fbcx_qy_inl1, 'outlet1_qy', dm%idom, niter, dm%ndbfre, dm%dxpc)
      call read_instantanous_plane(dm%fbcx_qy_inl2, 'outlet2_qy', dm%idom, niter, dm%ndbfre, dm%dxpc)
      call read_instantanous_plane(dm%fbcx_qz_inl1, 'outlet1_qz', dm%idom, niter, dm%ndbfre, dm%dxcp)
      call read_instantanous_plane(dm%fbcx_qz_inl2, 'outlet2_qz', dm%idom, niter, dm%ndbfre, dm%dxcp)
      call read_instantanous_plane(dm%fbcx_pr_inl1, 'outlet1_pr', dm%idom, niter, dm%ndbfre, dm%dxcc)
      call read_instantanous_plane(dm%fbcx_pr_inl2, 'outlet2_pr', dm%idom, niter, dm%ndbfre, dm%dxcc)
    end if

    call assign_instantanous_xinlet(fl, dm) ! every iteration



    return
  end subroutine
!==========================================================================================================
end module 