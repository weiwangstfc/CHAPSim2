module io_monitor_mod
  use precision_mod
  implicit none

  public :: write_monitor_ini
  public :: write_monitor_total
  public :: write_monitor_flow
  public :: write_monitor_thermo
  
contains
  subroutine write_monitor_ini(dm)
    use typeconvert_mod
    use wtformat_mod
    use udf_type_mod
    use files_io_mod
    use io_tools_mod
    use parameters_constant_mod
    implicit none 
    type(t_domain),  intent(inout) :: dm

    integer :: myunit
    integer :: i, j
    logical :: exist
    character(len=120) :: flname
    character(len=120) :: keyword

    integer :: idgb(3)
    integer :: nplc
    logical :: is_y, is_z
    integer, allocatable :: probeid(:, :)

    
    allocate( dm%probe_is_in(dm%proben) )
    dm%probe_is_in(:) = .false.

    allocate( probeid(3, dm%proben) )
    nplc = 0
    do i = 1, dm%proben
!----------------------------------------------------------------------------------------------------------
! probe points find the nearest cell centre, global index info, then convert to local index in x-pencil
!----------------------------------------------------------------------------------------------------------
      idgb(1:3) = 0

      idgb(1) = ceiling ( dm%probexyz(1, i) / dm%h(1) )
      idgb(3) = ceiling ( dm%probexyz(3, i) / dm%h(3) )

      if(dm%is_periodic(2)) then 
        if( dm%probexyz(2, i) >= dm%yp(dm%np(2)) .and. dm%probexyz(2, i) < dm%lyt) idgb(2) = dm%nc(2)
      end if
      do j = 1, dm%np(2) - 1
        if (dm%probexyz(2, i) >= dm%yp(j) .and. &
            dm%probexyz(2, i) < dm%yp(j+1)) then
          idgb(2) = j
        end if
      end do
!----------------------------------------------------------------------------------------------------------
! convert global id to local, based on x-pencil
!----------------------------------------------------------------------------------------------------------
      is_y = .false.
      is_z = .false.
      if( idgb(2) >= dm%dccc%xst(2) .and. idgb(2) <= dm%dccc%xen(2) ) is_y = .true.
      if( idgb(3) >= dm%dccc%xst(3) .and. idgb(3) <= dm%dccc%xen(3) ) is_z = .true.
      if(is_y .and. is_z) then 
        dm%probe_is_in(i) = .true.
        nplc = nplc + 1
        probeid(1, nplc) = idgb(1)
        probeid(2, nplc) = idgb(2) - dm%dccc%xst(2) + 1
        probeid(3, nplc) = idgb(3) - dm%dccc%xst(3) + 1
      end if
    end do

    if(nplc > 0) allocate(dm%probexid(3, nplc))

    do i = 1, nplc 
      dm%probexid(1:3, i) = probeid(1:3, i)
    end do

    deallocate (probeid)
!----------------------------------------------------------------------------------------------------------
! create probe history file for flow
!----------------------------------------------------------------------------------------------------------
    do i = 1, dm%proben
      if(.not. dm%probe_is_in(i)) cycle 

      keyword = "monitor_pt"//trim(int2str(i))//"_flow"
      call generate_pathfile_name(flname, dm%idom, keyword, dir_moni, 'dat')

      inquire(file = trim(flname), exist = exist)
      if (exist) then
        !open(newunit = myunit, file = trim(flname), status="old", position="append", action="write")
      else
        open(newunit = myunit, file = trim(flname), status="new", action="write")
        write(myunit, *) "# domain-id : ", dm%idom, "pt-id : ", i
        write(myunit, *) "# probe pts location ",  dm%probexyz(1:3, i)
        write(myunit, *) "# t, u, v, w" ! to add more instantanous or statistics
        close(myunit)
      end if
    end do
!----------------------------------------------------------------------------------------------------------
! create probe history file for thermo
!----------------------------------------------------------------------------------------------------------
    if(dm%is_thermo) then
      do i = 1, dm%proben
        if(.not. dm%probe_is_in(i)) cycle 
  
        keyword = "monitor_pt"//trim(int2str(i))//"_thermo"
        call generate_pathfile_name(flname, dm%idom, keyword, dir_moni, 'dat')
  
        inquire(file = trim(flname), exist = exist)
        if (exist) then
          !open(newunit = myunit, file = trim(flname), status="old", position="append", action="write")
        else
          open(newunit = myunit, file = trim(flname), status="new", action="write")
          write(myunit, *) "# domain-id : ", dm%idom, "pt-id : ", i
          write(myunit, *) "# probe pts location ",  dm%probexyz(1:3, i)
          write(myunit, *) "# t, T" ! to add more instantanous or statistics
          close(myunit)
        end if
      end do
    end if

!----------------------------------------------------------------------------------------------------------
! create history file for total variables
!----------------------------------------------------------------------------------------------------------
    keyword = "monitor_bulk"
    call generate_pathfile_name(flname, dm%idom, keyword, dir_moni, 'dat')
    inquire(file = trim(flname), exist = exist)
    if (exist) then
      !open(newunit = myunit, file = trim(flname), status="old", position="append", action="write")
    else
      open(newunit = myunit, file = trim(flname), status="new", action="write")
      write(myunit, *) "# domain-id : ", dm%idom, "pt-id : ", i
      write(myunit, *) "# energy" ! to add more instantanous or statistics
      close(myunit)
    end if


    return
  end subroutine
!==========================================================================================================
  subroutine write_monitor_total(fl, dm)
    use typeconvert_mod
    use parameters_constant_mod
    use wtformat_mod
    use udf_type_mod
    use files_io_mod
    use io_tools_mod
    use solver_tools_mod
    use operations
    implicit none 

    type(t_domain),  intent(in) :: dm
    type(t_flow), intent(in) :: fl

    character(len=120) :: flname
    character(len=120) :: keyword
    character(200) :: iotxt
    integer :: ioerr, myunit

    real(WP) :: bulk_energy
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc1
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc2
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc3
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: fenergy
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil

!----------------------------------------------------------------------------------------------------------
! open file
!----------------------------------------------------------------------------------------------------------
    keyword = "monitor_bulk"
    call generate_pathfile_name(flname, dm%idom, keyword, dir_moni, 'dat')

    open(newunit = myunit, file = trim(flname), status = "old", action = "write", position = "append", &
        iostat = ioerr, iomsg = iotxt)
    if(ioerr /= 0) then
      write (*, *) 'Problem openning probing file'
      write (*, *) 'Message: ', trim (iotxt)
      stop
    end if       
!----------------------------------------------------------------------------------------------------------
!   ux
!----------------------------------------------------------------------------------------------------------
    call Get_x_midp_P2C_3D(fl%qx, accc1, dm, dm%ibcx(:, 1))
!----------------------------------------------------------------------------------------------------------
!   uy
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qy, acpc_ypencil, dm%dcpc)
    call Get_y_midp_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%ibcy(:, 2))
    call transpose_y_to_x(accc_ypencil, accc2, dm%dccc)
!----------------------------------------------------------------------------------------------------------
!   uz
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qz, accp_ypencil, dm%dccp)
    call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
    call Get_z_midp_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%ibcz(:, 3))
    call transpose_z_to_y(accc_zpencil, accc_ypencil, dm%dccc)
    call transpose_y_to_x(accc_ypencil, accc3, dm%dccc)
!----------------------------------------------------------------------------------------------------------
!   x-pencil, 1/2*(uu+vv+ww) - calculation
!----------------------------------------------------------------------------------------------------------
    fenergy = HALF * (accc1 * accc1 + accc2 * accc2 + accc3 * accc3)
    call Get_volumetric_average_3d(.false., dm%ibcy(:, 1), dm%fbcy(:, 1), dm, dm%dccc, fenergy, bulk_energy)
!----------------------------------------------------------------------------------------------------------
!   write data out
!----------------------------------------------------------------------------------------------------------
    write(myunit, *) fl%time, bulk_energy
    close(myunit)

    return
  end subroutine

!==========================================================================================================
  subroutine write_monitor_flow(fl, dm)
    use typeconvert_mod
    use parameters_constant_mod
    use wtformat_mod
    use udf_type_mod
    use files_io_mod
    use io_tools_mod
    implicit none 

    type(t_domain),  intent(in) :: dm
    type(t_flow), intent(in) :: fl

    character(len=120) :: flname
    character(len=120) :: keyword
    character(200) :: iotxt
    integer :: ioerr, myunit
    integer :: ix, iy, iz
    integer :: i, j
!----------------------------------------------------------------------------------------------------------
! based on x-pencil
!----------------------------------------------------------------------------------------------------------
    j = 0
    do i = 1, dm%proben
      if( dm%probe_is_in(i) ) j = j + 1
      if(j > 0) then
!----------------------------------------------------------------------------------------------------------
! open file
!----------------------------------------------------------------------------------------------------------
        keyword = "monitor_pt"//trim(int2str(i))//"_flow"
        call generate_pathfile_name(flname, dm%idom, keyword, dir_moni, 'dat')

        open(newunit = myunit, file = trim(flname), status = "old", action = "write", position = "append", &
            iostat = ioerr, iomsg = iotxt)
        if(ioerr /= 0) then
          write (*, *) 'Problem openning probing file'
          write (*, *) 'Message: ', trim (iotxt)
          stop
        end if
!----------------------------------------------------------------------------------------------------------
! write out local data
!----------------------------------------------------------------------------------------------------------
        ix = dm%probexid(1, j)
        iy = dm%probexid(2, j)
        iz = dm%probexid(3, j)
        write(myunit, *) fl%time, fl%qx(ix, iy, iz), fl%qy(ix, iy, iz), fl%qz(ix, iy, iz)
        close(myunit)
      end if
    end do

    return
  end subroutine 

!==========================================================================================================
  subroutine write_monitor_thermo(tm, dm)
    use typeconvert_mod
    use parameters_constant_mod
    use wtformat_mod
    use udf_type_mod
    use files_io_mod
    use io_tools_mod
    implicit none 

    type(t_domain),  intent(in) :: dm
    type(t_thermo), intent(in) :: tm

    character(len=120) :: flname
    character(len=120) :: keyword
    character(200) :: iotxt
    integer :: ioerr, myunit
    integer :: ix, iy, iz
    integer :: i, j
!----------------------------------------------------------------------------------------------------------
! based on x-pencil
!----------------------------------------------------------------------------------------------------------
    j = 0
    do i = 1, dm%proben
      if( dm%probe_is_in(i) ) j = j + 1
      if(j > 0) then
!----------------------------------------------------------------------------------------------------------
! open file
!----------------------------------------------------------------------------------------------------------
        keyword = "monitor_pt"//trim(int2str(i))//"_thermo"
        call generate_pathfile_name(flname, dm%idom, keyword, dir_moni, 'dat')
        open(newunit = myunit, file = trim(flname), status = "old", action = "write", position = "append", &
            iostat = ioerr, iomsg = iotxt)
        if(ioerr /= 0) then
          write (*, *) 'Problem openning probing file'
          write (*, *) 'Message: ', trim (iotxt)
          stop
        end if
!----------------------------------------------------------------------------------------------------------
! write out local data
!----------------------------------------------------------------------------------------------------------
        ix = dm%probexid(1, j)
        iy = dm%probexid(2, j)
        iz = dm%probexid(3, j)
        write(myunit, *) tm%time, tm%tTemp(ix, iy, iz)
        close(myunit)
      end if
    end do

    return
  end subroutine 

!==========================================================================================================
end module