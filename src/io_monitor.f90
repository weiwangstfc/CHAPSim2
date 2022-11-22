module io_monitor_mod
  implicit none

  public :: write_monitor_ini
  public :: write_monitor
  
contains
  subroutine write_monitor_ini(dm)
    use typeconvert_mod
    use wtformat_mod
    use udf_type_mod
    use files_io_mod
    implicit none 
    type(t_domain),  intent(inout) :: dm

    integer :: myunit
    integer :: i, j
    logical :: exist
    character(len=100) :: flname

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
! create probe history file
!----------------------------------------------------------------------------------------------------------
    do i = 1, dm%proben
      if(.not. dm%probe_is_in(i)) cycle 
      flname = trim(dir_moni)//'/domain'//trim(int2str(dm%idom))//'_probe_pt'//trim(int2str(i))//'.log'
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

    return
  end subroutine

!==========================================================================================================
  subroutine write_monitor(dm, fl, tm)
    use typeconvert_mod
    use parameters_constant_mod
    use wtformat_mod
    use udf_type_mod
    use files_io_mod
    implicit none 

    type(t_domain),  intent(in) :: dm
    type(t_flow), intent(in) :: fl
    type(t_thermo), optional, intent(in) :: tm

    character(len=100) :: flname
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
        flname = trim(dir_moni)//'/domain'//trim(int2str(dm%idom))//'_probe_pt'//trim(int2str(i))//'.log'
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
end module