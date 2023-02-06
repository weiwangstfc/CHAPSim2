module boundary_conditions_mod


  public  :: configure_bc
  private :: map_bc_1d_uprofile
  public  :: buildup_bc_profile1d

  public  :: update_bc_interface_flow
  public  :: update_bc_interface_thermo
  public  :: apply_bc_const
  public  :: apply_convective_outlet
  


contains
!==========================================================================================================
!==========================================================================================================
  subroutine configure_bc(dm)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm


    integer :: m, n

!----------------------------------------------------------------------------------------------------------
! to set up real bc for calculation from given nominal b.c.
!----------------------------------------------------------------------------------------------------------
    do n = 1, 2
      do m = 1, 5
        
        if(dm%ibcx_nominal(n, m) == IBC_PROFILE1D)   then
            dm%ibcx(n, m) = IBC_DIRICHLET
        else if(dm%ibcx_nominal(n, m) == IBC_TURBGEN  .or. &
                dm%ibcx_nominal(n, m) == IBC_DATABASE )   then
            dm%ibcx(n, m) = IBC_INTERIOR
        else
            dm%ibcx(n, m) = dm%ibcx_nominal(n, m)   
        end if

        if(dm%ibcy_nominal(n, m) == IBC_PROFILE1D)   then
            dm%ibcy(n, m) = IBC_DIRICHLET
        else if(dm%ibcy_nominal(n, m) == IBC_TURBGEN  .or. &
                dm%ibcy_nominal(n, m) == IBC_DATABASE )   then
            dm%ibcy(n, m) = IBC_INTERIOR
        else
            dm%ibcy(n, m) = dm%ibcy_nominal(n, m)   
        end if

        if(dm%ibcz_nominal(n, m) == IBC_PROFILE1D)   then
            dm%ibcz(n, m) = IBC_DIRICHLET
        else if(dm%ibcz_nominal(n, m) == IBC_TURBGEN  .or. &
                dm%ibcz_nominal(n, m) == IBC_DATABASE )   then
            dm%ibcz(n, m) = IBC_INTERIOR
        else
            dm%ibcz(n, m) = dm%ibcz_nominal(n, m)   
        end if

      end do
    end do

!----------------------------------------------------------------------------------------------------------
! to exclude non-resonable input
!----------------------------------------------------------------------------------------------------------
    do m = 1, 5
      if(dm%ibcx_nominal(2, m) == IBC_PROFILE1D) call Print_error_msg(" This BC is not supported.")
      do n = 1, 2
        if(dm%ibcx(n, m)         >  IBC_OTHERS   ) dm%ibcx(n, m) = IBC_INTRPL
        if(dm%ibcy(n, m)         >  IBC_OTHERS   ) dm%ibcy(n, m) = IBC_INTRPL
        if(dm%ibcz(n, m)         >  IBC_OTHERS   ) dm%ibcz(n, m) = IBC_INTRPL
        if(dm%ibcy(n, m)         == IBC_INTERIOR ) call Print_error_msg(" This BC is not supported.")
        if(dm%ibcz(n, m)         == IBC_INTERIOR ) call Print_error_msg(" This BC is not supported.")
        if(dm%ibcy_nominal(n, m) == IBC_PROFILE1D) call Print_error_msg(" This BC is not supported.")
        if(dm%ibcz_nominal(n, m) == IBC_PROFILE1D) call Print_error_msg(" This BC is not supported.")
      end do
    end do

!----------------------------------------------------------------------------------------------------------
! to set up real bc values for calculation from given nominal b.c. values
! np, not nc, is used to allocate to provide enough space
!----------------------------------------------------------------------------------------------------------
    allocate( dm%fbcx_var(4,        dm%np(2), dm%np(3), 5) )
    allocate( dm%fbcy_var(dm%np(1), 4,        dm%np(3), 5) )
    allocate( dm%fbcz_var(dm%np(1), dm%np(2), 4,        5) )

    do m = 1, 5
      do n = 1, 2
        dm%fbcx_var(n, :, :, m) = dm%fbcx_const(n, m)
        dm%fbcy_var(:, n, :, m) = dm%fbcy_const(n, m)
        dm%fbcz_var(:, :, n, m) = dm%fbcz_const(n, m)
      end do
      do n = 3, 4
        dm%fbcx_var(n, :, :, m) = dm%fbcx_const(n - 2, m)
        dm%fbcy_var(:, n, :, m) = dm%fbcy_const(n - 2, m)
        dm%fbcz_var(:, :, n, m) = dm%fbcz_const(n - 2, m)
      end do

    end do

!----------------------------------------------------------------------------------------------------------
! to build up bc
!----------------------------------------------------------------------------------------------------------
    call buildup_bc_profile1d(dm)


    return
  end subroutine 

!==========================================================================================================
!==========================================================================================================
  subroutine  map_bc_1d_uprofile(filename, n, y, u)
    use parameters_constant_mod 
    use udf_type_mod
    implicit none 
    character(*), intent(in) :: filename
    integer, intent(in) :: n
    integer, intent(in) :: y(n)
    integer, intent(out) :: u(n)

    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit
    character(len = 80) :: str
    
    real(WP) :: rtmp
    integer :: i, nn
    real(WP), allocatable :: uprofile(:)
    real(WP), allocatable :: yy(:)

    !----------------------------------------------------------------------------------------------------------
    ! to read given u-velocity profile, dimensionless, x/H, U/Umean
    !----------------------------------------------------------------------------------------------------------

    open ( newunit = inputUnit,     &
           file    = trim(filename), &
           status  = 'old',         &
           action  = 'read',        &
           iostat  = ioerr,         &
           iomsg   = iotxt)
    if(ioerr /= 0) then
      write (*, *) 'Problem openning : '//trim(filename)
      write (*, *) 'Message: ', trim (iotxt)
      stop 4
    end if

    nn = 0
    read(inputUnit, *, iostat = ioerr) str

    do
      read(inputUnit, *, iostat = ioerr) rtmp, rtmp
      if(ioerr /= 0) exit
      nn = nn + 1
    end do
    rewind(inputUnit)
    !----------------------------------------------------------------------------------------------------------
    ! to read given u-velocity profile, dimensionless, x/H, U/Umean
    !----------------------------------------------------------------------------------------------------------
    allocate ( uprofile (nn) )
    allocate ( yy (nn) )

    read(inputUnit, *, iostat = ioerr) str
    do i = 1, nn
      read(inputUnit, *, iostat = ioerr) yy(i), uprofile(i)
    end do
    close(inputUnit)


    call map_1d_profile_to_case(nn, yy, uprofile, n, y, u)


    deallocate(uprofile)
    deallocate(yy)

    return
  end subroutine
  
!==========================================================================================================
!==========================================================================================================
  subroutine buildup_bc_profile1d(dm)
    implicit none
    type(t_domain), intent(inout) :: dm

    real(WP) :: var1y(dm%np(2))
    character(len=10) :: filename(5)

    filename(1) = 'bc_u1y.dat'
    filename(2) = 'bc_v1y.dat'
    filename(3) = 'bc_w1y.dat'
    filename(4) = 'bc_p1y.dat'
    filename(5) = 'bc_T1y.dat'

!   x-inlet profile, u(y), w(y), p(y), T(y) profile
    do m = 1, 5
      if(m == 2) cycle
      if(dm%ibcx_nominal(1, m) == IBC_PROFILE1D) then
        call map_bc_1d_uprofile(filename(i), dm%nc(2), dm%yc, var1y(1:dm%nc(2)))
        do k = 1, dm%nc(3)
          do j = 1, dm%nc(2)
            dm%fbcx_var(1, j, k, m) = var1y(j)
          end do
        end do
      end if
    end do
!   x-inlet profile, v(y)
    m = 2
    if(dm%ibcx_nominal(1, m) == IBC_PROFILE1D) then
      call map_bc_1d_uprofile(filename(i), dm%np(2), dm%yp, var1y(1:dm%np(2)))
      do k = 1, dm%nc(3)
        do j = 1, dm%np(2)
          dm%fbcx_var(1, j, k, m) = var1y(j)
        end do
      end do
    end if


    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine update_bc_interface_flow(dm0, fl0, dm1, fl1)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm0, dm1
    type(t_flow), intent(in)      :: fl0, fl1
    
    integer :: m
!----------------------------------------------------------------------------------------------------------
!   all in x-pencil
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
!   bc in x - direction
!----------------------------------------------------------------------------------------------------------
    ! ux, dm0-dm1, np
    m = 1
    if(dm1%ibcx(1, m) == IBC_INTERIOR) then
      dm1%fbcx_var(1, :, :, m) = fl0%qx(dm0%np(1),     :, :)
      dm1%fbcx_var(3, :, :, m) = fl0%qx(dm0%np(1) - 1, :, :)
    end if
    if(dm0%ibcx(2, m) == IBC_INTERIOR) then
      dm0%fbcx_var(2, :, :, m) = fl1%qx(1, :, :)
      dm0%fbcx_var(4, :, :, m) = fl1%qx(2, :, :)
    end if

    ! uy, dm0-dm1, nc
    m = 2
    if(dm1%ibcx(1, m) == IBC_INTERIOR) then
      dm1%fbcx_var(1, :, :, m) = fl0%qy(dm0%nc(1),     :, :)
      dm1%fbcx_var(3, :, :, m) = fl0%qy(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx(2, m) == IBC_INTERIOR) then
      dm0%fbcx_var(2, :, :, m) = fl1%qy(1, :, :)
      dm0%fbcx_var(4, :, :, m) = fl1%qy(2, :, :)
    end if

    ! uz, dm0-dm1, nc
    m = 3
    if(dm1%ibcx(1, m) == IBC_INTERIOR) then
      dm1%fbcx_var(1, :, :, m) = fl0%qz(dm0%nc(1),     :, :)
      dm1%fbcx_var(3, :, :, m) = fl0%qz(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx(2, m) == IBC_INTERIOR) then
      dm0%fbcx_var(2, :, :, m) = fl1%qz(1, :, :)
      dm0%fbcx_var(4, :, :, m) = fl1%qz(2, :, :)
    end if

    ! p, dm0-dm1, nc
    m = 4
    if(dm1%ibcx(1, m) == IBC_INTERIOR) then
      dm1%fbcx_var(1, :, :, m) = fl0%pres(dm0%nc(1),     :, :)
      dm1%fbcx_var(3, :, :, m) = fl0%pres(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx(2, m) == IBC_INTERIOR) then
      dm0%fbcx_var(2, :, :, m) = fl1%pres(1, :, :)
      dm0%fbcx_var(4, :, :, m) = fl1%pres(2, :, :)
    end if


    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine update_bc_interface_thermo(dm0, tm0, dm1, tm1)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm0, dm1
    type(t_thermo), intent(in)    :: tm0, tm1

    integer :: m
    
!----------------------------------------------------------------------------------------------------------
!   all in x-pencil, bc in x - direction
!----------------------------------------------------------------------------------------------------------
    ! T, dm0-dm1, nc
    m = 5
    if(dm1%ibcx(1, m) == IBC_INTERIOR) then
      dm1%fbcx_var(1, :, :, m) = tm0%tTemp(dm0%nc(1),     :, :)
      dm1%fbcx_var(3, :, :, m) = tm0%tTemp(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx(2, m) == IBC_INTERIOR) then
      dm0%fbcx_var(2, :, :, m) = tm1%tTemp(1, :, :)
      dm0%fbcx_var(4, :, :, m) = tm1%tTemp(2, :, :)
    end if

    return
  end subroutine
!==========================================================================================================
!> \brief Apply b.c. conditions 
!---------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   public
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!----------------------------------------------------------------------------------------------------------
!> \param[in]     d             domain
!> \param[out]    f             flow
!==========================================================================================================
  subroutine apply_bc_const (dm, fl) ! check, necessary?
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent( in    )   :: dm
    type(t_flow),   intent( inout )   :: fl

    integer :: m, s
    type(DECOMP_INFO) :: dtmp


     if(dm%ibcx(1, 1) /= IBC_DIRICHLET .and. &
        dm%ibcx(2, 1) /= IBC_DIRICHLET .and. &
        dm%ibcy(1, 2) /= IBC_DIRICHLET .and. &
        dm%ibcy(2, 2) /= IBC_DIRICHLET .and. &
        dm%ibcz(1, 3) /= IBC_DIRICHLET .and. &
        dm%ibcz(2, 3) /= IBC_DIRICHLET ) return
!----------------------------------------------------------------------------------------------------------
!   all in x-pencil
!----------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------
!   ux at x-direction. BC of others at x-direction are given in oeprations directly.
!----------------------------------------------------------------------------------------------------------
    m = 1
    dtmp = dm%dpcc
    ! constant value bc.
    do s = 1, 2
      if(dm%ibcx_nominal(s, m) == IBC_DIRICHLET) then
        if(dtmp%xst(m) == 1) then
          fl%qx(1, :, :) = dm%fbcx(s, m)
          if(dm%is_thermo) fl%gx(1, :, :) = fl%qx(1, :, :) * dm%fbc_dend(s, m)
        end if
        if(dtmp%xen(m) == dm%np(m)) then
          fl%qx(dtmp%xsz(m), :, :) = dm%fbcx(s, m)
          if(dm%is_thermo) fl%gx(dtmp%xsz(m), :, :) = fl%qx(dtmp%xsz(m), :, :) * dm%fbc_dend(s, m)
        end if
      end if
    end do

    ! variation value bc
    if(dm%ibcx_nominal(1, 1) == IBC_PROFILE1D .and. dtmp%xst(1) == 1) then
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        do k = 1, dtmp%xsz(3)
          kk = dtmp%xsz(3) + k - 1
          fl%qx(1, j, k) = dm%fbcxinlet(jj, kk, 1)
          if(dm%is_thermo) fl%gx(1, :, :) = fl%qx(1, :, :) * dm%fbc_dend(1, 1)
        end do
      end do
    end if

!----------------------------------------------------------------------------------------------------------
!   uy at y-direction. BC of others at y-direction are given in oeprations directly.
!----------------------------------------------------------------------------------------------------------
    m = 2
    dtmp = dm%dcpc
    do s = 1, 2
      if(dm%ibcy_nominal(s, m) == IBC_DIRICHLET) then
        if(dtmp%xst(m) == 1) then
          fl%qy(:, 1, :) = dm%fbcy(s, m)
          if(dm%is_thermo) fl%gy(:, 1, :) = fl%qy(:, 1, :) * dm%fbc_dend(s, m)
        end if
        if(dtmp%xen(m) == dm%np(m)) then
          fl%qy(:, dtmp%xsz(m), :) = dm%fbcy(s, m)
          if(dm%is_thermo) fl%gy(:, dtmp%xsz(m), :) = fl%qy(:, dtmp%xsz(m), :)  * dm%fbc_dend(s, m)
        end if
      end if
    end do

!----------------------------------------------------------------------------------------------------------
!   uz at z-direction. BC of others at z-direction are given in oeprations directly.
!----------------------------------------------------------------------------------------------------------
    m = 3
    dtmp = dm%dccp
    do s = 1, 2
      if(dm%ibcz_nominal(s, m) == IBC_DIRICHLET) then
        if(dtmp%xst(m) == 1) then
          fl%qz(:, :, 1) = dm%fbcz(s, m)
          if(dm%is_thermo) fl%gz(:, :, 1) = fl%qz(:, :, 1) * dm%fbc_dend(s, m)
        end if
        if(dtmp%xen(m) == dm%np(m)) then
          fl%qz(:, :, dtmp%xsz(m)) = dm%fbcz(s, m)
          if(dm%is_thermo)  fl%gz(:, :, dtmp%xsz(m)) = fl%qz(:, :, dtmp%xsz(m)) *  dm%fbc_dend(s, m)
        end if
      end if
    end do

    return
  end subroutine



!==========================================================================================================
!==========================================================================================================
  subroutine Apply_convective_outlet


  end subroutine



end module
