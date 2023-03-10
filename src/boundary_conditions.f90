module boundary_conditions_mod

  character(12) :: filename(5)
  
  private :: map_bc_1d_uprofile
  private :: apply_flow_bc_const
  public  :: configure_bc_type
  public  :: buildup_flow_bc_const
  public  :: update_flow_bc_interface
  public  :: update_thermo_bc_interface

  !public  :: apply_bc_const
  !public  :: apply_convective_outlet
  
contains

!==========================================================================================================
!==========================================================================================================
  subroutine  map_bc_1d_uprofile(filename, n, y, u)
    use parameters_constant_mod 
    use udf_type_mod
    implicit none 
    character(*), intent(in) :: filename
    integer,  intent(in)  :: n
    real(WP), intent(in)  :: y(n)
    real(WP), intent(out) :: u(n)

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
      write (*, *) ' File does not exsit : '//trim(filename)
      return
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
  subroutine apply_flow_bc_const (fl, dm) ! apply once only
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent( in)   :: dm
    type(t_flow), intent(inout)   :: fl

    real(WP) :: var1y(1:dm%np(2))
    
    integer :: m, n, k, ny
!==========================================================================================================
! to build up bc with constant values
! -3-1-||||-2-4
! for constant bc, 3=1= geometric bc, side 1;
!                  2=4= geometric bc, side 2 
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
! x-bc, qx
!----------------------------------------------------------------------------------------------------------
    do k = 1, dm%dpcc%xsz(3)
      do j = 1, dm%dpcc%xsz(2)
        do n = 1, 2
          fl%fbcx_qx(n,   j, k) =  dm%fbcx_const(n, 1)
          fl%fbcx_qx(n+2, j, k) =  fl%fbcx_qx(n, j, k)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! x-bc, qy
!----------------------------------------------------------------------------------------------------------
    do k = 1, dm%dcpc%xsz(3)
      do j = 1, dm%dcpc%xsz(2)
        jj = dm%dcpc%xst(2) + j - 1
        do n = 1, 2
          fl%fbcx_qy(n,   j, k) =  dm%fbcx_const(n, 2) / dm%rpi(jj)
          fl%fbcx_qy(n+2, j, k) =  fl%fbcx_qy(n, j, k)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! x-bc, qz, gz
!----------------------------------------------------------------------------------------------------------
    do k = 1, dm%dccp%xsz(3)
      do j = 1, dm%dccp%xsz(2)
        jj = dm%dccp%xst(2) + j - 1
        do n = 1, 2
          fl%fbcx_qz(n,   j, k) =  dm%fbcx_const(n, 3) / dm%rci(jj)
          fl%fbcx_qz(n+2, j, k) =  fl%fbcx_qz(n, j, k)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! x-bc, pr
!----------------------------------------------------------------------------------------------------------
    do k = 1, dm%dccc%xsz(3)
      do j = 1, dm%dccc%xsz(2)
        do n = 1, 2
          fl%fbcx_pr(n,   j, k) =  dm%fbcx_const(n, 4)
          fl%fbcx_pr(n+2, j, k) =  fl%fbcx_pr(n, j, k)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! y-bc, qx
!----------------------------------------------------------------------------------------------------------
    do k = 1, dm%dpcc%ysz(3)
      do i = 1, dm%dpcc%ysz(1)
        do n = 1, 2
          fl%fbcy_qx(i, n,   k) =  dm%fbcy_const(n, 1)
          fl%fbcy_qx(i, n+2, k) =  fl%fbcy_qx(i, n, k)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! y-bc, qy, gy
!----------------------------------------------------------------------------------------------------------
    do k = 1, dm%dcpc%ysz(3)
      do i = 1, dm%dcpc%ysz(1)
        do n = 1, 2
          if (n == 1) jj = 1
          if (n == 2) jj = dm%np_geo(2)
          fl%fbcy_qy(i, n,   k) =  dm%fbcy_const(n, 2) / dm%rpi(jj) ! check  rpi=maxp?
          fl%fbcy_qy(i, n+2, k) =  fl%fbcy_qy(i, n, k)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! y-bc, qz, gz
!----------------------------------------------------------------------------------------------------------
    do k = 1, dm%dccp%ysz(3)
      do i = 1, dm%dccp%ysz(1)
        do n = 1, 2
          if (n == 1) jj = 1
          if (n == 2) jj = dm%np_geo(2)
          fl%fbcy_qz(i, n,   k) =  dm%fbcy_const(n, 3) / dm%rpi(jj) ! check, rpi not rci
          fl%fbcy_qz(i, n+2, k) =  fl%fbcy_qz(i, n, k)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! y-bc, pr
!----------------------------------------------------------------------------------------------------------
    do k = 1, dm%dccc%ysz(3)
      do i = 1, dm%dccc%ysz(1)
        do n = 1, 2
          fl%fbcy_pr(i, n,   k) =  dm%fbcy_const(n, 4)
          fl%fbcy_pr(i, n+2, k) =  fl%fbcy_pr(i, n, k)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! z-bc, qx, gx
!----------------------------------------------------------------------------------------------------------
    do j = 1, dm%dpcc%zsz(2)
      do i = 1, dm%dpcc%zsz(1)
        do n = 1, 2
          fl%fbcz_qx(i, j, n  ) =  dm%fbcz_const(n, 1)
          fl%fbcz_qx(i, j, n+2) =  fl%fbcz_qx(i, j, n)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! z-bc, qy, gy
!----------------------------------------------------------------------------------------------------------
    do j = 1, dm%dcpc%zsz(2)
      jj = dm%dcpc%zst(2) + j - 1
      do i = 1, dm%dcpc%zsz(1)
        do n = 1, 2
          fl%fbcz_qy(i, j, n  ) =  dm%fbcz_const(n, 2) / dm%rpi(jj)
          fl%fbcz_qy(i, j, n+2) =  fl%fbcz_qx(i, j, n)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! z-bc, qz, gz
!----------------------------------------------------------------------------------------------------------
    do j = 1, dm%dccp%zsz(2)
      jj = dm%dccp%zst(2) + j - 1
      do i = 1, dm%dccp%zsz(1)
        do n = 1, 2
          fl%fbcz_qz(i, j, n  ) =  dm%fbcz_const(n, 3) / dm%rci(jj)
          fl%fbcz_qz(i, j, n+2) =  fl%fbcz_qz(i, j, n)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! z-bc, pr
!----------------------------------------------------------------------------------------------------------
    do j = 1, dm%dccc%zsz(2)
      do i = 1, dm%dccc%zsz(1)
        do n = 1, 2
          fl%fbcz_pr(i, j, n, ) =  dm%fbcz_const(n, 4)
          fl%fbcz_pr(i, j, n+2) =  fl%fbcz_pr(i, j, n)
        end do
      end do
    end do

!==========================================================================================================
! to build up bc for var(x_const, y, z)
!==========================================================================================================
    filename(1) = 'pf1d_u1y.dat' !(undim)
    filename(2) = 'pf1d_v1y.dat' !(undim)
    filename(3) = 'pf1d_w1y.dat' !(undim)
    filename(4) = 'pf1d_p1y.dat' !(undim)
    filename(5) = 'pf1d_T1y.dat' !(dim  )
!----------------------------------------------------------------------------------------------------------
! x-bc1, T(x_c, y, z), dim
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 5) == IBC_PROFILE1D) then

      var1y = ZERO
      ny = dm%nc(2)
      call map_bc_1d_uprofile( filename(5), ny, dm%yc, var1y(1:ny) )

      do k = 1, dm%dccc%xsz(3)
        do j = 1, dm%dccc%xsz(2)
          jj = dm%dccc%xst(2) + j - 1

          tm%fbcx_ftp(1, j, k)%t = var1y(jj) / fluidparam%ftp0ref%t
          tm%fbcx_ftp(3, j, k)%t = tm%fbcx_ftp(1, j, k)%t
          call ftp_refresh_thermal_properties_from_T_undim( tm%fbcx_ftp(1, j, k)%t )
          call ftp_refresh_thermal_properties_from_T_undim( tm%fbcx_ftp(3, j, k)%t )
        end do
      end do

    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, qx(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 1) == IBC_PROFILE1D) then

      var1y = ZERO
      ny = dm%nc(2)
      call map_bc_1d_uprofile( filename(1), ny, dm%yc, var1y(1:ny) )

      do k = 1, dm%dpcc%xsz(3)
        do j = 1, dm%dpcc%xsz(2)
          jj = dm%dpcc%xst(2) + j - 1
          fl%fbcx_qx(1, j, k) = var1y(jj)
          fl%fbcx_qx(3, j, k) = fl%fbcx_qx(1, j, k)
        end do
      end do

    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, qy(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 2) == IBC_PROFILE1D) then

      var1y = ZERO
      ny = dm%np(2)
      call map_bc_1d_uprofile( filename(2), ny, dm%yp, var1y(1:ny) )

      do k = 1, dm%dcpc%xsz(3)
        do j = 1, dm%dcpc%xsz(2)
          jj = dm%dpcc%xst(2) + j - 1
          fl%fbcx_qy(1, j, k) = var1y(jj) / dm%rpi(jj)
          fl%fbcx_qy(3, j, k) = fl%fbcx_qy(1, j, k)
        end do
      end do
      
    end if

!----------------------------------------------------------------------------------------------------------
! x-bc1, qz(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 3) == IBC_PROFILE1D) then

      var1y = ZERO
      ny = dm%nc(2)
      call map_bc_1d_uprofile( filename(3), ny, dm%yc, var1y(1:ny) )

      do k = 1, dm%dccp%xsz(3)
        do j = 1, dm%dccp%xsz(2)
          jj = dm%dccp%xst(2) + j - 1
          fl%fbcx_qz(1, j, k) = var1y(jj) / dm%rci(jj)
          fl%fbcx_qz(3, j, k) = fl%fbcx_qz(1, j, k)
        end do
      end do
      
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, pr(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 4) == IBC_PROFILE1D) then

      var1y = ZERO
      ny = dm%nc(2)
      call map_bc_1d_uprofile( filename(4), ny, dm%yc, var1y(1:ny) )

      do k = 1, dm%dccc%xsz(3)
        do j = 1, dm%dccc%xsz(2)
          jj = dm%dccc%xst(2) + j - 1
          fl%fbcx_qr(1, j, k) = var1y(jj)
          fl%fbcx_qr(3, j, k) = fl%fbcx_pr(1, j, k)
        end do
      end do

    end if

    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine configure_bc_type(dm)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm


    integer :: m, n

!----------------------------------------------------------------------------------------------------------
! to set up real bc for calculation from given nominal b.c.
!----------------------------------------------------------------------------------------------------------
    do n = 1, 2
      do m = 1, NBC
!----------------------------------------------------------------------------------------------------------
! x-direction
!----------------------------------------------------------------------------------------------------------
        if (dm%ibcx_nominal(n, m) == IBC_PROFILE1D)   then
          dm%ibcx(n, m) = IBC_DIRICHLET
        else if (dm%ibcx_nominal(n, m) == IBC_TURBGEN  .or. &
                 dm%ibcx_nominal(n, m) == IBC_DATABASE )   then
          if(m /=5) then
            dm%ibcx(n, m) = IBC_INTERIOR ! for temperature, default is no incoming thermal flow, check
          else 
            dm%ibcx(n, m) = IBC_DIRICHLET
          end if
        else if (dm%ibcx_nominal(n, m) == IBC_CONVECTIVE)   then
          dm%ibcx(n, m) = IBC_INTRPL
        else
          dm%ibcx(n, m) = dm%ibcx_nominal(n, m)   
        end if
!----------------------------------------------------------------------------------------------------------
! y-direction
!----------------------------------------------------------------------------------------------------------
        if (dm%ibcy_nominal(n, m) == IBC_PROFILE1D)   then
          dm%ibcy(n, m) = IBC_DIRICHLET
        else if (dm%ibcy_nominal(n, m) == IBC_TURBGEN  .or. &
                 dm%ibcy_nominal(n, m) == IBC_DATABASE )   then
          if(m /=5) then
            dm%ibcy(n, m) = IBC_INTERIOR
          else 
            dm%ibcy(n, m) = IBC_DIRICHLET
          end if
        else if (dm%ibcy_nominal(n, m) == IBC_CONVECTIVE)   then
          dm%ibcy(n, m) = IBC_INTRPL
        else
          dm%ibcy(n, m) = dm%ibcy_nominal(n, m)   
        end if
!----------------------------------------------------------------------------------------------------------
! z-direction
!----------------------------------------------------------------------------------------------------------
        if (dm%ibcz_nominal(n, m) == IBC_PROFILE1D)   then
          dm%ibcz(n, m) = IBC_DIRICHLET
        else if (dm%ibcz_nominal(n, m) == IBC_TURBGEN  .or. &
                 dm%ibcz_nominal(n, m) == IBC_DATABASE )   then
          if(m /=5) then
            dm%ibcz(n, m) = IBC_INTERIOR
          else 
            dm%ibcz(n, m) = IBC_DIRICHLET
          end if
        else if (dm%ibcz_nominal(n, m) == IBC_CONVECTIVE)   then
          dm%ibcz(n, m) = IBC_INTRPL
        else
          dm%ibcz(n, m) = dm%ibcz_nominal(n, m)   
        end if

      end do
    end do
!----------------------------------------------------------------------------------------------------------
! to exclude non-resonable input
!----------------------------------------------------------------------------------------------------------
    do m = 1, NBC
      if(dm%ibcx_nominal(2, m) == IBC_PROFILE1D) call Print_error_msg(" This BC IBC_PROFILE1D is not supported.")
      do n = 1, 2
        if(dm%ibcx(n, m)         >  IBC_OTHERS   ) dm%ibcx(n, m) = IBC_INTRPL
        if(dm%ibcy(n, m)         >  IBC_OTHERS   ) dm%ibcy(n, m) = IBC_INTRPL
        if(dm%ibcz(n, m)         >  IBC_OTHERS   ) dm%ibcz(n, m) = IBC_INTRPL
        if(dm%ibcy_nominal(n, m) == IBC_PROFILE1D) call Print_error_msg(" This BC IBC_PROFILE1D is not supported.")
        if(dm%ibcz_nominal(n, m) == IBC_PROFILE1D) call Print_error_msg(" This BC IBC_PROFILE1D is not supported.")
      end do
    end do

    return
  end subroutine 
!==========================================================================================================
!==========================================================================================================
  subroutine buildup_flow_bc_const(fl, dm)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in)  :: dm
    type(t_flow), intent(inout) :: fl

!----------------------------------------------------------------------------------------------------------
! to set up real bc values for calculation from given nominal b.c. values
!----------------------------------------------------------------------------------------------------------
      allocate( fl%fbcx_qx(4,              dm%dpcc%xsz(2), dm%dpcc%xsz(3)) )
      allocate( fl%fbcy_qx(dm%dpcc%ysz(1), 4,              dm%dpcc%ysz(3)) )
      allocate( fl%fbcz_qx(dm%dpcc%zsz(1), dm%dpcc%zsz(2), 4             ) )

      allocate( fl%fbcx_qy(4,              dm%dcpc%xsz(2), dm%dcpc%xsz(3)) )
      allocate( fl%fbcy_qy(dm%dcpc%ysz(1), 4,              dm%dcpc%ysz(3)) )
      allocate( fl%fbcz_qy(dm%dcpc%zsz(1), dm%dcpc%zsz(2), 4             ) )

      allocate( fl%fbcx_qz(4,              dm%dccp%xsz(2), dm%dccp%xsz(3)) )
      allocate( fl%fbcy_qz(dm%dccp%ysz(1), 4,              dm%dccp%ysz(3)) )
      allocate( fl%fbcz_qz(dm%dccp%zsz(1), dm%dccp%zsz(2), 4             ) )

      allocate( fl%fbcx_pr(4,              dm%dccc%xsz(2), dm%dccc%xsz(3)) )
      allocate( fl%fbcy_pr(dm%dccc%ysz(1), 4,              dm%dccc%ysz(3)) )
      allocate( fl%fbcz_pr(dm%dccc%zsz(1), dm%dccc%zsz(2), 4             ) )
    
    if(dm%is_thermo) then

      allocate( fl%fbcx_gx(4,              dm%dpcc%xsz(2), dm%dpcc%xsz(3)) )
      allocate( fl%fbcy_gx(dm%dpcc%ysz(1), 4,              dm%dpcc%ysz(3)) )
      allocate( fl%fbcz_gx(dm%dpcc%zsz(1), dm%dpcc%zsz(2), 4             ) )

      allocate( fl%fbcx_gy(4,              dm%dcpc%xsz(2), dm%dcpc%xsz(3)) )
      allocate( fl%fbcy_gy(dm%dcpc%ysz(1), 4,              dm%dcpc%ysz(3)) )
      allocate( fl%fbcz_gy(dm%dcpc%zsz(1), dm%dcpc%zsz(2), 4             ) )

      allocate( fl%fbcx_gz(4,              dm%dccp%xsz(2), dm%dccp%xsz(3)) )
      allocate( fl%fbcy_gz(dm%dccp%ysz(1), 4,              dm%dccp%ysz(3)) )
      allocate( fl%fbcz_gz(dm%dccp%zsz(1), dm%dccp%zsz(2), 4             ) )
      
    end if

    call apply_flow_bc_const(fl, dm)

    return
  end subroutine 

!==========================================================================================================
!==========================================================================================================
  subroutine update_flow_bc_interface(dm0, fl0, dm1, fl1)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm0, dm1
    type(t_flow), intent(in)      :: fl0, fl1
    
    integer :: m
!----------------------------------------------------------------------------------------------------------
!   all in x-pencil
!   no overlap of values
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
!   bc in x - direction
!----------------------------------------------------------------------------------------------------------
    ! ux, dm0-dm1, np
    if(dm1%ibcx(1, 1) == IBC_INTERIOR) then
      dm1%fbcx_qx(1, :, :) = fl0%qx(dm0%np_geo(1) - 1, :, :)
      dm1%fbcx_qx(3, :, :) = fl0%qx(dm0%np_geo(1) - 2, :, :)
    end if
    if(dm0%ibcx(2, 1) == IBC_INTERIOR) then
      dm0%fbcx_qx(2, :, :) = fl1%qx(2, :, :)
      dm0%fbcx_qx(4, :, :) = fl1%qx(3, :, :)
    end if

    ! uy, dm0-dm1, nc
    if(dm1%ibcx(1, 2) == IBC_INTERIOR) then
      dm1%fbcx_qy(1, :, :) = fl0%qy(dm0%nc(1),     :, :)
      dm1%fbcx_qy(3, :, :) = fl0%qy(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx(2, 2) == IBC_INTERIOR) then
      dm0%fbcx_qy(2, :, :) = fl1%qy(1, :, :)
      dm0%fbcx_qy(4, :, :) = fl1%qy(2, :, :)
    end if

    ! uz, dm0-dm1, nc
    if(dm1%ibcx(1, 3) == IBC_INTERIOR) then
      dm1%fbcx_qz(1, :, :) = fl0%qz(dm0%nc(1),     :, :)
      dm1%fbcx_qz(3, :, :) = fl0%qz(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx(2, 3) == IBC_INTERIOR) then
      dm0%fbcx_qz(2, :, :, m) = fl1%qz(1, :, :)
      dm0%fbcx_qz(4, :, :, m) = fl1%qz(2, :, :)
    end if

    ! p, dm0-dm1, nc
    if(dm1%ibcx(1, 4) == IBC_INTERIOR) then
      dm1%fbcx_pr(1, :, :) = fl0%pres(dm0%nc(1),     :, :)
      dm1%fbcx_pr(3, :, :) = fl0%pres(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx(2, m) == IBC_INTERIOR) then
      dm0%fbcx_pr(2, :, :) = fl1%pres(1, :, :)
      dm0%fbcx_pr(4, :, :) = fl1%pres(2, :, :)
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in y - direction
!----------------------------------------------------------------------------------------------------------
    ! ux, dm0-dm1, nc
    if(dm1%ibcy(1, 1) == IBC_INTERIOR) then
      dm1%fbcy_qx(:, 1, :) = fl0%qx(:, dm0%nc(2),     :)
      dm1%fbcy_qx(:, 3, :) = fl0%qx(:, dm0%nc(2) - 1, :)
    end if
    if(dm0%ibcy(2, 1) == IBC_INTERIOR) then
      dm0%fbcy_qx(:, 2, :) = fl1%qx(:, 1, :)
      dm0%fbcy_qx(:, 4, :) = fl1%qx(:, 2, :)
    end if

    ! uy, dm0-dm1, np
    if(dm1%ibcy(1, 2) == IBC_INTERIOR) then
      dm1%fbcy_qy(:, 1, :) = fl0%qy(:, dm0%np_geo(2) - 1, :)
      dm1%fbcy_qy(:, 3, :) = fl0%qy(:, dm0%np_geo(2) - 2, :)
    end if
    if(dm0%ibcy(2, 2) == IBC_INTERIOR) then
      dm0%fbcy_qy(:, 2, :) = fl1%qy(:, 2, :)
      dm0%fbcy_qy(:, 4, :) = fl1%qy(:, 3, :)
    end if

    ! uz, dm0-dm1, nc
    if(dm1%ibcy(1, 3) == IBC_INTERIOR) then
      dm1%fbcy_qz(:, 1, :) = fl0%qz(:, dm0%nc(2),     :)
      dm1%fbcy_qz(:, 3, :) = fl0%qz(:, dm0%nc(2) - 1, :)
    end if
    if(dm0%ibcy(2, 3) == IBC_INTERIOR) then
      dm0%fbcy_qz(:, 2, :) = fl1%qz(:, 1, :)
      dm0%fbcy_qz(:, 4, :) = fl1%qz(:, 2, :)
    end if

    ! p, dm0-dm1, nc
    if(dm1%ibcy(1, 4) == IBC_INTERIOR) then
      dm1%fbcy_pr(:, 1, :) = fl0%pres(:, dm0%nc(2),     :)
      dm1%fbcy_pr(:, 3, :) = fl0%pres(:, dm0%nc(2) - 1, :)
    end if
    if(dm0%ibcy(2, 4) == IBC_INTERIOR) then
      dm0%fbcy_pr(:, 2, :) = fl1%pres(:, 1, :)
      dm0%fbcy_pr(:, 4, :) = fl1%pres(:, 2, :)
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in z - direction
!----------------------------------------------------------------------------------------------------------
    ! ux, dm0-dm1, nc
    if(dm1%ibcz(1, 1) == IBC_INTERIOR) then
      dm1%fbcz_qx(:, :, 1) = fl0%qx(:, :, dm0%nc(2)    )
      dm1%fbcz_qx(:, :, 3) = fl0%qx(:, :, dm0%nc(2) - 1)
    end if
    if(dm0%ibcz(2, 1) == IBC_INTERIOR) then
      dm0%fbcz_qx(:, :, 2) = fl1%qx(:, :, 1)
      dm0%fbcz_qx(:, :, 4) = fl1%qx(:, :, 2)
    end if

    ! uy, dm0-dm1, nc
    if(dm1%ibcz(1, 2) == IBC_INTERIOR) then
      dm1%fbcz_qy(:, :, 1) = fl0%qy(:, :, dm0%nc(2)    )
      dm1%fbcz_qy(:, :, 3) = fl0%qy(:, :, dm0%nc(2) - 1)
    end if
    if(dm0%ibcz(2, 2) == IBC_INTERIOR) then
      dm0%fbcz_qy(:, :, 2) = fl1%qy(:, :, 1)
      dm0%fbcz_qy(:, :, 4) = fl1%qy(:, :, 2)
    end if

    ! uz, dm0-dm1, np
    if(dm1%ibcz(1, 3) == IBC_INTERIOR) then
      dm1%fbcz_qz(:, :, 1) = fl0%qz(:, :, dm0%np_geo(2) - 1)
      dm1%fbcz_qz(:, :, 3) = fl0%qz(:, :, dm0%np_geo(2) - 2)
    end if
    if(dm0%ibcz(2, 3) == IBC_INTERIOR) then
      dm0%fbcz_qz(:, :, 2) = fl1%qz(:, :, 2)
      dm0%fbcz_qz(:, :, 4) = fl1%qz(:, :, 3)
    end if

    ! p, dm0-dm1, nc
    if(dm1%ibcz(1, 4) == IBC_INTERIOR) then
      dm1%fbcz_pr(:, :, 1) = fl0%pres(:, :, dm0%nc(2)    )
      dm1%fbcz_pr(:, :, 3) = fl0%pres(:, :, dm0%nc(2) - 1)
    end if
    if(dm0%ibcz(2, 4) == IBC_INTERIOR) then
      dm0%fbcz_pr(:, :, 2) = fl1%pres(:, :, 1)
      dm0%fbcz_pr(:, :, 4) = fl1%pres(:, :, 2)
    end if
    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine update_thermo_bc_interface(dm0, fl0, tm0, dm1, fl1, tm1)
    use parameters_constant_mod
    use thermo_info_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm0, dm1
    type(t_flow), intent(in)      :: fl0, fl1
    type(t_thermo), intent(in)    :: tm0, tm1
    
    integer :: m, i, j, k
!----------------------------------------------------------------------------------------------------------
!   bc in x - direction
!----------------------------------------------------------------------------------------------------------
    ! thermal field, dm0-dm1
    if(dm1%ibcx(1, 5) == IBC_INTERIOR) then
      tm1%fbcx_ftp(1, :, :)%t = tm0%tTemp(dm0%nc(1),     :, :)
      tm1%fbcx_ftp(3, :, :)%t = tm0%tTemp(dm0%nc(1) - 1, :, :)
      do k = 1, dm1%dccc%xsz(3)
        do j = 1, dm1%dccc%xsz(2)
          call ftp_refresh_thermal_properties_from_T_undim(tm1%fbcx_ftp(1, j, k)%t)
          call ftp_refresh_thermal_properties_from_T_undim(tm1%fbcx_ftp(3, j, k)%t)
        end do
      end do
    end if
    if(dm0%ibcx(2, 5) == IBC_INTERIOR) then
      tm0%fbcx_ftp(2, :, :)%t = tm1%tTemp(1, :, :)
      tm0%fbcx_ftp(4, :, :)%t = tm1%tTemp(2, :, :)
      do k = 1, dm0%dccc%xsz(3)
        do j = 1, dm0%dccc%xsz(2)
          call ftp_refresh_thermal_properties_from_T_undim(dm0%fbcx_ftp(2, j, k)%t)
          call ftp_refresh_thermal_properties_from_T_undim(dm0%fbcx_ftp(4, j, k)%t)
        end do
      end do
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in y - direction
!----------------------------------------------------------------------------------------------------------
    ! thermal field, dm0-dm1
    if(dm1%ibcy(1, 5) == IBC_INTERIOR) then
      tm1%fbcy_ftp(:, 1, :)%t = tm0%tTemp(:, dm0%nc(1),     :)
      tm1%fbcy_ftp(:, 3, :)%t = tm0%tTemp(:, dm0%nc(1) - 1, :)
      do k = 1, dm1%dccc%ysz(3)
        do i = 1, dm1%dccc%ysz(1)
          call ftp_refresh_thermal_properties_from_T_undim(dm1%fbcy_ftp(i, 1, k)%t)
          call ftp_refresh_thermal_properties_from_T_undim(dm1%fbcy_ftp(i, 3, k)%t)
        end do
      end do
    end if
    if(dm0%ibcy(2, 5) == IBC_INTERIOR) then
      tm0%fbcy_ftp(:, 2, :)%t = tm1%tTemp(:, 1, :)
      tm0%fbcy_ftp(:, 4, :)%t = tm1%tTemp(:, 2, :)
      do k = 1, dm0%dccc%ysz(3)
        do i = 1, dm0%dccc%ysz(1)
          call ftp_refresh_thermal_properties_from_T_undim(tm0%fbcy_ftp(i, 2, k)%t)
          call ftp_refresh_thermal_properties_from_T_undim(tm0%fbcy_ftp(i, 4, k)%t)
        end do
      end do
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in z - direction
!----------------------------------------------------------------------------------------------------------
    ! thermal field, dm0-dm1
    if(dm1%ibcz(1, 5) == IBC_INTERIOR) then
      tm1%fbcz_ftp(:, :, 1)%t = tm0%tTemp(:, :, dm0%nc(1)    )
      tm1%fbcz_ftp(:, :, 3)%t = tm0%tTemp(:, :, dm0%nc(1) - 1)
      do j = 1, dm1%dccc%zsz(2)
        do i = 1, dm1%dccc%zsz(1)
          call ftp_refresh_thermal_properties_from_T_undim(tm1%fbcz_ftp(i, j, 1)%t)
          call ftp_refresh_thermal_properties_from_T_undim(tm1%fbcz_ftp(i, j, 3)%t)
        end do
      end do
    end if
    if(dm0%ibcz(2, 5) == IBC_INTERIOR) then
      tm0%fbcz_ftp(:, :, 2)%t = tm1%tTemp(:, :, 1)
      tm0%fbcz_ftp(:, :, 4)%t = tm1%tTemp(:, :, 2)
      do j = 1, dm0%dccc%zsz(2)
        do i = 1, dm0%dccc%zsz(1)
          call ftp_refresh_thermal_properties_from_T_undim(tm0%fbcz_ftp(i, j, 2)%t)
          call ftp_refresh_thermal_properties_from_T_undim(tm0%fbcz_ftp(i, j, 4)%t)
        end do
      end do
    end if
    
    return
  end subroutine

end module
