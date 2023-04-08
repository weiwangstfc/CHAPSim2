module boundary_conditions_mod

  character(12) :: filename(5)
  
  public  :: configure_bc_type

  private :: map_bc_1d_uprofile
  private :: apply_flow_bc_geo   ! constant given bc, apply once only  
  private :: apply_gxgygz_bc_geo ! constant given bc, apply once only
  public  :: buildup_thermo_bc_geo ! constant given bc, apply once only
  public  :: buildup_flow_bc_geo   ! constant given bc, apply once only
  
  public  :: update_thermo_bc_1dm_halo ! for pipe only, update every NS
  public  :: update_flow_bc_1dm_halo   ! for pipe only, update every NS 
  

  public  :: update_flow_bc_2dm_halo   ! for multiple domains, update every NS 
  public  :: update_thermo_bc_2dm_halo ! for multiple domains, update every NS
  
contains
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
! treatment for cylindrical coordinates
!----------------------------------------------------------------------------------------------------------
    if(dm(i)%icoordinate == ICYLINDRICAL) then
      dm%ibcz(:, :) = IBC_PERIODIC
    end if

    if (dm(i)%icase == ICASE_PIPE) then
      dm%ibcy(1, :) = IBC_INTERIOR
      dm%ibcy(1, :) = IBC_INTERIOR
      ! uy=0 at y=0
      dm%ibcy(1, 2) = IBC_DIRICHLET
      dm%fbcy_const(1, 2) = ZERO
    end if
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
!==========================================================================================================
  subroutine buildup_thermo_bc_geo(dm, tm)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(inout) :: tm

    integer :: i, j, k, n

    if(.not. dm%is_thermo) return

    !----------------------------------------------------------------------------------------------------------
    !   for x-pencil
    !   scale the given thermo b.c. in dimensional to undimensional
    !----------------------------------------------------------------------------------------------------------
    allocate( tm%fbcx_ftp(4,              dm%dccc%xsz(2), dm%dccc%xsz(3)) )
    allocate( tm%fbcy_ftp(dm%dccc%ysz(1), 4,              dm%dccc%ysz(3)) )
    allocate( tm%fbcz_ftp(dm%dccc%zsz(1), dm%dccc%zsz(2), 4             ) )

    allocate( tm%fbcx_heatflux(4,              dm%dccc%xsz(2), dm%dccc%xsz(3)) )
    allocate( tm%fbcy_heatflux(dm%dccc%ysz(1), 4,              dm%dccc%ysz(3)) )
    allocate( tm%fbcz_heatflux(dm%dccc%zsz(1), dm%dccc%zsz(2), 4             ) )

    !----------------------------------------------------------------------------------------------------------
    ! x-bc
    !----------------------------------------------------------------------------------------------------------
    do k = 1, dm%dccc%xsz(3)
      do j = 1, dm%dccc%xsz(2)
        do n = 1, 2
          if( dm%ibcx(n, 5) == IBC_DIRICHLET ) then
            ! dimensional T --> undimensional T
            tm%fbcx_ftp(n,   j, k)%t = dm%fbcx_const(n, 5) / fluidparam%ftp0ref%t
            tm%fbcx_ftp(n+2, j, k)%t = tm%fbcx_ftp(n, j, k)%t
            call ftp_refresh_thermal_properties_from_T_undim( tm%fbcx_ftp(n,   j, k)%t )
            call ftp_refresh_thermal_properties_from_T_undim( tm%fbcx_ftp(n+2, j, k)%t )

          else if (dm%ibcx(n, 5) == IBC_NEUMANN) then
            ! dimensional heat flux (k*dT/dx) --> undimensional heat flux (k*dT/dx)
            tm%fbcx_heatflux(n,   j, k)%t = dm%fbcx_const(n, 5) * tm%ref_l0 / fluidparam%ftp0ref%k / fluidparam%ftp0ref%t 
            tm%fbcx_heatflux(n+2, j, k)%t = tm%fbcx_heatflux(n, j, k)%t
          else
          end if

        end do 
      end do
    end do 

    !----------------------------------------------------------------------------------------------------------
    ! y-bc
    !----------------------------------------------------------------------------------------------------------
    do k = 1, dm%dccc%ysz(3)
      do i = 1, dm%dccc%ysz(1)
        do n = 1, 2
          if( dm%ibcy(n, 5) == IBC_DIRICHLET ) then
            ! dimensional T --> undimensional T
            tm%fbcy_ftp(i, n,   k)%t = dm%fbcy_const(n, 5) / fluidparam%ftp0ref%t
            tm%fbcy_ftp(i, n+2, k)%t = tm%fbcy_ftp(i, n, k)%t
            call ftp_refresh_thermal_properties_from_T_undim( tm%fbcy_ftp(i, n,   k)%t )
            call ftp_refresh_thermal_properties_from_T_undim( tm%fbcy_ftp(i, n+2, k)%t )

          else if (dm%ibcy(n, 5) == IBC_NEUMANN) then
            ! dimensional heat flux (k*dT/dx) --> undimensional heat flux (k*dT/dx)
            tm%fbcy_heatflux(i, n,   k)%t = dm%fbcy_const(n, 5) * tm%ref_l0 / fluidparam%ftp0ref%k / fluidparam%ftp0ref%t 
            tm%fbcy_heatflux(i, n+2, k)%t = tm%fbcy_heatflux(i, n, k)%t
          else
          end if

        end do 
      end do
    end do 

    !----------------------------------------------------------------------------------------------------------
    ! z-bc
    !----------------------------------------------------------------------------------------------------------
    do j = 1, dm%dccc%zsz(2)
      do i = 1, dm%dccc%zsz(1)
        do n = 1, 2
          if( dm%ibcz(n, 5) == IBC_DIRICHLET ) then
            ! dimensional T --> undimensional T
            tm%fbcz_ftp(i, j, n  )%t = dm%fbcz_const(n, 5) / fluidparam%ftp0ref%t
            tm%fbcz_ftp(i, j, n+2)%t = tm%fbcz_ftp(i, j, n)%t
            call ftp_refresh_thermal_properties_from_T_undim( tm%fbcz_ftp(i, j, n  )%t )
            call ftp_refresh_thermal_properties_from_T_undim( tm%fbcz_ftp(i, j, n+2)%t )

          else if (dm%ibcz(n, 5) == IBC_NEUMANN) then
            ! dimensional heat flux (k*dT/dx) --> undimensional heat flux (k*dT/dx)
            tm%fbcz_heatflux(i, j, n  )%t = dm%fbcz_const(n, 5) * tm%ref_l0 / fluidparam%ftp0ref%k / fluidparam%ftp0ref%t 
            tm%fbcz_heatflux(i, j, n+2)%t = tm%fbcz_heatflux(i, j, n)%t
          else
          end if

        end do 
      end do
    end do 

!----------------------------------------------------------------------------------------------------------
! x-bc1, T(x_c, y, z), dim
!----------------------------------------------------------------------------------------------------------
    filename(5) = 'pf1d_T1y.dat' !(dim  )
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
  subroutine apply_flow_bc_geo (dm, fl) ! apply once only
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
! x-bc, qz
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
! y-bc, qy
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
! y-bc, qz
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
! z-bc, qx
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
! z-bc, qy
!----------------------------------------------------------------------------------------------------------
    do j = 1, dm%dcpc%zsz(2)
      jj = dm%dcpc%zst(2) + j - 1
      do i = 1, dm%dcpc%zsz(1)
        do n = 1, 2
          fl%fbcz_qy(i, j, n  ) =  dm%fbcz_const(n, 2) / dm%rpi(jj)
          fl%fbcz_qy(i, j, n+2) =  fl%fbcz_qy(i, j, n)
        end do
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! z-bc, qz
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
! z-bc
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
! y bc for primary variables, uy, uz
!==========================================================================================================
    if(dm%icoordinate == ICYLINDRICAL) then
!----------------------------------------------------------------------------------------------------------
! y-bc, qyr = qy/r = uy
!----------------------------------------------------------------------------------------------------------
      do k = 1, dm%dcpc%ysz(3)
        do i = 1, dm%dcpc%ysz(1)
          do n = 1, 2
            if (n == 1) jj = 1
            if (n == 2) jj = dm%np_geo(2)
            fl%fbcy_qyr(i, n,   k) =  dm%fbcy_const(n, 2)
            fl%fbcy_qyr(i, n+2, k) =  fl%fbcy_qyr(i, n, k)
          end do
        end do
      end do
!----------------------------------------------------------------------------------------------------------
! y-bc, qzr = qz/r = uz
!----------------------------------------------------------------------------------------------------------
      do k = 1, dm%dccp%ysz(3)
        do i = 1, dm%dccp%ysz(1)
          do n = 1, 2
            if (n == 1) jj = 1
            if (n == 2) jj = dm%np_geo(2)
            fl%fbcy_qzr(i, n,   k) =  dm%fbcy_const(n, 3)
            fl%fbcy_qzr(i, n+2, k) =  fl%fbcy_qzr(i, n, k)
          end do
        end do
      end do
!----------------------------------------------------------------------------------------------------------
! z-bc, qyr = qy/r = uy
!----------------------------------------------------------------------------------------------------------
      do j = 1, dm%dcpc%zsz(2)
        jj = dm%dcpc%zst(2) + j - 1
        do i = 1, dm%dcpc%zsz(1)
          do n = 1, 2
            fl%fbcz_qyr(i, j, n  ) =  dm%fbcz_const(n, 2)
            fl%fbcz_qyr(i, j, n+2) =  fl%fbcz_qyr(i, j, n)
          end do
        end do
      end do
!----------------------------------------------------------------------------------------------------------
! z-bc, qzr = qz/r = uz
!----------------------------------------------------------------------------------------------------------
      do j = 1, dm%dccp%zsz(2)
        jj = dm%dccp%zst(2) + j - 1
        do i = 1, dm%dccp%zsz(1)
          do n = 1, 2
            fl%fbcz_qzr(i, j, n  ) =  dm%fbcz_const(n, 3)
            fl%fbcz_qzr(i, j, n+2) =  fl%fbcz_qzr(i, j, n)
          end do
        end do
      end do

    end if

!==========================================================================================================
! to build up bc for var(x_const, y, z)
!==========================================================================================================
    filename(1) = 'pf1d_u1y.dat' !(undim)
    filename(2) = 'pf1d_v1y.dat' !(undim)
    filename(3) = 'pf1d_w1y.dat' !(undim)
    filename(4) = 'pf1d_p1y.dat' !(undim)
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
          fl%fbcx_pr(1, j, k) = var1y(jj)
          fl%fbcx_pr(3, j, k) = fl%fbcx_pr(1, j, k)
        end do
      end do

    end if

    return
  end subroutine


!==========================================================================================================
!==========================================================================================================
  subroutine buildup_flow_bc_geo(dm, fl)
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

    if(dm%icoordinate == ICYLINDRICAL) then 
      allocate( fl%fbcy_qyr(dm%dcpc%ysz(1), 4,              dm%dcpc%ysz(3)) )
      allocate( fl%fbcz_qyr(dm%dcpc%zsz(1), dm%dcpc%zsz(2), 4             ) )
      allocate( fl%fbcy_qzr(dm%dccp%ysz(1), 4,              dm%dccp%ysz(3)) )
      allocate( fl%fbcz_qzr(dm%dccp%zsz(1), dm%dccp%zsz(2), 4             ) )
    end if

    call apply_flow_bc_geo(dm, fl)
    
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

      if(dm%icoordinate == ICYLINDRICAL) then 
        allocate( fl%fbcy_gyr(dm%dcpc%ysz(1), 4,              dm%dcpc%ysz(3)) )
        allocate( fl%fbcz_gyr(dm%dcpc%zsz(1), dm%dcpc%zsz(2), 4             ) )
        allocate( fl%fbcy_gzr(dm%dccp%ysz(1), 4,              dm%dccp%ysz(3)) )
        allocate( fl%fbcz_gzr(dm%dccp%zsz(1), dm%dccp%zsz(2), 4             ) )
      end if

      call apply_gxgygz_bc_geo(dm, fl)
      
    end if

    return
  end subroutine 



!==========================================================================================================
!==========================================================================================================
  subroutine update_flow_bc_1dm_halo(dm, fl) ! for cylindrical only
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    type(t_flow), intent(in)      :: fl

    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: apcc_ypencil
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: apcc_zpencil
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: apcc1_zpencil

    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp1_zpencil
    
    integer :: k

    if(dm%icase /= ICASE_PIPE) return
!----------------------------------------------------------------------------------------------------------
!   all in z-pencil
!   no overlap of values
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
!   (qx, qy, qz, pr) bc in y - direction
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy(1, 1) == IBC_INTERIOR) then
      call transpose_x_to_y(fl%qx,        apcc_ypencil, dm%apcc)
      call transpose_y_to_z(apcc_ypencil, apcc_zpencil, dm%apcc)

      do k = 1, dm%dpcc%zsz(3)
        apcc1_zpencil(:, :, k) = apcc_zpencil(:, :, dm%knc_sym(k))
      end do
      call transpose_z_to_y(apcc1_zpencil, apcc_ypencil, dm%apcc)

      fl%fbcy_qx(:, 1, :) = apcc_ypencil(:, 1, :)
      fl%fbcy_qx(:, 3, :) = apcc_ypencil(:, 2, :)
    end if

    ! uy, dm0-dm1, np
    fl%fbcy_qy(:, 1, :) = ZERO
    fl%fbcy_ur(:, 1, :) = ZERO

    ! uz, dm0-dm1, nc
    if(dm%ibcy(1, 3) == IBC_INTERIOR) then
      call transpose_x_to_y(fl%qz,        accp_ypencil, dm%accp)
      call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%accp)

      do k = 1, dm%dpcc%zsz(3)
        accp1_zpencil(:, :, k) = accp_zpencil(:, :, dm%knc_sym(k))
      end do
      call transpose_z_to_y(accp1_zpencil, accp_ypencil, dm%accp)

      fl%fbcy_qz(:, 1, :) = accp_ypencil(:, 1, :)
      fl%fbcy_qz(:, 3, :) = accp_ypencil(:, 2, :)

      fl%fbcy_uz(:, 1, :) = fl%fbcy_qz(:, 1, :) * dm%rci(1)
      fl%fbcy_uz(:, 3, :) = fl%fbcy_qz(:, 3, :) * dm%rci(2)

    end if

    ! p, dm0-dm1, nc
    if(dm%ibcy(1, 4) == IBC_INTERIOR) then
      call transpose_x_to_y(fl%pres,      accc_ypencil, dm%accc)
      call transpose_y_to_z(accc_ypencil, accc_zpencil, dm%accc)

      do k = 1, dm%dpcc%zsz(3)
        accc1_zpencil(:, :, k) = accc_zpencil(:, :, dm%knc_sym(k))
      end do
      call transpose_z_to_y(accc1_zpencil, accc_ypencil, dm%accc)

      fl%fbcy_pr(:, 1, :) = accp_ypencil(:, 1, :)
      fl%fbcy_pr(:, 3, :) = accp_ypencil(:, 2, :)

    end if
!----------------------------------------------------------------------------------------------------------
!   (gx, gy, gz) bc in y - direction
!----------------------------------------------------------------------------------------------------------
    if( dm%is_thermo) then
      if(dm%ibcy(1, 1) == IBC_INTERIOR) then
        call transpose_x_to_y(fl%gx,        apcc_ypencil, dm%apcc)
        call transpose_y_to_z(apcc_ypencil, apcc_zpencil, dm%apcc)

        do k = 1, dm%dpcc%zsz(3)
          apcc1_zpencil(:, :, k) = apcc_zpencil(:, :, dm%knc_sym(k))
        end do
        call transpose_z_to_y(apcc1_zpencil, apcc_ypencil, dm%apcc)

        fl%fbcy_gx(:, 1, :) = apcc_ypencil(:, 1, :)
        fl%fbcy_gx(:, 3, :) = apcc_ypencil(:, 2, :)
      end if

      ! uy, dm0-dm1, np
      fl%fbcy_gy(:, 1, :) = ZERO

      ! uz, dm0-dm1, nc
      if(dm%ibcy(1, 3) == IBC_INTERIOR) then
        call transpose_x_to_y(fl%gz,        accp_ypencil, dm%accp)
        call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%accp)

        do k = 1, dm%dpcc%zsz(3)
          accp1_zpencil(:, :, k) = accp_zpencil(:, :, dm%knc_sym(k))
        end do
        call transpose_z_to_y(accp1_zpencil, accp_ypencil, dm%accp)

        fl%fbcy_gz(:, 1, :) = accp_ypencil(:, 1, :)
        fl%fbcy_gz(:, 3, :) = accp_ypencil(:, 2, :)

      end if
    end if


    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine update_flow_bc_2dm_halo(dm0, fl0, dm1, fl1)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm0, dm1
    type(t_flow), intent(in)      :: fl0, fl1
    
    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: apcc_ypencil
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: apcc_zpencil
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: acpc_xpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: acpc_zpencil
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: accp_xpencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil
    real(WP), dimension( dm%dccc%ysz(1), dm%cpcc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%cpcc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil

!----------------------------------------------------------------------------------------------------------
!   default in x-pencil
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
!   qx, qy, qz, pres, bc in x - direction,  no overlap of values
!----------------------------------------------------------------------------------------------------------
    ! qx - x, dm0<->dm1, np
    if(dm1%ibcx(1, 1) == IBC_INTERIOR) then
      fl1%fbcx_qx(1, :, :) = fl0%qx(dm0%np_geo(1) - 1, :, :)
      fl1%fbcx_qx(3, :, :) = fl0%qx(dm0%np_geo(1) - 2, :, :)
    end if
    if(dm0%ibcx(2, 1) == IBC_INTERIOR) then
      fl0%fbcx_qx(2, :, :) = fl1%qx(2, :, :)
      fl0%fbcx_qx(4, :, :) = fl1%qx(3, :, :)
    end if
    ! qy - x, dm0<->dm1, nc
    if(dm1%ibcx(1, 2) == IBC_INTERIOR) then
      fl1%fbcx_qy(1, :, :) = fl0%qy(dm0%nc(1),     :, :)
      fl1%fbcx_qy(3, :, :) = fl0%qy(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx(2, 2) == IBC_INTERIOR) then
      fl0%fbcx_qy(2, :, :) = fl1%qy(1, :, :)
      fl0%fbcx_qy(4, :, :) = fl1%qy(2, :, :)
    end if
    ! qz - x, dm0<->dm1, nc
    if(dm1%ibcx(1, 3) == IBC_INTERIOR) then
      fl1%fbcx_qz(1, :, :) = fl0%qz(dm0%nc(1),     :, :)
      fl1%fbcx_qz(3, :, :) = fl0%qz(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx(2, 3) == IBC_INTERIOR) then
      fl0%fbcx_qz(2, :, :, m) = fl1%qz(1, :, :)
      fl0%fbcx_qz(4, :, :, m) = fl1%qz(2, :, :)
    end if
    ! p - x, dm0-dm1, nc
    if(dm1%ibcx(1, 4) == IBC_INTERIOR) then
      fl1%fbcx_pr(1, :, :) = fl0%pres(dm0%nc(1),     :, :)
      fl1%fbcx_pr(3, :, :) = fl0%pres(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx(2, m) == IBC_INTERIOR) then
      fl0%fbcx_pr(2, :, :) = fl1%pres(1, :, :)
      fl0%fbcx_pr(4, :, :) = fl1%pres(2, :, :)
    end if
!----------------------------------------------------------------------------------------------------------
!   qx, qy, qz, pres, bc in y - direction,  no overlap of values
!----------------------------------------------------------------------------------------------------------
    ! qx - y, dm0<->dm1, nc
    if(dm1%ibcy(1, 1) == IBC_INTERIOR) then
      call transpose_x_to_y(fl0%qx, apcc_ypencil, dm%apcc)
      fl1%fbcy_qx(:, 1, :) = apcc_ypencil(:, dm0%nc(2),     :)
      fl1%fbcy_qx(:, 3, :) = apcc_ypencil(:, dm0%nc(2) - 1, :)
    end if
    if(dm0%ibcy(2, 1) == IBC_INTERIOR) then
      call transpose_x_to_y(fl1%qx, apcc_ypencil, dm%apcc)
      fl0%fbcy_qx(:, 2, :) = apcc_ypencil(:, 1, :)
      fl0%fbcy_qx(:, 4, :) = apcc_ypencil(:, 2, :)
    end if
    ! qy - y, dm0<->dm1, np
    if(dm1%ibcy(1, 2) == IBC_INTERIOR) then
      call transpose_x_to_y(fl0%qy, acpc_ypencil, dm%acpc)
      fl1%fbcy_qy(:, 1, :) = acpc_ypencil(:, dm0%np_geo(2) - 1, :)
      fl1%fbcy_qy(:, 3, :) = acpc_ypencil(:, dm0%np_geo(2) - 2, :)
    end if
    if(dm0%ibcy(2, 2) == IBC_INTERIOR) then
      call transpose_x_to_y(fl1%qy, acpc_ypencil, dm%acpc)
      fl0%fbcy_qy(:, 2, :) = acpc_ypencil(:, 2, :)
      fl0%fbcy_qy(:, 4, :) = acpc_ypencil(:, 3, :)
    end if
    ! qz - y, dm0<->dm1, nc
    if(dm1%ibcy(1, 3) == IBC_INTERIOR) then
      call transpose_x_to_y(fl0%qz, accp_ypencil, dm%accp)
      fl1%fbcy_qz(:, 1, :) = accp_ypencil(:, dm0%nc(2),     :)
      fl1%fbcy_qz(:, 3, :) = accp_ypencil(:, dm0%nc(2) - 1, :)
    end if
    if(dm0%ibcy(2, 3) == IBC_INTERIOR) then
      call transpose_x_to_y(fl1%qz, accp_ypencil, dm%accp)
      fl0%fbcy_qz(:, 2, :) = accp_ypencil(:, 1, :)
      fl0%fbcy_qz(:, 4, :) = accp_ypencil(:, 2, :)
    end if
    ! p - y, dm0<->dm1, nc
    if(dm1%ibcy(1, 4) == IBC_INTERIOR) then
      call transpose_x_to_y(fl0%pres, accc_ypencil, dm%accc)
      fl1%fbcy_pr(:, 1, :) = accc_ypencil(:, dm0%nc(2),     :)
      fl1%fbcy_pr(:, 3, :) = accc_ypencil(:, dm0%nc(2) - 1, :)
    end if
    if(dm0%ibcy(2, 4) == IBC_INTERIOR) then
      call transpose_x_to_y(fl1%pres, accc_ypencil, dm%accc)
      fl0%fbcy_pr(:, 2, :) = accc_ypencil(:, 1, :)
      fl0%fbcy_pr(:, 4, :) = accc_ypencil(:, 2, :)
    end if

!----------------------------------------------------------------------------------------------------------
!   qx, qy, qz, pres, bc in z - direction,  no overlap of values
!----------------------------------------------------------------------------------------------------------
    ! qx - z, dm0<->dm1, nc
    if(dm1%ibcz(1, 1) == IBC_INTERIOR) then
      call transpose_x_to_y(fl0%qx, apcc_ypencil, dm%apcc)
      call transpose_y_to_z(apcc_ypencil, apcc_zpencil, dm%apcc)
      fl1%fbcz_qx(:, :, 1) = apcc_zpencil(:, :, dm0%nc(2)    )
      fl1%fbcz_qx(:, :, 3) = apcc_zpencil(:, :, dm0%nc(2) - 1)
    end if
    if(dm0%ibcz(2, 1) == IBC_INTERIOR) then
      call transpose_x_to_y(fl1%qx, apcc_ypencil, dm%apcc)
      call transpose_y_to_z(apcc_ypencil, apcc_zpencil, dm%apcc)
      fl0%fbcz_qx(:, :, 2) = apcc_zpencil(:, :, 1)
      fl0%fbcz_qx(:, :, 4) = apcc_zpencil(:, :, 2)
    end if
    ! qy - z, dm0-dm1, np
    if(dm1%ibcz(1, 2) == IBC_INTERIOR) then
      call transpose_x_to_y(fl0%qy, acpc_ypencil, dm%acpc)
      call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%acpc)
      fl1%fbcz_qy(:, :, 1) = acpc_zpencil(:, :, dm0%nc(2)    )
      fl1%fbcz_qy(:, :, 3) = acpc_zpencil(:, :, dm0%nc(2) - 1)
    end if
    if(dm0%ibcz(2, 2) == IBC_INTERIOR) then
      call transpose_x_to_y(fl1%qy, acpc_ypencil, dm%acpc)
      call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%acpc)
      fl0%fbcz_qy(:, :, 2) = acpc_zpencil(:, :, 1)
      fl0%fbcz_qy(:, :, 4) = acpc_zpencil(:, :, 2)
    end if
    ! qz - z, dm0-dm1, nc
    if(dm1%ibcz(1, 3) == IBC_INTERIOR) then
      call transpose_x_to_y(fl0%qz, accp_ypencil, dm%accp)
      call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%accp)
      fl1%fbcz_qz(:, :, 1) = accp_zpencil(:, :, dm0%np_geo(2) - 1)
      fl1%fbcz_qz(:, :, 3) = accp_zpencil(:, :, dm0%np_geo(2) - 2)
    end if
    if(dm0%ibcz(2, 3) == IBC_INTERIOR) then
      call transpose_x_to_y(fl1%qz, accp_ypencil, dm%accp)
      call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%accp)
      fl0%fbcz_qz(:, :, 2) = accp_zpencil(:, :, 2)
      fl0%fbcz_qz(:, :, 4) = accp_zpencil(:, :, 3)
    end if
    ! p - z, dm0-dm1, nc
    if(dm1%ibcz(1, 4) == IBC_INTERIOR) then
      call transpose_x_to_y(fl0%pres, accc_ypencil, dm%accc)
      call transpose_y_to_z(accc_ypencil, accc_zpencil, dm%accc)
      fl1%fbcz_pr(:, :, 1) = accc_zpencil(:, :, dm0%nc(2)    )
      fl1%fbcz_pr(:, :, 3) = accc_zpencil(:, :, dm0%nc(2) - 1)
    end if
    if(dm0%ibcz(2, 4) == IBC_INTERIOR) then
      call transpose_x_to_y(fl1%pres, accc_ypencil, dm%accc)
      call transpose_y_to_z(accc_ypencil, accc_zpencil, dm%accc)
      fl0%fbcz_pr(:, :, 2) = accc_zpencil(:, :, 1)
      fl0%fbcz_pr(:, :, 4) = accc_zpencil(:, :, 2)
    end if

!----------------------------------------------------------------------------------------------------------
!   thermal only
!----------------------------------------------------------------------------------------------------------
    ! gx - x, dm0<->dm1, np
    if(dm1%is_thermo .and. dm1%ibcx(1, 1) == IBC_INTERIOR) then
      if (dm0%is_thermo) then
        fl1%fbcx_gx(1, :, :) = fl0%gx(dm0%np_geo(1) - 1, :, :)
        fl1%fbcx_gx(3, :, :) = fl0%gx(dm0%np_geo(1) - 2, :, :)
      else
        fl1%fbcx_gx(1, :, :) = fl0%fbcx_qx(1, :, :) * tm1%fbcx_ftp(1, :, :)%d
        fl1%fbcx_gx(3, :, :) = fl0%fbcx_qx(3, :, :) * tm1%fbcx_ftp(3, :, :)%d
      end if
    end if
    if(dm0%is_thermo .and. dm0%ibcx(2, 1) == IBC_INTERIOR) then
      if (dm1%is_thermo) then
        fl1%fbcx_gx(2, :, :) = fl0%gx(2, :, :)
        fl1%fbcx_gx(4, :, :) = fl0%gx(3, :, :)
      else
        call Print_warning_msg("Error in setting up thermo b.c.")
      end if
    end if
    ! gy - x, dm0<->dm1, nc
    if(dm1%is_thermo .and. dm1%ibcx(1, 2) == IBC_INTERIOR) then
      if (dm0%is_thermo) then
        fl1%fbcx_gy(1, :, :) = fl0%gy(dm0%nc(1),     :, :)
        fl1%fbcx_gy(3, :, :) = fl0%gy(dm0%nc(1) - 1, :, :)
      else
        fl1%fbcx_gy(1, :, :) = fl0%fbcx_qy(1, :, :) * tm1%fbcx_ftp(1, :, :)%d
        fl1%fbcx_gy(3, :, :) = fl0%fbcx_qy(3, :, :) * tm1%fbcx_ftp(3, :, :)%d
      end if
    end if
    if(dm0%is_thermo .and. dm0%ibcx(2, 2) == IBC_INTERIOR) then
      if (dm1%is_thermo) then
        fl0%fbcx_gy(2, :, :) = fl1%gy(1, :, :)
        fl0%fbcx_gy(4, :, :) = fl1%gy(2, :, :)
      else
        call Print_warning_msg("Error in setting up thermo b.c.")
      end if
    end if
    !gz - x, dm0<->dm1, nc
    if(dm1%is_thermo .and. dm1%ibcx(1, 3) == IBC_INTERIOR) then
      if (dm0%is_thermo) then
        fl1%fbcx_gz(1, :, :) = fl0%gz(dm0%nc(1),     :, :)
        fl1%fbcx_gz(3, :, :) = fl0%gz(dm0%nc(1) - 1, :, :)
      else
        fl1%fbcx_gz(1, :, :) = fl0%fbcx_qz(1, :, :) * tm1%fbcx_ftp(1, :, :)%d
        fl1%fbcx_gz(3, :, :) = fl0%fbcx_qz(3, :, :) * tm1%fbcx_ftp(3, :, :)%d
      end if
    end if
    if(dm0%is_thermo .and. dm0%ibcx(2, 3) == IBC_INTERIOR) then
      if (dm1%is_thermo) then
        fl0%fbcx_gz(2, :, :) = fl1%gz(1, :, :)
        fl0%fbcx_gz(4, :, :) = fl1%gz(2, :, :)
      else
        call Print_warning_msg("Error in setting up thermo b.c.")
      end if
    end if
    ! gx - y, dm0<->dm1, nc
    if(dm1%is_thermo .and. dm1%ibcy(1, 1) == IBC_INTERIOR) then
      if (dm0%is_thermo) then
        call transpose_x_to_y(fl0%gx, apcc_ypencil, dm%apcc)
        fl1%fbcy_gx(:, 1, :) = apcc_ypencil(:, dm0%nc(2),     :)
        fl1%fbcy_gx(:, 3, :) = apcc_ypencil(:, dm0%nc(2) - 1, :)
      else
        fl1%fbcy_gx(:, 1, :) = fl1%fbcy_qx(:, 1, :) * tm1%fbcy_ftp(:, 1, :)%d
        fl1%fbcy_gx(:, 3, :) = fl1%fbcy_qx(:, 3, :) * tm1%fbcy_ftp(:, 1, :)%d
      end if
    end if
    if(dm0%is_thermo .and. dm0%ibcy(2, 1) == IBC_INTERIOR) then
      if (dm1%is_thermo) then
        call transpose_x_to_y(fl1%gx, apcc_ypencil, dm%apcc)
        fl0%fbcy_gx(:, 2, :) = apcc_ypencil(:, 1, :)
        fl0%fbcy_gx(:, 4, :) = apcc_ypencil(:, 2, :)
      else
        call Print_warning_msg("Error in setting up thermo b.c.")
      end if
    end if
    ! gy - y, dm0<->dm1, np
    if(dm1%is_thermo .and. dm1%ibcy(1, 2) == IBC_INTERIOR) then
      if (dm0%is_thermo) then
        call transpose_x_to_y(fl0%gy, acpc_ypencil, dm%acpc)
        fl1%fbcy_gy(:, 1, :) = acpc_ypencil(:, dm0%np_geo(2) - 1, :)
        fl1%fbcy_gy(:, 3, :) = acpc_ypencil(:, dm0%np_geo(2) - 2, :)
      else
        fl1%fbcy_gy(:, 1, :) = fl1%fbcy_qy(:, 1, :) * tm1%fbcy_ftp(:, 1, :)%d
        fl1%fbcy_gy(:, 3, :) = fl1%fbcy_qy(:, 3, :) * tm1%fbcy_ftp(:, 1, :)%d
      end if
    end if
    if(dm0%is_thermo .and. dm0%ibcy(2, 2) == IBC_INTERIOR) then
      if (dm1%is_thermo) then
        call transpose_x_to_y(fl1%gy, acpc_ypencil, dm%acpc)
        fl0%fbcy_gy(:, 2, :) = acpc_ypencil(:, 2, :)
        fl0%fbcy_gy(:, 4, :) = acpc_ypencil(:, 3, :)
      else
        call Print_warning_msg("Error in setting up thermo b.c.")
      end if
    end if
    ! gz - y, dm0<->dm1, nc
    if(dm1%is_thermo .and. dm1%ibcy(1, 3) == IBC_INTERIOR) then
      if (dm0%is_thermo) then
        call transpose_x_to_y(fl0%gz, accp_ypencil, dm%accp)
        fl1%fbcy_gz(:, 1, :) = accp_ypencil(:, dm0%nc(2),     :)
        fl1%fbcy_gz(:, 3, :) = accp_ypencil(:, dm0%nc(2) - 1, :)
      else 
        fl1%fbcy_gz(:, 1, :) = fl1%fbcy_qz(:, 1, :) * tm1%fbcy_ftp(:, 1, :)%d
        fl1%fbcy_gz(:, 3, :) = fl1%fbcy_qz(:, 3, :) * tm1%fbcy_ftp(:, 1, :)%d
      end if
    end if
    if(dm0%is_thermo .and. dm0%ibcy(2, 3) == IBC_INTERIOR) then
      if (dm1%is_thermo) then
        call transpose_x_to_y(fl1%gz, accp_ypencil, dm%accp)
        fl0%fbcy_gz(:, 2, :) = accp_ypencil(:, 1, :)
        fl0%fbcy_gz(:, 4, :) = accp_ypencil(:, 2, :)
      else
        call Print_warning_msg("Error in setting up thermo b.c.")
      end if
    end if
    ! gx - z, dm0<->dm1, nc
    if(dm1%is_thermo .and. dm1%ibcz(1, 1) == IBC_INTERIOR) then
      if (dm0%is_thermo) then
        call transpose_x_to_y(fl0%gx, apcc_ypencil, dm%apcc)
        call transpose_y_to_z(apcc_ypencil, apcc_zpencil, dm%apcc)
        fl1%fbcz_gx(:, :, 1) = apcc_zpencil(:, :, dm0%nc(2)    )
        fl1%fbcz_gx(:, :, 3) = apcc_zpencil(:, :, dm0%nc(2) - 1)
      else
        fl1%fbcz_gx(:, :, 1) = fl1%fbcz_qx(:, :, 1) * tm1%fbcz_ftp(:, :, 1)%d
        fl1%fbcz_gx(:, :, 3) = fl1%fbcz_qx(:, :, 3) * tm1%fbcz_ftp(:, :, 1)%d
      end if
    end if
    if(dm0%is_thermo .and. dm0%ibcz(2, 1) == IBC_INTERIOR) then
      if (dm1%is_thermo) then
        call transpose_x_to_y(fl1%gx, apcc_ypencil, dm%apcc)
        call transpose_y_to_z(apcc_ypencil, apcc_zpencil, dm%apcc)
        fl0%fbcz_gx(:, :, 2) = apcc_zpencil(:, :, 1)
        fl0%fbcz_gx(:, :, 4) = apcc_zpencil(:, :, 2)
      else 
        call Print_warning_msg("Error in setting up thermo b.c.")
      end if
    end if
    ! gy - z, dm0-dm1, np
    if(dm1%is_thermo .and. dm1%ibcz(1, 2) == IBC_INTERIOR) then
      if (dm0%is_thermo) then
        call transpose_x_to_y(fl0%gy, acpc_ypencil, dm%acpc)
        call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%acpc)
        fl1%fbcz_gy(:, :, 1) = acpc_zpencil(:, :, dm0%nc(2)    )
        fl1%fbcz_gy(:, :, 3) = acpc_zpencil(:, :, dm0%nc(2) - 1)
      else
        fl1%fbcz_gy(:, :, 1) = fl1%fbcz_qy(:, :, 1) * tm1%fbcz_ftp(:, :, 1)%d
        fl1%fbcz_gy(:, :, 3) = fl1%fbcz_qy(:, :, 3) * tm1%fbcz_ftp(:, :, 1)%d
      end if
    end if
    if(dm0%is_thermo .and. dm0%ibcz(2, 2) == IBC_INTERIOR) then
      if (dm1%is_thermo) then
        call transpose_x_to_y(fl1%qy, acpc_ypencil, dm%acpc)
        call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%acpc)
        fl0%fbcz_qy(:, :, 2) = acpc_zpencil(:, :, 1)
        fl0%fbcz_qy(:, :, 4) = acpc_zpencil(:, :, 2)
      else
        call Print_warning_msg("Error in setting up thermo b.c.")
      end if
    end if
    ! gz - z, dm0-dm1, nc
    if(dm1%is_thermo .and. dm1%ibcz(1, 3) == IBC_INTERIOR) then
      if (dm0%is_thermo) then
        call transpose_x_to_y(fl0%gz, accp_ypencil, dm%accp)
        call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%accp)
        fl1%fbcz_gz(:, :, 1) = accp_zpencil(:, :, dm0%np_geo(2) - 1)
        fl1%fbcz_gz(:, :, 3) = accp_zpencil(:, :, dm0%np_geo(2) - 2)
      else
        fl1%fbcz_gz(:, :, 1) = fl1%fbcz_qz(:, :, 1) * tm1%fbcz_ftp(:, :, 1)%d
        fl1%fbcz_gz(:, :, 3) = fl1%fbcz_qz(:, :, 3) * tm1%fbcz_ftp(:, :, 1)%d
      end if
    end if
    if(dm0%is_thermo .and. dm0%ibcz(2, 3) == IBC_INTERIOR) then
      if (dm1%is_thermo) then
        call transpose_x_to_y(fl1%qz, accp_ypencil, dm%accp)
        call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%accp)
        fl0%fbcz_qz(:, :, 2) = accp_zpencil(:, :, 2)
        fl0%fbcz_qz(:, :, 4) = accp_zpencil(:, :, 3)
      else 
        call Print_warning_msg("Error in setting up thermo b.c.")
      end if
    end if
!----------------------------------------------------------------------------------------------------------
!   clyindrical only
!----------------------------------------------------------------------------------------------------------
   if(dm1%icoordinate == ICYLINDRICAL .or. dm0%icoordinate == ICYLINDRICAL) then
      ! uy = qr/r - y, dm0-dm1, np
      if(dm1%ibcy(1, 2) == IBC_INTERIOR) then
        acpc_xpencil = fl0%qy
        call multiple_cylindrical_rn(acpc_xpencil, dm%dcpc, dm%rpi, 1, IPENCIL(1))
        call transpose_x_to_y(acpc_xpencil, acpc_ypencil, dm%acpc)
        fl1%fbcy_uy(:, 1, :) = acpc_ypencil(:, dm0%np_geo(2) - 1, :)
        fl1%fbcy_uy(:, 3, :) = acpc_ypencil(:, dm0%np_geo(2) - 2, :)
      end if
      if(dm0%ibcy(2, 2) == IBC_INTERIOR) then
        acpc_xpencil = fl1%qy
        call multiple_cylindrical_rn(acpc_xpencil, dm%dcpc, dm%rpi, 1, IPENCIL(1))
        call transpose_x_to_y(acpc_xpencil, acpc_ypencil, dm%acpc)
        fl0%fbcy_uy(:, 2, :) = acpc_ypencil(:, 2, :)
        fl0%fbcy_uy(:, 4, :) = acpc_ypencil(:, 3, :)
      end if
      ! uz = qz/r - y, dm0-dm1, nc
      if(dm1%ibcy(1, 3) == IBC_INTERIOR) then
        accp_xpencil = fl0%qz
        call multiple_cylindrical_rn(accp_xpencil, dm%dccp, dm%rci, 1, IPENCIL(1))
        call transpose_x_to_y(accp_xpencil, accp_ypencil, dm%accp)
        fl1%fbcy_uz(:, 1, :) = accp_ypencil(:, dm0%nc(2),     :)
        fl1%fbcy_uz(:, 3, :) = accp_ypencil(:, dm0%nc(2) - 1, :)
      end if
      if(dm0%ibcy(2, 3) == IBC_INTERIOR) then
        accp_xpencil = fl1%qz
        call multiple_cylindrical_rn(accp_xpencil, dm%dccp, dm%rci, 1, IPENCIL(1))
        call transpose_x_to_y(accp_xpencil, accp_ypencil, dm%accp)
        fl0%fbcy_uz(:, 2, :) = accp_ypencil(:, 1, :)
        fl0%fbcy_uz(:, 4, :) = accp_ypencil(:, 2, :)
      end if
      ! uy = qr/r - z, dm0-dm1, np
      if(dm1%ibcz(1, 2) == IBC_INTERIOR) then
        acpc_xpencil = fl0%qy
        call multiple_cylindrical_rn(acpc_xpencil, dm%dcpc, dm%rpi, 1, IPENCIL(1))
        call transpose_x_to_y(acpc_xpencil, acpc_ypencil, dm%acpc)
        call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%acpc)
        fl1%fbcz_uy(:, :, 1) = acpc_zpencil(:, :, dm0%nc(2)    )
        fl1%fbcz_uy(:, :, 3) = acpc_zpencil(:, :, dm0%nc(2) - 1)
      end if
      if(dm0%ibcz(2, 2) == IBC_INTERIOR) then
        acpc_xpencil = fl1%qy
        call multiple_cylindrical_rn(acpc_xpencil, dm%dcpc, dm%rpi, 1, IPENCIL(1))
        call transpose_x_to_y(acpc_xpencil, acpc_ypencil, dm%acpc)
        call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%acpc)
        fl0%fbcz_uy(:, :, 2) = acpc_zpencil(:, :, 1)
        fl0%fbcz_uy(:, :, 4) = acpc_zpencil(:, :, 2)
      end if

      ! uz = qz/r - z, dm0-dm1, nc
      if(dm1%ibcz(1, 3) == IBC_INTERIOR) then
        accp_xpencil = fl0%qz
        call multiple_cylindrical_rn(accp_xpencil, dm%dccp, dm%rci, 1, IPENCIL(1))
        call transpose_x_to_y(accp_xpencil, accp_ypencil, dm%accp)
        call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%accp)
        fl1%fbcz_uz(:, :, 1) = accp_zpencil(:, :, dm0%np_geo(2) - 1)
        fl1%fbcz_uz(:, :, 3) = accp_zpencil(:, :, dm0%np_geo(2) - 2)
      end if
      if(dm0%ibcz(2, 3) == IBC_INTERIOR) then
        accp_xpencil = fl1%qz
        call multiple_cylindrical_rn(accp_xpencil, dm%dccp, dm%rci, 1, IPENCIL(1))
        call transpose_x_to_y(accp_xpencil, accp_ypencil, dm%accp)
        call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%accp)
        fl0%fbcz_uz(:, :, 2) = accp_zpencil(:, :, 2)
        fl0%fbcz_uz(:, :, 4) = accp_zpencil(:, :, 3)
      end if
    end if
    
    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine update_thermo_bc_2dm_halo(dm0, tm0, dm1, tm1)
    use parameters_constant_mod
    use thermo_info_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in)    :: dm0, dm1
    type(t_thermo), intent(inout) :: tm0, tm1
    
    integer :: m, i, j, k
    real(WP), dimension( dm%dccc%ysz(1), dm%cpcc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%cpcc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
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
          call ftp_refresh_thermal_properties_from_T_undim(fl0%fbcx_ftp(2, j, k)%t)
          call ftp_refresh_thermal_properties_from_T_undim(fl0%fbcx_ftp(4, j, k)%t)
        end do
      end do
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in y - direction
!----------------------------------------------------------------------------------------------------------
    ! thermal field, dm0-dm1
    if(dm1%ibcy(1, 5) == IBC_INTERIOR) then
      call transpose_x_to_y(tm0%tTemp, accc_ypencil, dm%dccc)
      tm1%fbcy_ftp(:, 1, :)%t = accc_ypencil(:, dm0%nc(1),     :)
      tm1%fbcy_ftp(:, 3, :)%t = accc_ypencil(:, dm0%nc(1) - 1, :)
      do k = 1, dm1%dccc%ysz(3)
        do i = 1, dm1%dccc%ysz(1)
          call ftp_refresh_thermal_properties_from_T_undim(fl1%fbcy_ftp(i, 1, k)%t)
          call ftp_refresh_thermal_properties_from_T_undim(fl1%fbcy_ftp(i, 3, k)%t)
        end do
      end do
    end if

    if(dm0%ibcy(2, 5) == IBC_INTERIOR) then
      call transpose_x_to_y(tm1%tTemp, accc_ypencil, dm%dccc)
      tm0%fbcy_ftp(:, 2, :)%t = accc_ypencil(:, 1, :)
      tm0%fbcy_ftp(:, 4, :)%t = accc_ypencil(:, 2, :)
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
      call transpose_x_to_y(tm0%tTemp, accc_ypencil, dm%dccc)
      call transpose_y_to_z(accc_ypencil, accc_zpencil, dm%dccc)
      tm1%fbcz_ftp(:, :, 1)%t = accc_zpencil(:, :, dm0%nc(1)    )
      tm1%fbcz_ftp(:, :, 3)%t = accc_zpencil(:, :, dm0%nc(1) - 1)
      do j = 1, dm1%dccc%zsz(2)
        do i = 1, dm1%dccc%zsz(1)
          call ftp_refresh_thermal_properties_from_T_undim(tm1%fbcz_ftp(i, j, 1)%t)
          call ftp_refresh_thermal_properties_from_T_undim(tm1%fbcz_ftp(i, j, 3)%t)
        end do
      end do
    end if

    if(dm0%ibcz(2, 5) == IBC_INTERIOR) then
      call transpose_x_to_y(tm1%tTemp, accc_ypencil, dm%dccc)
      call transpose_y_to_z(accc_ypencil, accc_zpencil, dm%dccc)
      tm0%fbcz_ftp(:, :, 2)%t = accc_zpencil(:, :, 1)
      tm0%fbcz_ftp(:, :, 4)%t = accc_zpencil(:, :, 2)
      do j = 1, dm0%dccc%zsz(2)
        do i = 1, dm0%dccc%zsz(1)
          call ftp_refresh_thermal_properties_from_T_undim(tm0%fbcz_ftp(i, j, 2)%t)
          call ftp_refresh_thermal_properties_from_T_undim(tm0%fbcz_ftp(i, j, 4)%t)
        end do
      end do
    end if

    return
  end subroutine


!==========================================================================================================
!==========================================================================================================
  subroutine update_thermo_bc_1dm_halo(dm, tm)
    use parameters_constant_mod
    use thermo_info_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in)    :: dm
    type(t_thermo), intent(inout) :: tm
    
    integer :: i, k
    real(WP), dimension( dm%dccc%ysz(1), dm%cpcc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%cpcc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
    real(WP), dimension( dm%dccc%zsz(1), dm%cpcc%zsz(2), dm%dccc%zsz(3) ) :: accc1_zpencil

    if(dm%icase /= ICASE_PIPE) return
    if(.not. dm%is_thermo) return
!----------------------------------------------------------------------------------------------------------
!   all in z-pencil
!   no overlap of values
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
!   bc in y - direction
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy(1, 5) == IBC_INTERIOR) then
      call transpose_x_to_y(tm%tTemp, accc_ypencil, dm%dccc)
      call transpose_y_to_z(accc_ypencil, accc_zpencil, dm%accc)
      do k = 1, dm%dccc%zsz(3)
        accc1_zpencil(:, :, k) = accc_zpencil(:, :, dm%knc_sym(k))
      end do
      call transpose_z_to_y(accc1_zpencil, accc_ypencil, dm%apcc)

      tm%fbcy_ftp(:, 1, :)%t = accc_ypencil(:, 1, :)
      tm%fbcy_ftp(:, 3, :)%t = accc_ypencil(:, 2, :)
      do k = 1, dm%dccc%ysz(3)
        do i = 1, dm%dccc%ysz(1)
          call ftp_refresh_thermal_properties_from_T_undim(tm%fbcy_ftp(i, 1, k)%t)
          call ftp_refresh_thermal_properties_from_T_undim(tm%fbcy_ftp(i, 3, k)%t)
        end do
      end do
    end if

    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine apply_gxgygz_bc_geo (fl, tm) ! 
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_flow), intent( inout) :: fl
    type(t_thermo), intent(in)   :: tm

    if(.not. dm%is_thermo) return

!----------------------------------------------------------------------------------------------------------
!   get bc gx, gy, gz (at bc not cell centre)
!----------------------------------------------------------------------------------------------------------
    fl%fbcx_gx(:, :, :) = fl%fbcx_qx(:, :, :) * tm%fbcx_ftp(:, :, :)%d
    fl%fbcx_gy(:, :, :) = fl%fbcx_qy(:, :, :) * tm%fbcx_ftp(:, :, :)%d
    fl%fbcx_gz(:, :, :) = fl%fbcx_qz(:, :, :) * tm%fbcx_ftp(:, :, :)%d

    fl%fbcy_gx(:, :, :) = fl%fbcy_qx(:, :, :) * tm%fbcy_ftp(:, :, :)%d
    fl%fbcy_gy(:, :, :) = fl%fbcy_qy(:, :, :) * tm%fbcy_ftp(:, :, :)%d
    fl%fbcy_gz(:, :, :) = fl%fbcy_qz(:, :, :) * tm%fbcy_ftp(:, :, :)%d

    fl%fbcz_gx(:, :, :) = fl%fbcz_qx(:, :, :) * tm%fbcz_ftp(:, :, :)%d
    fl%fbcz_gy(:, :, :) = fl%fbcz_qy(:, :, :) * tm%fbcz_ftp(:, :, :)%d
    fl%fbcz_gz(:, :, :) = fl%fbcz_qz(:, :, :) * tm%fbcz_ftp(:, :, :)%d

    if(dm%icoordinate == ICYLINDRICAL) then 
      fl%fbcy_gyr(:, :, :) = fl%fbcy_qyr(:, :, :) * tm%fbcy_ftp(:, :, :)%d
      fl%fbcz_gyr(:, :, :) = fl%fbcz_qyr(:, :, :) * tm%fbcz_ftp(:, :, :)%d
      fl%fbcy_gzr(:, :, :) = fl%fbcy_qzr(:, :, :) * tm%fbcy_ftp(:, :, :)%d
      fl%fbcz_gzr(:, :, :) = fl%fbcz_qzr(:, :, :) * tm%fbcz_ftp(:, :, :)%d
    end if

    return
  end subroutine

end module
