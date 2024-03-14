module boundary_conditions_mod
  use precision_mod
  use parameters_constant_mod
  character(12) :: filename(5)
  
  private :: refresh_bc_type
  public  :: configure_bc_type

  private :: map_bc_1d_uprofile  ! for bc with a given profile

  private :: apply_x_bc_geo
  private :: apply_y_bc_geo
  private :: apply_z_bc_geo

  private :: apply_qxqyqzpr_bc_geo   ! constant given bc, apply once only  
  public  :: apply_gxgygz_bc_geo ! constant given bc, apply once only

  public  :: buildup_thermo_bc_geo ! constant given bc, apply once only
  public  :: buildup_flow_bc_geo   ! constant given bc, apply once only
  
  private :: get_cyindrical_x_axial_symmetry_bc
  public  :: update_thermo_bc_1dm_halo ! for pipe only, update every NS
  public  :: update_flow_bc_1dm_halo   ! for pipe only, update every NS 
  
  private :: apply_x_bc_2dm_halo
  private :: apply_y_bc_2dm_halo
  private :: apply_z_bc_2dm_halo
  public  :: update_flow_bc_2dm_halo   ! for multiple domains only, update every NS 
  public  :: update_thermo_bc_2dm_halo ! for multiple domains only, update every NS

  public  :: get_dirichlet_geo_bcx
  public  :: get_dirichlet_geo_bcy
  public  :: get_dirichlet_geo_bcz
  
contains

  subroutine refresh_bc_type(bc_nominal, bc_real)
    implicit none 
    integer, intent(in) :: bc_nominal(2, NBC)
    integer, intent(inout) :: bc_real(2, NBC)
    integer :: n, m

    do n = 1, 2
      do m = 1, NBC
        if (bc_nominal(n, m) == IBC_PROFILE1D)   then
          bc_real(n, m) = IBC_DIRICHLET
        else if (bc_nominal(n, m) == IBC_TURBGEN  .or. &
                 bc_nominal(n, m) == IBC_DATABASE )   then
          if(m /=5) then
            bc_real(n, m) = IBC_INTERIOR  ! for u, v, w, p
          else 
            bc_real(n, m) = IBC_DIRICHLET ! for temperature, default is no incoming thermal flow, check
          end if
        else if (bc_nominal(n, m) == IBC_CONVECTIVE)   then
          bc_real(n, m) = IBC_INTRPL
        else
          bc_real(n, m) = bc_nominal(n, m)   
        end if
      end do
    end do

    return
  end subroutine 


!==========================================================================================================
!==========================================================================================================
  subroutine configure_bc_type(dm)
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    integer :: m, n

!----------------------------------------------------------------------------------------------------------
! to set up real bc for calculation from given nominal b.c.
!----------------------------------------------------------------------------------------------------------
    call refresh_bc_type(dm%ibcx_nominal, dm%ibcx)
    call refresh_bc_type(dm%ibcy_nominal, dm%ibcy)
    call refresh_bc_type(dm%ibcz_nominal, dm%ibcz)
!----------------------------------------------------------------------------------------------------------
! treatment for cylindrical coordinates
!----------------------------------------------------------------------------------------------------------
    if(dm%icoordinate == ICYLINDRICAL) then
      dm%ibcz(:, :) = IBC_PERIODIC
    end if

    if (dm%icase == ICASE_PIPE) then
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
    use udf_type_mod
    use thermo_info_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(inout) :: tm
    real(WP) :: var1y(1:dm%np(2))

    integer :: i, j, k, n, jj, ny

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
            tm%fbcx_ftp(n, j, k)%t = dm%fbcx_const(n, 5) / fluidparam%ftp0ref%t
            call ftp_refresh_thermal_properties_from_T_undim( tm%fbcx_ftp(n,   j, k))
            tm%fbcx_ftp(n+2, j, k) = tm%fbcx_ftp(n, j, k)
          else if (dm%ibcx(n, 5) == IBC_NEUMANN) then
            ! dimensional heat flux (k*dT/dx) --> undimensional heat flux (k*dT/dx)
            tm%fbcx_heatflux(n,   j, k) = dm%fbcx_const(n, 5) * tm%ref_l0 / fluidparam%ftp0ref%k / fluidparam%ftp0ref%t 
            tm%fbcx_heatflux(n+2, j, k) = tm%fbcx_heatflux(n, j, k)
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
            call ftp_refresh_thermal_properties_from_T_undim( tm%fbcy_ftp(i, n,   k))
            call ftp_refresh_thermal_properties_from_T_undim( tm%fbcy_ftp(i, n+2, k))

          else if (dm%ibcy(n, 5) == IBC_NEUMANN) then
            ! dimensional heat flux (k*dT/dx) --> undimensional heat flux (k*dT/dx)
            tm%fbcy_heatflux(i, n,   k) = dm%fbcy_const(n, 5) * tm%ref_l0 / fluidparam%ftp0ref%k / fluidparam%ftp0ref%t 
            tm%fbcy_heatflux(i, n+2, k) = tm%fbcy_heatflux(i, n, k)
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
            call ftp_refresh_thermal_properties_from_T_undim( tm%fbcz_ftp(i, j, n  ))
            call ftp_refresh_thermal_properties_from_T_undim( tm%fbcz_ftp(i, j, n+2))

          else if (dm%ibcz(n, 5) == IBC_NEUMANN) then
            ! dimensional heat flux (k*dT/dx) --> undimensional heat flux (k*dT/dx)
            tm%fbcz_heatflux(i, j, n  ) = dm%fbcz_const(n, 5) * tm%ref_l0 / fluidparam%ftp0ref%k / fluidparam%ftp0ref%t 
            tm%fbcz_heatflux(i, j, n+2) = tm%fbcz_heatflux(i, j, n)
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
          call ftp_refresh_thermal_properties_from_T_undim( tm%fbcx_ftp(1, j, k))
          call ftp_refresh_thermal_properties_from_T_undim( tm%fbcx_ftp(3, j, k))
        end do
      end do
    end if



    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine apply_x_bc_profile(fbcx, var1y, jst, ri)
    use precision_mod
    implicit none
    real(WP), intent(inout) :: fbcx(:, :, :)
    real(WP), intent(in)    :: var1y(:)
    integer,  intent(in), optional :: jst
    real(WP), intent(in), optional :: ri(:)

    integer :: k, j, jj

    jj = 0
    do k = 1, size(fbcx, 3)
      do j = 1, size(fbcx, 2)
        jj = jst + j - 1
        if(present(ri)) then
            fbcx(1, j, k) = var1y(jj) / ri(jj)
        else
          fbcx(1, j, k) = var1y(jj)
        end if
        fbcx(3, j, k) = fbcx(1, j, k)
      end do
    end do

    return
  end subroutine
!==========================================================================================================
  subroutine apply_x_bc_geo(fbcx, fbcx_const, ri, jst)
    implicit none
    real(WP), intent(inout) :: fbcx(:, :, :)
    real(WP), intent(in)    :: fbcx_const(4)
    integer,  intent(in), optional :: jst
    real(WP), intent(in), optional :: ri(:)

    integer :: k, j, n, jj
    
    jj = 0
    do k = 1, size(fbcx, 3)
      do j = 1, size(fbcx, 2)
        if(present(jst)) jj = jst + j - 1
        do n = 1, 4
          if(present(ri)) then
            fbcx(n,   j, k) =  fbcx_const(n) / ri(jj)
          else
            fbcx(n,   j, k) =  fbcx_const(n)
          end if
        end do
      end do
    end do

    return
  end subroutine
!==========================================================================================================
  subroutine apply_y_bc_geo(fbcy, fbcy_const, ri)
    implicit none
    real(WP), intent(inout) :: fbcy(:, :, :)
    real(WP), intent(in)    :: fbcy_const(2)
    real(WP), intent(in), optional :: ri(:)

    integer :: k, i, n
    real(WP) :: ri_new(2)

    if(present(ri)) then
      ri_new(1) = ri(1)
      ri_new(2) = ri(size(ri))
    else
      ri_new = ONE
    end if
    

    do k = 1, size(fbcy, 3)
      do i = 1, size(fbcy, 1)
        do n = 1, 2
          fbcy(i, n,   k) =  fbcy_const(n) / ri_new(n) ! check  rpi=maxp?
          fbcy(i, n+2, k) =  fbcy(i, n, k)
        end do
      end do
    end do

    return
  end subroutine
!==========================================================================================================
  subroutine apply_z_bc_geo(fbcz, fbcz_const, ri, jst)
    implicit none
    real(WP), intent(inout) :: fbcz(:, :, :)
    real(WP), intent(in)    :: fbcz_const(2)
    integer,  intent(in), optional :: jst
    real(WP), intent(in), optional :: ri(:)


    integer :: i, j, n, jj

    jj = 0
    do j = 1, 1, size(fbcz, 2)
      if(present(jst)) jj = jst + j - 1
      do i = 1, size(fbcz, 1)
        do n = 1, 2
        if(present(ri)) then
            fbcz(i, j, n  ) =  fbcz_const(n) / ri(jj)
          else
            fbcz(i, j, n  ) =  fbcz_const(n)
          end if
          fbcz(i, j, n+2) =  fbcz(i, j, n)
        end do
      end do
    end do

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
  subroutine apply_qxqyqzpr_bc_geo (dm, fl) ! apply once only
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent( in)   :: dm
    type(t_flow), intent(inout)   :: fl

    real(WP) :: var1y(1:dm%np(2))
    
    integer :: ny
!==========================================================================================================
! to build up bc with constant values
! -3-1-||||-2-4
! for constant bc, 3=1= geometric bc, side 1;
!                  2=4= geometric bc, side 2 
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
! x-bc in x-pencil, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
    call apply_x_bc_geo(fl%fbcx_qx, dm%fbcx_const(:, 1))
    call apply_x_bc_geo(fl%fbcx_qy, dm%fbcx_const(:, 2), dm%rpi, dm%dcpc%xst(2))
    call apply_x_bc_geo(fl%fbcx_qz, dm%fbcx_const(:, 3), dm%rci, dm%dccp%xst(2))
    call apply_x_bc_geo(fl%fbcx_pr, dm%fbcx_const(:, 4))
!----------------------------------------------------------------------------------------------------------
! y-bc in y-pencil, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
    call apply_y_bc_geo(fl%fbcy_qx, dm%fbcy_const(:, 1))
    call apply_y_bc_geo(fl%fbcy_qy, dm%fbcy_const(:, 2), dm%rpi)
    call apply_y_bc_geo(fl%fbcy_qz, dm%fbcy_const(:, 3), dm%rpi) ! geo_bc, rpi, not rci
    call apply_y_bc_geo(fl%fbcy_pr, dm%fbcy_const(:, 4))
!----------------------------------------------------------------------------------------------------------
! z-bc in z-pencil, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
    call apply_z_bc_geo(fl%fbcz_qx, dm%fbcz_const(:, 1))
    call apply_z_bc_geo(fl%fbcz_qy, dm%fbcz_const(:, 2), dm%rpi, dm%dcpc%zst(2))
    call apply_z_bc_geo(fl%fbcz_qz, dm%fbcz_const(:, 3), dm%rci, dm%dccp%zst(2))
    call apply_z_bc_geo(fl%fbcz_pr, dm%fbcz_const(:, 4))
!==========================================================================================================
! y bc for primary variables, uy, uz
!==========================================================================================================
    if(dm%icoordinate == ICYLINDRICAL) then
!----------------------------------------------------------------------------------------------------------
! y-bc in y-pencil, qyr = qy/r = uy; qzr = qz/r = uz
!----------------------------------------------------------------------------------------------------------
      call apply_y_bc_geo(fl%fbcy_qyr, dm%fbcy_const(:, 2))
      call apply_y_bc_geo(fl%fbcy_qzr, dm%fbcy_const(:, 3))
!----------------------------------------------------------------------------------------------------------
! z-bc in z-pencil, qyr = qy/r = uy; qzr = qz/r = uz
!----------------------------------------------------------------------------------------------------------
      call apply_z_bc_geo(fl%fbcz_qyr, dm%fbcz_const(:, 2))
      call apply_z_bc_geo(fl%fbcz_qzr, dm%fbcz_const(:, 3))
    end if

!==========================================================================================================
! to build up bc for var(x_const, y, z)
!==========================================================================================================!----------------------------------------------------------------------------------------------------------
! x-bc1, qx(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 1) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%nc(2)
      filename(1) = 'pf1d_u1y.dat' !(undim)
      call map_bc_1d_uprofile( filename(1), ny, dm%yc, var1y(1:ny) )
      call apply_x_bc_profile(fl%fbcx_qx, var1y, dm%dpcc%xst(2))
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, qy(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 2) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%np(2)
      filename(2) = 'pf1d_v1y.dat' !(undim)
      call map_bc_1d_uprofile( filename(2), ny, dm%yp, var1y(1:ny) )
      call apply_x_bc_profile(fl%fbcx_qy, var1y, dm%dcpc%xst(2), dm%rpi)
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, qz(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 3) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%nc(2)
      filename(3) = 'pf1d_w1y.dat' !(undim)
      call map_bc_1d_uprofile( filename(3), ny, dm%yc, var1y(1:ny) )
      call apply_x_bc_profile(fl%fbcx_qz, var1y, dm%dccp%xst(2), dm%rci)
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, pr(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 4) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%nc(2)
      filename(4) = 'pf1d_p1y.dat' !(undim)
      call map_bc_1d_uprofile( filename(4), ny, dm%yc, var1y(1:ny) )
      call apply_x_bc_profile(fl%fbcx_pr, var1y, dm%dccc%xst(2))
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

    call apply_qxqyqzpr_bc_geo(dm, fl)
    
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
      
    end if

    return
  end subroutine 
!==========================================================================================================
!==========================================================================================================
  subroutine get_cyindrical_x_axial_symmetry_bc(var_xpencil, fbcy, ksym, dtmp)
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), intent(in) :: var_xpencil(:, :, :)
    real(WP), intent(inout) :: fbcy(:, :, :)
    integer, intent(in) :: ksym(:)

    real(WP), dimension( dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3) ) :: var_ypencil
    real(WP), dimension( dtmp%zsz(1), dtmp%zsz(2), dtmp%zsz(3) ) :: var_zpencil
    real(WP), dimension( dtmp%zsz(1), dtmp%zsz(2), dtmp%zsz(3) ) :: var_zpencil1

    integer :: k
!----------------------------------------------------------------------------------------------------------
!   all in z-pencil
!   no overlap of values
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(var_xpencil, var_ypencil, dtmp)
    call transpose_y_to_z(var_ypencil, var_zpencil, dtmp)

    do k = 1, dtmp%zsz(3)
      var_zpencil1(:, :, k) = var_zpencil(:, :, ksym(k))
    end do
    call transpose_z_to_y(var_zpencil1, var_ypencil, dtmp)

    fbcy(:, 1, :) = var_ypencil(:, 1, :)
    fbcy(:, 3, :) = var_ypencil(:, 2, :)

    return
  end subroutine 


!==========================================================================================================
!==========================================================================================================
  subroutine update_flow_bc_1dm_halo(dm, fl) ! for cylindrical only
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    type(t_flow), intent(inout)    :: fl

    if(dm%icase /= ICASE_PIPE) return

!----------------------------------------------------------------------------------------------------------
!   (qx, qy, qz, pr) bc in y - direction, interior
!----------------------------------------------------------------------------------------------------------
    ! ux, dm0-dm1, nc
    if(dm%ibcy(1, 1) == IBC_INTERIOR) then
      call get_cyindrical_x_axial_symmetry_bc(fl%qx, fl%fbcy_qx, dm%knc_sym, dm%dpcc)
      if( dm%is_thermo) &
      call get_cyindrical_x_axial_symmetry_bc(fl%gx, fl%fbcy_gx, dm%knc_sym, dm%dpcc)
    end if

    ! uy, dm0-dm1, np
    if(dm%ibcy(1, 2) /= IBC_DIRICHLET) call Print_error_msg('BC error.') ! check, axial-symmetric at y=2?
    fl%fbcy_qy (:, 1, :) = ZERO
    fl%fbcy_qyr(:, 1, :) = ZERO
    if( dm%is_thermo) then
      fl%fbcy_gy (:, 1, :) = ZERO
      fl%fbcy_gyr(:, 1, :) = ZERO
    end if

    ! uz, dm0-dm1, nc
    if(dm%ibcy(1, 3) == IBC_INTERIOR) then
      call get_cyindrical_x_axial_symmetry_bc(fl%qz, fl%fbcy_qz, dm%knc_sym, dm%dcpc)
      fl%fbcy_qzr(:, 1, :) = fl%fbcy_qz(:, 1, :) * dm%rci(1)
      fl%fbcy_qzr(:, 3, :) = fl%fbcy_qz(:, 3, :) * dm%rci(2)
      if( dm%is_thermo) then
        call get_cyindrical_x_axial_symmetry_bc(fl%gz, fl%fbcy_gz, dm%knc_sym, dm%dcpc)
        fl%fbcy_gzr(:, 1, :) = fl%fbcy_gz(:, 1, :) * dm%rci(1)
        fl%fbcy_gzr(:, 3, :) = fl%fbcy_gz(:, 3, :) * dm%rci(2)
      end if
    end if

    ! p, dm0-dm1, nc
    if(dm%ibcy(1, 4) == IBC_INTERIOR) then
      call get_cyindrical_x_axial_symmetry_bc(fl%pres, fl%fbcy_pr, dm%knc_sym, dm%dccc)
    end if


    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine apply_x_bc_2dm_halo(fbcx, var, iside, n)

    real(WP), intent(inout) :: fbcx(:, :, :)
    real(WP), intent(in)    :: var(:, :, :)
    integer, intent(in)     :: iside
    integer, intent(in)     :: n

    integer :: isign
    
    isign = 0
    if(iside == 1) then
      isign = 1
    else if(iside == 2) then
      isign = -1
    else
      call Print_error_msg('Error input.')
    end if

    fbcx(iside,   :, :) = var(n,       :, :)
    fbcx(iside+2, :, :) = var(n-isign, :, :)

    return
  end subroutine
!==========================================================================================================
  subroutine apply_y_bc_2dm_halo(fbcy, var, iside, n, dtmp)

    real(WP), intent(inout) :: fbcy(:, :, :)
    real(WP), intent(in)    :: var(:, :, :)
    integer, intent(in)     :: iside
    integer, intent(in)     :: n
    type(DECOMP_INFO), intent(in) :: dtmp


    real(WP), dimension( dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3) ) :: var_ypencil
    integer :: isign

    isign = 0
    if(iside == 1) then
      isign = 1
    else if(iside == 2) then
      isign = -1
    else
      call Print_error_msg('Error input.')
    end if

    call transpose_x_to_y(var, var_ypencil, dtmp)
    fbcy(:, iside,   :) = var_ypencil(:, n,       :)
    fbcy(:, iside+2, :) = var_ypencil(:, n-isign, :)

    return
  end subroutine
!==========================================================================================================
  subroutine apply_z_bc_2dm_halo(fbcz, var, iside, n, dtmp)

    real(WP), intent(inout) :: fbcz(:, :, :)
    real(WP), intent(in)    :: var(:, :, :)
    integer, intent(in)     :: iside
    integer, intent(in)     :: n
    type(DECOMP_INFO), intent(in) :: dtmp

    real(WP), dimension( dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3) ) :: var_ypencil
    real(WP), dimension( dtmp%zsz(1), dtmp%zsz(2), dtmp%zsz(3) ) :: var_zpencil
    integer :: isign

    isign = 0
    if(iside == 1) then
      isign = 1
    else if(iside == 2) then
      isign = -1
    else
      call Print_error_msg('Error input.')
    end if

    call transpose_x_to_y(var,         var_ypencil, dtmp)
    call transpose_y_to_z(var_ypencil, var_zpencil, dtmp)
    fbcz(:, :, iside  ) = var_zpencil(:, :, n      )
    fbcz(:, :, iside+2) = var_zpencil(:, :, n-isign)

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine update_flow_bc_2dm_halo(dm0, fl0, dm1, fl1, tm1)
    !use solver_tools_mod
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in)           :: dm0, dm1
    type(t_flow),   intent(inout)        :: fl0, fl1
    type(t_thermo), intent(in), optional :: tm1

    integer :: iside, n0
    real(WP), dimension( dm0%dcpc%xsz(1), dm0%dcpc%xsz(2), dm0%dcpc%xsz(3) ) :: acpc0_xpencil
    real(WP), dimension( dm1%dcpc%xsz(1), dm1%dcpc%xsz(2), dm1%dcpc%xsz(3) ) :: acpc1_xpencil
    real(WP), dimension( dm0%dccp%xsz(1), dm0%dccp%xsz(2), dm0%dccp%xsz(3) ) :: accp0_xpencil
    real(WP), dimension( dm1%dccp%xsz(1), dm1%dccp%xsz(2), dm1%dccp%xsz(3) ) :: accp1_xpencil
    
!----------------------------------------------------------------------------------------------------------
!   x-bc for qx, qy, qz, pres,  no overlap of values
!   iside = 1, 0-->1  : |-domain-0----3-1|-domain-1-----|
!   iside = 2, 0<--1  : |-domain-0---|-2-4-domain-1-----|
!----------------------------------------------------------------------------------------------------------
    iside = 1
    n0 = dm0%nc(1)
    if(dm1%ibcx(iside, 1) == IBC_INTERIOR) &
    call apply_x_bc_2dm_halo(fl1%fbcx_qx, fl0%qx,   iside, n0)

    if(dm1%ibcx(iside, 2) == IBC_INTERIOR) &
    call apply_x_bc_2dm_halo(fl1%fbcx_qy, fl0%qy,   iside, n0)

    if(dm1%ibcx(iside, 3) == IBC_INTERIOR) &
    call apply_x_bc_2dm_halo(fl1%fbcx_qz, fl0%qz,   iside, n0)

    if(dm1%ibcx(iside, 4) == IBC_INTERIOR) &
    call apply_x_bc_2dm_halo(fl1%fbcx_pr, fl0%pres, iside, n0)
    
    iside = 2
    n0 = 1
    if(dm0%ibcx(iside, 1) == IBC_INTERIOR) &
    call apply_x_bc_2dm_halo(fl0%fbcx_qx, fl1%qx,   iside, n0+1)

    if(dm0%ibcx(iside, 2) == IBC_INTERIOR) &
    call apply_x_bc_2dm_halo(fl0%fbcx_qy, fl1%qy,   iside, n0)

    if(dm0%ibcx(iside, 3) == IBC_INTERIOR) &
    call apply_x_bc_2dm_halo(fl0%fbcx_qz, fl1%qz,   iside, n0)

    if(dm0%ibcx(iside, 4) == IBC_INTERIOR) &
    call apply_x_bc_2dm_halo(fl0%fbcx_pr, fl1%pres, iside, n0)
!----------------------------------------------------------------------------------------------------------
!   y-bc for qx, qy, qz, pres,  no overlap of values
!   iside = 1, 0-->1  : |-domain-0----3-1|-domain-1-----|
!   iside = 2, 0<--1  : |-domain-0---|-2-4-domain-1-----|
!----------------------------------------------------------------------------------------------------------
    iside = 1
    n0 = dm0%nc(2)
    if(dm1%ibcy(iside, 1) == IBC_INTERIOR) &
    call apply_y_bc_2dm_halo(fl1%fbcy_qx, fl0%qx,   iside, n0, dm0%dpcc)

    if(dm1%ibcy(iside, 2) == IBC_INTERIOR) &
    call apply_y_bc_2dm_halo(fl1%fbcy_qy, fl0%qy,   iside, n0, dm0%dcpc)

    if(dm1%ibcy(iside, 3) == IBC_INTERIOR) &
    call apply_y_bc_2dm_halo(fl1%fbcy_qz, fl0%qz,   iside, n0, dm0%dccp)

    if(dm1%ibcy(iside, 4) == IBC_INTERIOR) &
    call apply_y_bc_2dm_halo(fl1%fbcy_pr, fl0%pres, iside, n0, dm0%dccc)

    iside = 2
    n0 = 1
    if(dm0%ibcy(iside, 1) == IBC_INTERIOR) &
    call apply_y_bc_2dm_halo(fl0%fbcy_qx, fl1%qx,   iside, n0,   dm1%dpcc)

    if(dm0%ibcy(iside, 2) == IBC_INTERIOR) &
    call apply_y_bc_2dm_halo(fl0%fbcy_qy, fl1%qy,   iside, n0+1, dm1%dcpc)

    if(dm0%ibcy(iside, 3) == IBC_INTERIOR) &
    call apply_y_bc_2dm_halo(fl0%fbcy_qz, fl1%qz,   iside, n0,   dm1%dccp)

    if(dm0%ibcy(iside, 4) == IBC_INTERIOR) &
    call apply_y_bc_2dm_halo(fl0%fbcy_pr, fl1%pres, iside, n0,   dm1%dccc)
!----------------------------------------------------------------------------------------------------------
!   z-bc for qx, qy, qz, pres,  no overlap of values
!   iside = 1, 0-->1  : |-domain-0----3-1|-domain-1-----|
!   iside = 2, 0<--1  : |-domain-0---|-2-4-domain-1-----|
!----------------------------------------------------------------------------------------------------------
    iside = 1
    n0 = dm0%nc(3)
    if(dm1%ibcz(iside, 1) == IBC_INTERIOR) &
    call apply_z_bc_2dm_halo(fl1%fbcz_qx, fl0%qx,   iside, n0, dm0%dpcc)

    if(dm1%ibcz(iside, 2) == IBC_INTERIOR) &
    call apply_z_bc_2dm_halo(fl1%fbcz_qy, fl0%qy,   iside, n0, dm0%dcpc)

    if(dm1%ibcz(iside, 3) == IBC_INTERIOR) &
    call apply_z_bc_2dm_halo(fl1%fbcz_qz, fl0%qz,   iside, n0, dm0%dccp)

    if(dm1%ibcz(iside, 4) == IBC_INTERIOR) &
    call apply_z_bc_2dm_halo(fl1%fbcz_pr, fl0%pres, iside, n0, dm0%dccc)
    
    iside = 2
    n0 = 1
    if(dm0%ibcz(iside, 1) == IBC_INTERIOR) &
    call apply_z_bc_2dm_halo(fl0%fbcz_qx, fl1%qx,   iside, n0,   dm1%dpcc)

    if(dm0%ibcz(iside, 2) == IBC_INTERIOR) &
    call apply_z_bc_2dm_halo(fl0%fbcz_qy, fl1%qy,   iside, n0,   dm1%dcpc)

    if(dm0%ibcz(iside, 3) == IBC_INTERIOR) &
    call apply_z_bc_2dm_halo(fl0%fbcz_qz, fl1%qz,   iside, n0+1, dm1%dccp)

    if(dm0%ibcz(iside, 4) == IBC_INTERIOR) &
    call apply_z_bc_2dm_halo(fl0%fbcz_pr, fl1%pres, iside, n0,   dm1%dccc)
!==========================================================================================================
!   thermal only
!==========================================================================================================
    if(dm1%is_thermo .and. dm0%is_thermo) then  
!----------------------------------------------------------------------------------------------------------
!   x-bc for gx, gy, gz, no overlap of values
!   iside = 1, 0-->1  : |-domain-0----3-1|-domain-1-----|
!   iside = 2, 0<--1  : |-domain-0---|-2-4-domain-1-----|
!----------------------------------------------------------------------------------------------------------
      iside = 1
      n0 = dm0%nc(1)
      if(dm1%ibcx(iside, 1) == IBC_INTERIOR) &
      call apply_x_bc_2dm_halo(fl1%fbcx_gx, fl0%gx,   iside, n0)

      if(dm1%ibcx(iside, 2) == IBC_INTERIOR) &
      call apply_x_bc_2dm_halo(fl1%fbcx_gy, fl0%gy,   iside, n0)

      if(dm1%ibcx(iside, 3) == IBC_INTERIOR) &
      call apply_x_bc_2dm_halo(fl1%fbcx_gz, fl0%gz,   iside, n0)
      
      iside = 2
      n0 = 1
      if(dm0%ibcx(iside, 1) == IBC_INTERIOR) &
      call apply_x_bc_2dm_halo(fl0%fbcx_gx, fl1%gx,   iside, n0+1)

      if(dm0%ibcx(iside, 2) == IBC_INTERIOR) &
      call apply_x_bc_2dm_halo(fl0%fbcx_gy, fl1%gy,   iside, n0)

      if(dm0%ibcx(iside, 3) == IBC_INTERIOR) &
      call apply_x_bc_2dm_halo(fl0%fbcx_gz, fl1%gz,   iside, n0)
!----------------------------------------------------------------------------------------------------------
!   y-bc for gx, gy, gz, no overlap of values
!   iside = 1, 0-->1  : |-domain-0----3-1|-domain-1-----|
!   iside = 2, 0<--1  : |-domain-0---|-2-4-domain-1-----|
!----------------------------------------------------------------------------------------------------------
      iside = 1
      n0 = dm0%nc(2)
      if(dm1%ibcy(iside, 1) == IBC_INTERIOR) &
      call apply_y_bc_2dm_halo(fl1%fbcy_gx, fl0%gx,   iside, n0, dm0%dpcc)

      if(dm1%ibcy(iside, 2) == IBC_INTERIOR) &
      call apply_y_bc_2dm_halo(fl1%fbcy_gy, fl0%gy,   iside, n0, dm0%dcpc)

      if(dm1%ibcy(iside, 3) == IBC_INTERIOR) &
      call apply_y_bc_2dm_halo(fl1%fbcy_gz, fl0%gz,   iside, n0, dm0%dccp)

      iside = 2
      n0 = 1
      if(dm0%ibcy(iside, 1) == IBC_INTERIOR) &
      call apply_y_bc_2dm_halo(fl0%fbcy_gx, fl1%gx,   iside, n0,   dm1%dpcc)

      if(dm0%ibcy(iside, 2) == IBC_INTERIOR) &
      call apply_y_bc_2dm_halo(fl0%fbcy_gy, fl1%gy,   iside, n0+1, dm1%dcpc)

      if(dm0%ibcy(iside, 3) == IBC_INTERIOR) &
      call apply_y_bc_2dm_halo(fl0%fbcy_gz, fl1%gz,   iside, n0,   dm1%dccp)
!----------------------------------------------------------------------------------------------------------
!   z-bc for gx, gy, gz,  no overlap of values
!   iside = 1, 0-->1  : |-domain-0----3-1|-domain-1-----|
!   iside = 2, 0<--1  : |-domain-0---|-2-4-domain-1-----|
!----------------------------------------------------------------------------------------------------------
      iside = 1
      n0 = dm0%nc(3)
      if(dm1%ibcz(iside, 1) == IBC_INTERIOR) &
      call apply_z_bc_2dm_halo(fl1%fbcz_gx, fl0%gx,   iside, n0, dm0%dpcc)

      if(dm1%ibcz(iside, 2) == IBC_INTERIOR) &
      call apply_z_bc_2dm_halo(fl1%fbcz_gy, fl0%gy,   iside, n0, dm0%dcpc)

      if(dm1%ibcz(iside, 3) == IBC_INTERIOR) &
      call apply_z_bc_2dm_halo(fl1%fbcz_gz, fl0%gz,   iside, n0, dm0%dccp)
      
      iside = 2
      n0 = 1
      if(dm0%ibcz(iside, 1) == IBC_INTERIOR) &
      call apply_z_bc_2dm_halo(fl0%fbcz_gx, fl1%gx,   iside, n0,   dm1%dpcc)

      if(dm0%ibcz(iside, 2) == IBC_INTERIOR) &
      call apply_z_bc_2dm_halo(fl0%fbcz_gy, fl1%gy,   iside, n0,   dm1%dcpc)

      if(dm0%ibcz(iside, 3) == IBC_INTERIOR) &
      call apply_z_bc_2dm_halo(fl0%fbcz_gz, fl1%gz,   iside, n0+1, dm1%dccp)

    else if (dm1%is_thermo .and. (.not.dm0%is_thermo)) then
      if(.not. present(tm1)) call Print_error_msg('Input Error.')
!----------------------------------------------------------------------------------------------------------
!   x, y, z-bc for gx, gy, gz, no overlap of values
!   iside = 1, 0-->1  : |-domain-0----3-1|-domain-1-----|
!----------------------------------------------------------------------------------------------------------
      iside = 1
      do n0 = 1, 3, 2
        if(dm1%ibcx(iside, 1) == IBC_INTERIOR) &
        fl1%fbcx_gx(n0, :, :) = fl0%fbcx_qx(n0, :, :) * tm1%fbcx_ftp(1, :, :)%d

        if(dm1%ibcx(iside, 2) == IBC_INTERIOR) &
        fl1%fbcx_gy(n0, :, :) = fl0%fbcx_qy(n0, :, :) * tm1%fbcx_ftp(1, :, :)%d

        if(dm1%ibcx(iside, 3) == IBC_INTERIOR) &
        fl1%fbcx_gz(n0, :, :) = fl0%fbcx_qz(n0, :, :) * tm1%fbcx_ftp(1, :, :)%d

        if(dm1%ibcy(iside, 1) == IBC_INTERIOR) &
        fl1%fbcy_gx(:, n0, :) = fl1%fbcy_qx(:, n0, :) * tm1%fbcy_ftp(:, 1, :)%d
        
        if(dm1%ibcy(iside, 2) == IBC_INTERIOR) &
        fl1%fbcy_gy(:, n0, :) = fl1%fbcy_qy(:, n0, :) * tm1%fbcy_ftp(:, 1, :)%d
        
        if(dm1%ibcy(iside, 3) == IBC_INTERIOR) &
        fl1%fbcy_gz(:, n0, :) = fl1%fbcy_qz(:, n0, :) * tm1%fbcy_ftp(:, 1, :)%d

        if(dm1%ibcz(iside, 1) == IBC_INTERIOR) &
        fl1%fbcz_gx(:, :, n0) = fl1%fbcz_qx(:, :, n0) * tm1%fbcz_ftp(:, :, 1)%d
        
        if(dm1%ibcz(iside, 2) == IBC_INTERIOR) &
        fl1%fbcz_gy(:, :, n0) = fl1%fbcz_qy(:, :, n0) * tm1%fbcz_ftp(:, :, 1)%d

        if(dm1%ibcz(iside, 3) == IBC_INTERIOR) &
        fl1%fbcz_gz(:, :, n0) = fl1%fbcz_qz(:, :, n0) * tm1%fbcz_ftp(:, :, 1)%d
      end do
    else 
      call Print_error_msg("Error in setting up thermo b.c.")
    end if
!==========================================================================================================
!   clyindrical only
!==========================================================================================================
    if(dm1%icoordinate == ICYLINDRICAL .or. dm0%icoordinate == ICYLINDRICAL) then
!----------------------------------------------------------------------------------------------------------
!   y-bc for qy/r, qz/r,  no overlap of values
!   iside = 1, 0-->1  : |-domain-0----3-1|-domain-1-----|
!   iside = 2, 0<--1  : |-domain-0---|-2-4-domain-1-----|
!----------------------------------------------------------------------------------------------------------
      iside = 1
      n0 = dm0%nc(2)
      if(dm1%ibcy(iside, 2) == IBC_INTERIOR) then
        acpc0_xpencil = fl0%qy
        call multiple_cylindrical_rn(acpc0_xpencil, dm0%dcpc, dm0%rpi, 1, IPENCIL(1))
        call apply_y_bc_2dm_halo(fl1%fbcy_qyr, acpc0_xpencil, iside, n0, dm0%dcpc)
      end if

      if(dm1%ibcy(iside, 3) == IBC_INTERIOR) then
        accp0_xpencil = fl0%qz
        call multiple_cylindrical_rn(accp0_xpencil, dm0%dccp, dm0%rci, 1, IPENCIL(1))
        call apply_y_bc_2dm_halo(fl1%fbcy_qzr, accp0_xpencil, iside, n0, dm0%dccp)
      end if

      iside = 2
      n0 = 1
      if(dm0%ibcy(iside, 2) == IBC_INTERIOR) then
        acpc1_xpencil = fl1%qy
        call multiple_cylindrical_rn(acpc1_xpencil, dm1%dcpc, dm1%rpi, 1, IPENCIL(1))
        call apply_y_bc_2dm_halo(fl0%fbcy_qyr, acpc1_xpencil, iside, n0+1, dm1%dcpc)
      end if

      if(dm0%ibcy(iside, 3) == IBC_INTERIOR) then
        accp1_xpencil = fl1%qz
        call multiple_cylindrical_rn(accp1_xpencil, dm1%dccp, dm1%rci, 1, IPENCIL(1))
        call apply_y_bc_2dm_halo(fl0%fbcy_qzr, accp1_xpencil, iside, n0,   dm1%dccp)
      end if
!----------------------------------------------------------------------------------------------------------
!   z-bc for qy/r, qz/r,  no overlap of values
!   iside = 1, 0-->1  : |-domain-0----3-1|-domain-1-----|
!   iside = 2, 0<--1  : |-domain-0---|-2-4-domain-1-----|
!----------------------------------------------------------------------------------------------------------
      iside = 1
      n0 = dm0%nc(3)
      if(dm1%ibcz(iside, 2) == IBC_INTERIOR) then
        acpc0_xpencil = fl0%qy
        call multiple_cylindrical_rn(acpc0_xpencil, dm0%dcpc, dm0%rpi, 1, IPENCIL(1))
        call apply_z_bc_2dm_halo(fl1%fbcz_qyr, acpc0_xpencil, iside, n0, dm0%dcpc)
      end if

      if(dm1%ibcy(iside, 3) == IBC_INTERIOR) then
        accp0_xpencil = fl0%qz
        call multiple_cylindrical_rn(accp0_xpencil, dm0%dccp, dm0%rci, 1, IPENCIL(1))
        call apply_z_bc_2dm_halo(fl1%fbcz_qzr, accp0_xpencil, iside, n0, dm0%dccp)
      end if

      iside = 2
      n0 = 1
      if(dm0%ibcy(iside, 2) == IBC_INTERIOR) then
        acpc1_xpencil = fl1%qy
        call multiple_cylindrical_rn(acpc1_xpencil, dm1%dcpc, dm1%rpi, 1, IPENCIL(1))
        call apply_z_bc_2dm_halo(fl0%fbcz_qyr, acpc1_xpencil, iside, n0, dm1%dcpc)
      end if

      if(dm0%ibcy(iside, 3) == IBC_INTERIOR) then
        accp1_xpencil = fl1%qz
        call multiple_cylindrical_rn(accp1_xpencil, dm1%dccp, dm1%rci, 1, IPENCIL(1))
        call apply_z_bc_2dm_halo(fl0%fbcz_qzr, accp1_xpencil, iside, n0+1, dm1%dccp)
      end if

      if(dm1%is_thermo .and. dm0%is_thermo) then  
!----------------------------------------------------------------------------------------------------------
!   y-bc for gy/r, gz/r,  no overlap of values
!   iside = 1, 0-->1  : |-domain-0----3-1|-domain-1-----|
!   iside = 2, 0<--1  : |-domain-0---|-2-4-domain-1-----|
!----------------------------------------------------------------------------------------------------------
        iside = 1
        n0 = dm0%nc(2)
        if(dm1%ibcy(iside, 2) == IBC_INTERIOR) then
          acpc0_xpencil = fl0%gy
          call multiple_cylindrical_rn(acpc0_xpencil, dm0%dcpc, dm0%rpi, 1, IPENCIL(1))
          call apply_y_bc_2dm_halo(fl1%fbcy_gyr, acpc0_xpencil, iside, n0, dm0%dcpc)
        end if

        if(dm1%ibcy(iside, 3) == IBC_INTERIOR) then
          accp0_xpencil = fl0%gz
          call multiple_cylindrical_rn(accp0_xpencil, dm0%dccp, dm0%rci, 1, IPENCIL(1))
          call apply_y_bc_2dm_halo(fl1%fbcy_gzr, accp0_xpencil, iside, n0, dm0%dccp)
        end if

        iside = 2
        n0 = 1
        if(dm0%ibcy(iside, 2) == IBC_INTERIOR) then
          acpc1_xpencil = fl1%gy
          call multiple_cylindrical_rn(acpc1_xpencil, dm1%dcpc, dm1%rpi, 1, IPENCIL(1))
          call apply_y_bc_2dm_halo(fl0%fbcy_gyr, acpc1_xpencil, iside, n0+1, dm1%dcpc)
        end if

        if(dm0%ibcy(iside, 3) == IBC_INTERIOR) then
          accp1_xpencil = fl1%gz
          call multiple_cylindrical_rn(accp1_xpencil, dm1%dccp, dm1%rci, 1, IPENCIL(1))
          call apply_y_bc_2dm_halo(fl0%fbcy_gzr, accp1_xpencil, iside, n0,   dm1%dccp)
        end if
!----------------------------------------------------------------------------------------------------------
!   z-bc for qy/r, qz/r,  no overlap of values
!   iside = 1, 0-->1  : |-domain-0----3-1|-domain-1-----|
!   iside = 2, 0<--1  : |-domain-0---|-2-4-domain-1-----|
!----------------------------------------------------------------------------------------------------------
        iside = 1
        n0 = dm0%nc(3)
        if(dm1%ibcz(iside, 2) == IBC_INTERIOR) then
          acpc0_xpencil = fl0%gy
          call multiple_cylindrical_rn(acpc0_xpencil, dm0%dcpc, dm0%rpi, 1, IPENCIL(1))
          call apply_z_bc_2dm_halo(fl1%fbcz_gyr, acpc0_xpencil, iside, n0, dm0%dcpc)
        end if

        if(dm1%ibcy(iside, 3) == IBC_INTERIOR) then
          accp0_xpencil = fl0%gz
          call multiple_cylindrical_rn(accp0_xpencil, dm0%dccp, dm0%rci, 1, IPENCIL(1))
          call apply_z_bc_2dm_halo(fl1%fbcz_gzr, accp0_xpencil, iside, n0, dm0%dccp)
        end if

        iside = 2
        n0 = 1
        if(dm0%ibcy(iside, 2) == IBC_INTERIOR) then
          acpc1_xpencil = fl1%gy
          call multiple_cylindrical_rn(acpc1_xpencil, dm1%dcpc, dm1%rpi, 1, IPENCIL(1))
          call apply_z_bc_2dm_halo(fl0%fbcz_gyr, acpc1_xpencil, iside, n0, dm1%dcpc)
        end if

        if(dm0%ibcy(iside, 3) == IBC_INTERIOR) then
          accp1_xpencil = fl1%gz
          call multiple_cylindrical_rn(accp1_xpencil, dm1%dccp, dm1%rci, 1, IPENCIL(1))
          call apply_z_bc_2dm_halo(fl0%fbcz_gzr, accp1_xpencil, iside, n0+1, dm1%dccp)
        end if

      else if (dm1%is_thermo .and. (.not.dm0%is_thermo)) then
        if(.not. present(tm1)) call Print_error_msg('Input Error.')
!----------------------------------------------------------------------------------------------------------
!    y, z-bc for gyr, gzr, no overlap of values
!   iside = 1, 0-->1  : |-domain-0----3-1|-domain-1-----|
!----------------------------------------------------------------------------------------------------------
        iside = 1
        do n0 = 1, 3, 2
          
          !if(dm1%ibcx(iside, 2) == IBC_INTERIOR) &
          !fl1%fbcx_gyr(n0, :, :) = fl0%fbcx_qyr(n0, :, :) * tm1%fbcx_ftp(1, :, :)%d

          !if(dm1%ibcx(iside, 3) == IBC_INTERIOR) &
          !fl1%fbcx_gzr(n0, :, :) = fl0%fbcx_qzr(n0, :, :) * tm1%fbcx_ftp(1, :, :)%d

          
          if(dm1%ibcy(iside, 2) == IBC_INTERIOR) &
          fl1%fbcy_gyr(:, n0, :) = fl1%fbcy_qyr(:, n0, :) * tm1%fbcy_ftp(:, 1, :)%d
          
          if(dm1%ibcy(iside, 3) == IBC_INTERIOR) &
          fl1%fbcy_gzr(:, n0, :) = fl1%fbcy_qzr(:, n0, :) * tm1%fbcy_ftp(:, 1, :)%d

          
          if(dm1%ibcz(iside, 2) == IBC_INTERIOR) &
          fl1%fbcz_gyr(:, :, n0) = fl1%fbcz_qyr(:, :, n0) * tm1%fbcz_ftp(:, :, 1)%d

          if(dm1%ibcz(iside, 3) == IBC_INTERIOR) &
          fl1%fbcz_gzr(:, :, n0) = fl1%fbcz_qzr(:, :, n0) * tm1%fbcz_ftp(:, :, 1)%d
        end do
      else 
        call Print_error_msg("Error in setting up thermo b.c.")
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
    
    integer :: i, j, k
    real(WP), dimension( dm0%dccc%ysz(1), dm0%dccc%ysz(2), dm0%dccc%ysz(3) ) :: accc0_ypencil
    real(WP), dimension( dm0%dccc%zsz(1), dm0%dccc%zsz(2), dm0%dccc%zsz(3) ) :: accc0_zpencil
    real(WP), dimension( dm1%dccc%ysz(1), dm1%dccc%ysz(2), dm1%dccc%ysz(3) ) :: accc1_ypencil
    real(WP), dimension( dm1%dccc%zsz(1), dm1%dccc%zsz(2), dm1%dccc%zsz(3) ) :: accc1_zpencil
!----------------------------------------------------------------------------------------------------------
!   bc in x - direction
!----------------------------------------------------------------------------------------------------------
    ! thermal field, dm0-dm1
    if(dm1%ibcx(1, 5) == IBC_INTERIOR) then
      tm1%fbcx_ftp(1, :, :)%t = tm0%tTemp(dm0%nc(1),     :, :)
      tm1%fbcx_ftp(3, :, :)%t = tm0%tTemp(dm0%nc(1) - 1, :, :)
      do k = 1, dm1%dccc%xsz(3)
        do j = 1, dm1%dccc%xsz(2)
          call ftp_refresh_thermal_properties_from_T_undim(tm1%fbcx_ftp(1, j, k))
          call ftp_refresh_thermal_properties_from_T_undim(tm1%fbcx_ftp(3, j, k))
        end do
      end do
    end if
    if(dm0%ibcx(2, 5) == IBC_INTERIOR) then
      tm0%fbcx_ftp(2, :, :)%t = tm1%tTemp(1, :, :)
      tm0%fbcx_ftp(4, :, :)%t = tm1%tTemp(2, :, :)
      do k = 1, dm0%dccc%xsz(3)
        do j = 1, dm0%dccc%xsz(2)
          call ftp_refresh_thermal_properties_from_T_undim(tm0%fbcx_ftp(2, j, k))
          call ftp_refresh_thermal_properties_from_T_undim(tm0%fbcx_ftp(4, j, k))
        end do
      end do
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in y - direction
!----------------------------------------------------------------------------------------------------------
    ! thermal field, dm0-dm1
    if(dm1%ibcy(1, 5) == IBC_INTERIOR) then
      call transpose_x_to_y(tm0%tTemp, accc0_ypencil, dm0%dccc)
      tm1%fbcy_ftp(:, 1, :)%t = accc0_ypencil(:, dm0%nc(1),     :)
      tm1%fbcy_ftp(:, 3, :)%t = accc0_ypencil(:, dm0%nc(1) - 1, :)
      do k = 1, dm1%dccc%ysz(3)
        do i = 1, dm1%dccc%ysz(1)
          call ftp_refresh_thermal_properties_from_T_undim(tm1%fbcy_ftp(i, 1, k))
          call ftp_refresh_thermal_properties_from_T_undim(tm1%fbcy_ftp(i, 3, k))
        end do
      end do
    end if

    if(dm0%ibcy(2, 5) == IBC_INTERIOR) then
      call transpose_x_to_y(tm1%tTemp, accc1_ypencil, dm1%dccc)
      tm0%fbcy_ftp(:, 2, :)%t = accc1_ypencil(:, 1, :)
      tm0%fbcy_ftp(:, 4, :)%t = accc1_ypencil(:, 2, :)
      do k = 1, dm0%dccc%ysz(3)
        do i = 1, dm0%dccc%ysz(1)
          call ftp_refresh_thermal_properties_from_T_undim(tm0%fbcy_ftp(i, 2, k))
          call ftp_refresh_thermal_properties_from_T_undim(tm0%fbcy_ftp(i, 4, k))
        end do
      end do
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in z - direction
!----------------------------------------------------------------------------------------------------------
    ! thermal field, dm0-dm1
    if(dm1%ibcz(1, 5) == IBC_INTERIOR) then
      call transpose_x_to_y(tm0%tTemp,     accc0_ypencil, dm0%dccc)
      call transpose_y_to_z(accc0_ypencil, accc0_zpencil, dm0%dccc)
      tm1%fbcz_ftp(:, :, 1)%t = accc0_zpencil(:, :, dm0%nc(1)    )
      tm1%fbcz_ftp(:, :, 3)%t = accc0_zpencil(:, :, dm0%nc(1) - 1)
      do j = 1, dm1%dccc%zsz(2)
        do i = 1, dm1%dccc%zsz(1)
          call ftp_refresh_thermal_properties_from_T_undim(tm1%fbcz_ftp(i, j, 1))
          call ftp_refresh_thermal_properties_from_T_undim(tm1%fbcz_ftp(i, j, 3))
        end do
      end do
    end if

    if(dm0%ibcz(2, 5) == IBC_INTERIOR) then
      call transpose_x_to_y(tm1%tTemp,     accc1_ypencil, dm1%dccc)
      call transpose_y_to_z(accc1_ypencil, accc1_zpencil, dm1%dccc)
      tm0%fbcz_ftp(:, :, 2)%t = accc1_zpencil(:, :, 1)
      tm0%fbcz_ftp(:, :, 4)%t = accc1_zpencil(:, :, 2)
      do j = 1, dm0%dccc%zsz(2)
        do i = 1, dm0%dccc%zsz(1)
          call ftp_refresh_thermal_properties_from_T_undim(tm0%fbcz_ftp(i, j, 2))
          call ftp_refresh_thermal_properties_from_T_undim(tm0%fbcz_ftp(i, j, 4))
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
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc1_zpencil

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
      call transpose_y_to_z(accc_ypencil, accc_zpencil, dm%dccc)
      do k = 1, dm%dccc%zsz(3)
        accc1_zpencil(:, :, k) = accc_zpencil(:, :, dm%knc_sym(k))
      end do
      call transpose_z_to_y(accc1_zpencil, accc_ypencil, dm%dpcc)

      tm%fbcy_ftp(:, 1, :)%t = accc_ypencil(:, 1, :)
      tm%fbcy_ftp(:, 3, :)%t = accc_ypencil(:, 2, :)
      do k = 1, dm%dccc%ysz(3)
        do i = 1, dm%dccc%ysz(1)
          call ftp_refresh_thermal_properties_from_T_undim(tm%fbcy_ftp(i, 1, k))
          call ftp_refresh_thermal_properties_from_T_undim(tm%fbcy_ftp(i, 3, k))
        end do
      end do
    end if

    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine apply_gxgygz_bc_geo (dm, fl, tm) ! 
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_flow), intent( inout) :: fl
    type(t_thermo), intent(in)   :: tm
    type(t_domain), intent(in)   :: dm

    if(.not. dm%is_thermo) return

!----------------------------------------------------------------------------------------------------------
!   get bc gx, gy, gz (at geo-bc, not cell centre)
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


!==========================================================================================================
!==========================================================================================================
  subroutine get_dirichlet_geo_bcx(ibcx, varx, fbcx)
    use parameters_constant_mod
    implicit none
    integer,  intent(in)    :: ibcx(2)
    real(WP), intent(in)    :: varx(:, :, :)
    real(WP), intent(inout) :: fbcx(:, :, :)

    integer :: n

    if( any(ibcx /= IBC_DIRICHLET) ) return

    if(size(varx, 3) /= size(fbcx, 3) .or. &
       size(varx, 2) /= size(fbcx, 2)) call Print_error_msg("Error input.")

    if (ibcx(1) == IBC_DIRICHLET) then
      fbcx(1, :, :) = varx(1, :, :)
      fbcx(3, :, :) = fbcx(1, :, :)
    end if

    if (ibcx(2) == IBC_DIRICHLET) then
      n = size(varx, 1)
      fbcx(2, :, :) = varx(n, :, :)
      fbcx(4, :, :) = fbcx(2, :, :)
    end if

    return
  end subroutine
!==========================================================================================================
  subroutine get_dirichlet_geo_bcy(ibcy, vary, fbcy)
    use parameters_constant_mod
    implicit none
    integer,  intent(in)    :: ibcy(2)
    real(WP), intent(in)    :: vary(:, :, :)
    real(WP), intent(inout) :: fbcy(:, :, :)

    integer :: n

    if( any(ibcy /= IBC_DIRICHLET) ) return

    if(size(vary, 1) /= size(fbcy, 1) .or. &
       size(vary, 3) /= size(fbcy, 3)) call Print_error_msg("Error input.")

    if (ibcy(1) == IBC_DIRICHLET) then
      fbcy(:, 1, :) = vary(:, 1, :)
      fbcy(:, 3, :) = fbcy(:, 1, :)
    end if

    if (ibcy(2) == IBC_DIRICHLET) then
      n = size(vary, 2)
      fbcy(:, 2, :) = vary(:, n, :)
      fbcy(:, 4, :) = fbcy(:, 2, :)
    end if

    return
  end subroutine
!==========================================================================================================
  subroutine get_dirichlet_geo_bcz(ibcz, varz, fbcz)
    use parameters_constant_mod
    implicit none
    integer,  intent(in)    :: ibcz(2)
    real(WP), intent(in)    :: varz(:, :, :)
    real(WP), intent(inout) :: fbcz(:, :, :)

    integer :: n

    if( any(ibcz /= IBC_DIRICHLET) ) return

    if ( ANY(ibcz /= IBC_DIRICHLET) ) return

    if(size(varz, 1) /= size(fbcz, 1) .or. &
       size(varz, 2) /= size(fbcz, 2)) call Print_error_msg("Error input.")

    if (ibcz(1) == IBC_DIRICHLET) then
      fbcz(:, :, 1) = varz(:, :, 1)
      fbcz(:, :, 3) = fbcz(:, :, 1)
    end if

    if (ibcz(2) == IBC_DIRICHLET) then
      n = size(varz, 3)
      fbcz(:, :, 2) = varz(:, :, n)
      fbcz(:, :, 4) = fbcz(:, :, 2)
    end if

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
end module
