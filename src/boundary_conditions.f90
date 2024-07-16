module boundary_conditions_mod
  use udf_type_mod
  use parameters_constant_mod
  use print_msg_mod
  implicit none
  character(18) :: filename(5)
  integer, parameter :: IFBC(1:2) = (/1, 2/)
  integer, save :: mbcx_cov1(2), &
                   mbcy_cov1(2), &
                   mbcz_cov1(2), &
                   mbcx_tau1(2), &
                   mbcy_tau1(2), &
                   mbcz_tau1(2), &
                   mbcx_cov2(2), &
                   mbcy_cov2(2), &
                   mbcz_cov2(2), &
                   mbcr_cov2(2), &
                   mbcy_tau2(2), &
                   mbcx_tau2(2), &
                   mbcz_tau2(2), &
                   mbcr_tau2(2), &
                   mbcx_cov3(2), &
                   mbcy_cov3(2), &
                   mbcz_cov3(2), &
                   mbcr_cov3(2), &
                   mbcy_tau3(2), &
                   mbcx_tau3(2), &
                   mbcz_tau3(2), &
                   mbcr_tau3(2), &
                   ebcx_conv(2), &
                   ebcy_conv(2), &
                   ebcz_conv(2), &
                   ebcx_difu(2), &
                   ebcy_difu(2), &
                   ebcz_difu(2)
  logical, save :: is_fbcx_velo_required, &
                   is_fbcy_velo_required, &
                   is_fbcz_velo_required

  private :: reassign_ibc
  public  :: config_calc_basic_ibc  ! applied once only, just before calculation

  public  :: allocate_fbc_flow   ! applied once only
  public  :: allocate_fbc_thermo ! applied once only

  private :: map_bc_1d_uprofile     
  private :: apply_fbcx_given_profile 
  private :: apply_fbcx_given_const
  private :: apply_fbcy_given_const
  private :: apply_fbcz_given_const
  public  :: apply_fbc_given_flow   ! applied once only, for bc of constant velocity
  public  :: apply_fbc_given_thermo ! applied once only, for bc of constant temperature
  
  

  private :: get_fbcy_circle_centre
  public  :: update_fbcy_cc_flow_halo   ! for pipe only, applied every NS, cc for circle central point and var stored in xcx
  public  :: update_fbcy_cc_thermo_halo ! for pipe only, applied every NS, cc for circle central point and var stored in xcx

  private :: apply_fbcx_2dm_halo
  private :: apply_fbcy_2dm_halo
  private :: apply_fbcz_2dm_halo
  private :: apply_fbc_2dm_flow_halo
  public  :: update_fbc_2dm_flow_halo   ! for multiple domains only, update every NS 
  public  :: update_fbc_2dm_thermo_halo ! for multiple domains only, update every NS
  
  public :: reconstruct_symmetry_ibc    ! applied if necessary
  public  :: config_calc_eqs_ibc

contains
!==========================================================================================================
!==========================================================================================================
! function: re-assign calcuation ibc and keep the nominal bc
  subroutine reassign_ibc(bc_nominal, ibc)
    integer, intent(in) :: bc_nominal(2, 5)
    integer, intent(out) :: ibc(2, 5)
    integer :: n, m

    do n = 1, 2
      do m = 1, 5
        if (bc_nominal(n, m) == IBC_PROFILE1D)   then
          ibc(n, m) = IBC_DIRICHLET
        else if (bc_nominal(n, m) == IBC_TURBGEN  .or. &
                 bc_nominal(n, m) == IBC_DATABASE )   then
          if(m == 5) then
            ibc(n, m) = IBC_DIRICHLET ! for temperature, default is no incoming thermal flow, it is initilazed temperature
          else 
            ibc(n, m) = IBC_INTERIOR  ! for u, v, w, p
          end if
        else if (bc_nominal(n, m) == IBC_CONVECTIVE)   then ! check for convetive outlet
          ibc(n, m) = IBC_INTRPL
        else
          ibc(n, m) = bc_nominal(n, m)   
        end if
      end do
    end do

    return
  end subroutine 
!==========================================================================================================
!==========================================================================================================
! to get all ibc for calculation
!==========================================================================================================
  subroutine config_calc_basic_ibc(dm)
    use wtformat_mod
    type(t_domain), intent(inout) :: dm
    integer :: n
    integer :: ibcx(2, 5), ibcy(2, 5), ibcz(2, 5)
!----------------------------------------------------------------------------------------------------------
! to check velocity symmetric and asymmetric
!----------------------------------------------------------------------------------------------------------
    do n = 1, 2
      if(dm%ibcx_nominal(n, 1) == IBC_SYMMETRIC) &
         dm%ibcx_nominal(n, 1) =  IBC_ASYMMETRIC
      if(dm%ibcy_nominal(n, 2) == IBC_SYMMETRIC) &
         dm%ibcy_nominal(n, 2) =  IBC_ASYMMETRIC
      if(dm%ibcz_nominal(n, 3) == IBC_SYMMETRIC) &
         dm%ibcz_nominal(n, 3) =  IBC_ASYMMETRIC
    end do
!----------------------------------------------------------------------------------------------------------
! to set up real bc for calculation from given nominal b.c.
!----------------------------------------------------------------------------------------------------------
    call reassign_ibc(dm%ibcx_nominal, ibcx(1:2, 1:5))
    call reassign_ibc(dm%ibcy_nominal, ibcy(1:2, 1:5))
    call reassign_ibc(dm%ibcz_nominal, ibcz(1:2, 1:5))
!----------------------------------------------------------------------------------------------------------
! allocate bc to variables
!----------------------------------------------------------------------------------------------------------
    do n = 1, 2
      dm%ibcx_qx(n) = ibcx(n, 1)
      dm%ibcx_qy(n) = ibcx(n, 2)
      dm%ibcx_qz(n) = ibcx(n, 3)
      dm%ibcx_pr(n) = ibcx(n, 4)
      dm%ibcx_Th(n) = ibcx(n, 5)

      dm%ibcy_qx(n) = ibcy(n, 1)
      dm%ibcy_qy(n) = ibcy(n, 2)
      dm%ibcy_qz(n) = ibcy(n, 3)
      dm%ibcy_pr(n) = ibcy(n, 4)
      dm%ibcy_Th(n) = ibcy(n, 5)

      dm%ibcz_qx(n) = ibcz(n, 1)
      dm%ibcz_qy(n) = ibcz(n, 2)
      dm%ibcz_qz(n) = ibcz(n, 3)
      dm%ibcz_pr(n) = ibcz(n, 4)
      dm%ibcz_Th(n) = ibcz(n, 5)
    end do 

    if(nrank == 0) then

      write (*, wrtfmt1s) '  Boundary type options : '
      write (*, wrtfmt1s) '   0  = IBC_INTERIOR'
      write (*, wrtfmt1s) '   1  = IBC_PERIODIC'
      write (*, wrtfmt1s) '   2  = IBC_SYMMETRIC'
      write (*, wrtfmt1s) '   3  = IBC_ASYMMETRIC'
      write (*, wrtfmt1s) '   4  = IBC_DIRICHLET'
      write (*, wrtfmt1s) '   5  = IBC_NEUMANN'
      write (*, wrtfmt1s) '   6  = IBC_INTRPL'
      write (*, wrtfmt1s) '   7  = IBC_CONVECTIVE'
      write (*, wrtfmt1s) '   8  = IBC_TURBGEN'
      write (*, wrtfmt1s) '   9  = IBC_PROFILE1D'
      write (*, wrtfmt1s) '   10 = IBC_DATABASE'

      write (*, *) 'is periodic in xyz? ', dm%is_periodic(1:3)
      write (*, wrtfmt1s) 'BC in the X direction: norminal BC, calc BC'
      write (*, wrtfmt4i2r) '  u-bc :', dm%ibcx_nominal(1:2, 1), dm%ibcx_qx(1:2), dm%fbcx_const(1:2, 1)
      write (*, wrtfmt4i2r) '  v-bc :', dm%ibcx_nominal(1:2, 2), dm%ibcx_qy(1:2), dm%fbcx_const(1:2, 2)
      write (*, wrtfmt4i2r) '  w-bc :', dm%ibcx_nominal(1:2, 3), dm%ibcx_qz(1:2), dm%fbcx_const(1:2, 3)
      write (*, wrtfmt4i2r) '  p-bc :', dm%ibcx_nominal(1:2, 4), dm%ibcx_pr(1:2), dm%fbcx_const(1:2, 4)
      write (*, wrtfmt4i2r) '  T-bc :', dm%ibcx_nominal(1:2, 5), dm%ibcx_Th(1:2), dm%fbcx_const(1:2, 5)
      write (*, wrtfmt1s) 'BC in the Y direction: norminal BC, calc BC'
      write (*, wrtfmt4i2r) '  u-bc :', dm%ibcy_nominal(1:2, 1), dm%ibcy_qx(1:2), dm%fbcy_const(1:2, 1)
      write (*, wrtfmt4i2r) '  v-bc :', dm%ibcy_nominal(1:2, 2), dm%ibcy_qy(1:2), dm%fbcy_const(1:2, 2)
      write (*, wrtfmt4i2r) '  w-bc :', dm%ibcy_nominal(1:2, 3), dm%ibcy_qz(1:2), dm%fbcy_const(1:2, 3)
      write (*, wrtfmt4i2r) '  p-bc :', dm%ibcy_nominal(1:2, 4), dm%ibcy_pr(1:2), dm%fbcy_const(1:2, 4)
      write (*, wrtfmt4i2r) '  T-bc :', dm%ibcy_nominal(1:2, 5), dm%ibcy_Th(1:2), dm%fbcy_const(1:2, 5)
      write (*, wrtfmt1s) 'BC in the Z direction: norminal BC, calc BC'
      write (*, wrtfmt4i2r) '  u-bc :', dm%ibcz_nominal(1:2, 1), dm%ibcz_qx(1:2), dm%fbcz_const(1:2, 1)
      write (*, wrtfmt4i2r) '  v-bc :', dm%ibcz_nominal(1:2, 2), dm%ibcz_qy(1:2), dm%fbcz_const(1:2, 2)
      write (*, wrtfmt4i2r) '  w-bc :', dm%ibcz_nominal(1:2, 3), dm%ibcz_qz(1:2), dm%fbcz_const(1:2, 3)
      write (*, wrtfmt4i2r) '  p-bc :', dm%ibcz_nominal(1:2, 4), dm%ibcz_pr(1:2), dm%fbcz_const(1:2, 4)
      write (*, wrtfmt4i2r) '  T-bc :', dm%ibcz_nominal(1:2, 5), dm%ibcz_Th(1:2), dm%fbcz_const(1:2, 5)
    end if

    return
  end subroutine 
!==========================================================================================================
!==========================================================================================================
  subroutine  map_bc_1d_uprofile(filename, n, y, u)
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

    call profile_interpolation(nn, yy, uprofile, n, y, u)

    deallocate(uprofile)
    deallocate(yy)

    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine apply_fbcx_given_profile(fbcx, var1y, jst)
    real(WP), intent(inout) :: fbcx(:, :, :)
    real(WP), intent(in)    :: var1y(:)
    integer,  intent(in)    :: jst

    integer :: k, j, jj

    do k = 1, size(fbcx, 3)
      do j = 1, size(fbcx, 2)
        jj = jst + j - 1
        fbcx(1, j, k) = var1y(jj)
        fbcx(3, j, k) = fbcx(1, j, k)
      end do
    end do

    return
  end subroutine
!==========================================================================================================
  subroutine apply_fbcx_given_const(fbcx, fbcx_const)
    real(WP), intent(inout) :: fbcx(:, :, :)
    real(WP), intent(in)    :: fbcx_const(4)

    integer :: k, j, n
    
    do k = 1, size(fbcx, 3)
      do j = 1, size(fbcx, 2)
        do n = 1, 4
          fbcx(n, j, k) =  fbcx_const(n)
        end do
      end do
    end do

    return
  end subroutine
!==========================================================================================================
  subroutine apply_fbcy_given_const(fbcy, fbcy_const, ri)
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
  subroutine apply_fbcz_given_const(fbcz, fbcz_const, ri, jst)
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
!==========================================================================================================
  subroutine allocate_fbc_flow(dm)
    type(t_domain), intent(inout)  :: dm
!----------------------------------------------------------------------------------------------------------
! to set up real bc values for calculation from given nominal b.c. values
! bc always saved on the boundar face centre 
! warning: this bc treatment is not proper for a inlet plane with field data.... to check and to update
!----------------------------------------------------------------------------------------------------------
    allocate( dm%fbcx_qx(             4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)) )! default x pencil
    allocate( dm%fbcy_qx(dm%dpcc%ysz(1),              4, dm%dpcc%ysz(3)) )! default y pencil
    allocate( dm%fbcz_qx(dm%dpcc%zsz(1), dm%dpcc%zsz(2),              4) )! default z pencil

    allocate( dm%fbcx_qy(             4, dm%dcpc%xsz(2), dm%dcpc%xsz(3)) )! default x pencil
    allocate( dm%fbcy_qy(dm%dcpc%ysz(1),              4, dm%dcpc%ysz(3)) )! default y pencil
    allocate( dm%fbcz_qy(dm%dcpc%zsz(1), dm%dcpc%zsz(2),              4) )! default z pencil

    allocate( dm%fbcx_qz(             4, dm%dccp%xsz(2), dm%dccp%xsz(3)) )! default x pencil
    allocate( dm%fbcy_qz(dm%dccp%ysz(1),              4, dm%dccp%ysz(3)) )! default y pencil
    allocate( dm%fbcz_qz(dm%dccp%zsz(1), dm%dccp%zsz(2),              4) )! default z pencil

    allocate( dm%fbcx_pr(             4, dm%dccc%xsz(2), dm%dccc%xsz(3)) )! default x pencil
    allocate( dm%fbcy_pr(dm%dccc%ysz(1),              4, dm%dccc%ysz(3)) )! default y pencil
    allocate( dm%fbcz_pr(dm%dccc%zsz(1), dm%dccc%zsz(2),              4) )! default z pencil

    if(dm%icoordinate == ICYLINDRICAL) then 
      allocate( dm%fbcy_qyr(dm%dcpc%ysz(1), 4,              dm%dcpc%ysz(3)) )
      allocate( dm%fbcz_qyr(dm%dcpc%zsz(1), dm%dcpc%zsz(2), 4             ) )
      allocate( dm%fbcy_qzr(dm%dccp%ysz(1), 4,              dm%dccp%ysz(3)) )
      allocate( dm%fbcz_qzr(dm%dccp%zsz(1), dm%dccp%zsz(2), 4             ) )
    end if

    return
  end subroutine 
!==========================================================================================================
!==========================================================================================================
  subroutine allocate_fbc_thermo(dm)
    type(t_domain), intent(inout) :: dm

    if( .not. dm%is_thermo) return

    allocate( dm%fbcx_gx(4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)) )! default x pencil
    allocate( dm%fbcx_gy(4, dm%dcpc%xsz(2), dm%dcpc%xsz(3)) )! default x pencil
    allocate( dm%fbcx_gz(4, dm%dccp%xsz(2), dm%dccp%xsz(3)) )! default x pencil

    allocate( dm%fbcy_gx(dm%dpcc%ysz(1), 4, dm%dpcc%ysz(3)) )! default y pencil
    allocate( dm%fbcy_gy(dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3)) )! default y pencil
    allocate( dm%fbcy_gz(dm%dccp%ysz(1), 4, dm%dccp%ysz(3)) )! default y pencil

    allocate( dm%fbcz_gx(dm%dpcc%zsz(1), dm%dpcc%zsz(2), 4) )! default z pencil
    allocate( dm%fbcz_gy(dm%dcpc%zsz(1), dm%dcpc%zsz(2), 4) )! default z pencil
    allocate( dm%fbcz_gz(dm%dccp%zsz(1), dm%dccp%zsz(2), 4) )! default z pencil

    if(dm%icoordinate == ICYLINDRICAL) then 
      allocate( dm%fbcy_gyr(dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3)) )
      allocate( dm%fbcy_gzr(dm%dccp%ysz(1), 4, dm%dccp%ysz(3)) )
      allocate( dm%fbcz_gyr(dm%dcpc%zsz(1), dm%dcpc%zsz(2), 4) )
      allocate( dm%fbcz_gzr(dm%dccp%zsz(1), dm%dccp%zsz(2), 4) )
    end if

    allocate( dm%fbcx_qw (4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)) )! default x pencil
    allocate( dm%fbcx_ftp(4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)) )! default x pencil

    allocate( dm%fbcy_qw (dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3)) )! default x pencil
    allocate( dm%fbcy_ftp(dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3)) )! default y pencil
    
    allocate( dm%fbcz_qw (dm%dccp%zsz(1), dm%dccp%zsz(2), 4)  )! default x pencil
    allocate( dm%fbcz_ftp(dm%dccp%zsz(1), dm%dccp%zsz(2), 4)  )! default z pencil

    return
  end subroutine 
!==========================================================================================================
!==========================================================================================================
  subroutine apply_fbc_given_flow (dm) ! apply once only
    type(t_domain), intent(inout)   :: dm

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
    call apply_fbcx_given_const(dm%fbcx_qx, dm%fbcx_const(:, 1))
    call apply_fbcx_given_const(dm%fbcx_qy, dm%fbcx_const(:, 2))
    call apply_fbcx_given_const(dm%fbcx_qz, dm%fbcx_const(:, 3))
    call apply_fbcx_given_const(dm%fbcx_pr, dm%fbcx_const(:, 4))
!----------------------------------------------------------------------------------------------------------
! y-bc in y-pencil, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
    call apply_fbcy_given_const(dm%fbcy_qx, dm%fbcy_const(:, 1))
    call apply_fbcy_given_const(dm%fbcy_qy, dm%fbcy_const(:, 2))
    call apply_fbcy_given_const(dm%fbcy_qz, dm%fbcy_const(:, 3)) ! geo_bc, rpi, not rci
    call apply_fbcy_given_const(dm%fbcy_pr, dm%fbcy_const(:, 4))
!----------------------------------------------------------------------------------------------------------
! z-bc in z-pencil, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
    call apply_fbcz_given_const(dm%fbcz_qx, dm%fbcz_const(:, 1))
    call apply_fbcz_given_const(dm%fbcz_qy, dm%fbcz_const(:, 2))
    call apply_fbcz_given_const(dm%fbcz_qz, dm%fbcz_const(:, 3))
    call apply_fbcz_given_const(dm%fbcz_pr, dm%fbcz_const(:, 4))
!==========================================================================================================
! y bc for primary variables, uy, uz
!==========================================================================================================
    if(dm%icoordinate == ICYLINDRICAL) then
!----------------------------------------------------------------------------------------------------------
! y-bc in y-pencil, qyr = qy/r = uy; qzr = qz/r = uz
!----------------------------------------------------------------------------------------------------------
      call apply_fbcy_given_const(dm%fbcy_qyr, dm%fbcy_const(:, 2), dm%rpi)
      call apply_fbcy_given_const(dm%fbcy_qzr, dm%fbcy_const(:, 3), dm%rci)
!----------------------------------------------------------------------------------------------------------
! z-bc in z-pencil, qyr = qy/r = uy; qzr = qz/r = uz
!----------------------------------------------------------------------------------------------------------
      call apply_fbcz_given_const(dm%fbcz_qyr, dm%fbcz_const(:, 2), dm%rpi, dm%dcpc%zst(2))
      call apply_fbcz_given_const(dm%fbcz_qzr, dm%fbcz_const(:, 3), dm%rci, dm%dccp%zst(2))
    end if

!==========================================================================================================
! to build up bc for var(x_const, y, z)
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
! x-bc1, qx(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 1) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%nc(2)
      filename(1) = trim('pf1d_u1y.dat') !(undim)
      call map_bc_1d_uprofile( filename(1), ny, dm%yc, var1y(1:ny) )
      call apply_fbcx_given_profile(dm%fbcx_qx, var1y, dm%dpcc%xst(2))
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, qy(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 2) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%np(2)
      filename(2) = trim('pf1d_v1y.dat') !(undim)
      call map_bc_1d_uprofile( filename(2), ny, dm%yp, var1y(1:ny) )
      if(dm%icoordinate == ICARTESIAN) var1y(1:ny) =  var1y(1:ny) / dm%rpi(1:ny)
      call apply_fbcx_given_profile(dm%fbcx_qy, var1y, dm%dcpc%xst(2))
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, qz(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 3) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%nc(2)
      filename(3) = trim('pf1d_w1y.dat') !(undim)
      call map_bc_1d_uprofile( filename(3), ny, dm%yc, var1y(1:ny) )
      if(dm%icoordinate == ICARTESIAN) var1y(1:ny) =  var1y(1:ny) / dm%rci(1:ny)
      call apply_fbcx_given_profile(dm%fbcx_qz, var1y, dm%dccp%xst(2))
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, pr(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 4) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%nc(2)
      filename(4) = trim('pf1d_p1y.dat') !(undim)
      call map_bc_1d_uprofile( filename(4), ny, dm%yc, var1y(1:ny) )
      call apply_fbcx_given_profile(dm%fbcx_pr, var1y, dm%dccc%xst(2))
    end if

    return
  end subroutine 

!==========================================================================================================
!==========================================================================================================
  subroutine apply_fbc_given_thermo(tm, dm) ! call this after scaling the fbc_ftp values
    use thermo_info_mod
    type(t_domain), intent(inout) :: dm
    type(t_thermo), intent(in)    :: tm

    real(WP) :: var1y(1:dm%np(2))
    
    integer :: ny, n, i, j, k
    real(WP) :: density
!==========================================================================================================
! to build up bc with constant values
! -3-1-||||-2-4
! for constant bc, 3=1= geometric bc, side 1;
!                  2=4= geometric bc, side 2 
!==========================================================================================================
!==========================================================================================================
! to build up bc for var(x_const, y, z)
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
! x-bc1, pr(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 5) == IBC_PROFILE1D) then !
      var1y = ZERO
      ny = dm%nc(2)
      filename(4) = trim('pf1d_T1y_undim.dat') !(undim)
      call map_bc_1d_uprofile( filename(4), ny, dm%yc, var1y(1:ny) )
      call apply_fbcx_given_profile(dm%fbcx_ftp(:,:,:)%t, var1y, dm%dccc%xst(2))
      call ftp_refresh_thermal_properties_from_T_undim_3d(dm%fbcx_ftp)
    end if
    
    do n = 1, 2 
!----------------------------------------------------------------------------------------------------------
! x-bc in x-pencil, ftp, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
      if(dm%ibcx_nominal(n, 5) == IBC_DIRICHLET) then
        dm%fbcx_ftp(n, :, :)%t   = dm%fbcx_const(n, 5)
        call ftp_refresh_thermal_properties_from_T_undim_3D(dm%fbcx_ftp)
        density = dm%fbcx_ftp(n, 1, 1)%d
      else if (dm%ibcx_nominal(n, 5) == IBC_NEUMANN) then
        dm%fbcx_qw(n, :, :) = dm%fbcx_const(n, 5) 
        dm%fbcx_ftp = tm%ftp_ini
        density = tm%ftp_ini%d
      else
        dm%fbcx_ftp = tm%ftp_ini
        density = tm%ftp_ini%d
      end if

      dm%fbcx_gx(n, :, :) = dm%fbcx_qx(n, :, :) * density
      dm%fbcx_gy(n, :, :) = dm%fbcx_qy(n, :, :) * density
      dm%fbcx_gz(n, :, :) = dm%fbcx_qz(n, :, :) * density
!----------------------------------------------------------------------------------------------------------
! y-bc in y-pencil, ftp, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
      if(dm%ibcy_nominal(n, 5) == IBC_DIRICHLET) then
        dm%fbcy_ftp(:, n, :)%t   = dm%fbcy_const(n, 5)
        call ftp_refresh_thermal_properties_from_T_undim_3D(dm%fbcy_ftp)
        density = dm%fbcy_ftp(1, n, 1)%d
      else if (dm%ibcy_nominal(n, 5) == IBC_NEUMANN) then
        dm%fbcy_qw(:, n, :) = dm%fbcy_const(n, 5) 
        dm%fbcy_ftp = tm%ftp_ini
        density = tm%ftp_ini%d
      else
        dm%fbcy_ftp = tm%ftp_ini
        density = tm%ftp_ini%d
      end if

      dm%fbcy_gx(:, n, :) = dm%fbcy_qx(:, n, :) * density
      dm%fbcy_gy(:, n, :) = dm%fbcy_qy(:, n, :) * density
      dm%fbcy_gz(:, n, :) = dm%fbcy_qz(:, n, :) * density
      !----------------------------------------------------------------------------------------------------------
      ! y-bc in y-pencil, qyr = qy/r = uy; qzr = qz/r = uz
      !----------------------------------------------------------------------------------------------------------
      if(dm%icoordinate == ICYLINDRICAL) then
        dm%fbcy_gyr(:, n, :) = dm%fbcy_qyr(:, n, :) * density
        dm%fbcy_gzr(:, n, :) = dm%fbcy_qzr(:, n, :) * density
      end if
!----------------------------------------------------------------------------------------------------------
! z-bc in z-pencil, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
      if( dm%ibcz_nominal(n, 5) == IBC_DIRICHLET ) then
        dm%fbcz_ftp(:, :, n)%t   = dm%fbcz_const(n, 5)
        call ftp_refresh_thermal_properties_from_T_undim_3D(dm%fbcz_ftp)
        density = dm%fbcz_ftp(1, 1, n)%d
      else if (dm%ibcz_nominal(n, 5) == IBC_NEUMANN) then
        dm%fbcz_qw(:, :, n) = dm%fbcz_const(n, 5) 
        dm%fbcz_ftp = tm%ftp_ini
        density = tm%ftp_ini%d
      else 
        dm%fbcz_ftp = tm%ftp_ini
        density = tm%ftp_ini%d
      end if
      dm%fbcz_gx(:, :, n) = dm%fbcz_qx(:, :, n) * density
      dm%fbcz_gy(:, :, n) = dm%fbcz_qy(:, :, n) * density
      dm%fbcz_gz(:, :, n) = dm%fbcz_qz(:, :, n) * density
!----------------------------------------------------------------------------------------------------------
! z-bc in z-pencil, qyr = qy/r = uy; qzr = qz/r = uz
!----------------------------------------------------------------------------------------------------------
      if(dm%icoordinate == ICYLINDRICAL) then
        dm%fbcz_gyr(:, :, n) = dm%fbcz_qyr(:, :, n) * density
        dm%fbcz_gzr(:, :, n) = dm%fbcz_qzr(:, :, n) * density
      end if
    end do

    return
  end subroutine 

!==========================================================================================================
!==========================================================================================================
  subroutine get_fbcy_circle_centre(var_xpencil, fbcy, ksym, dtmp)
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
  subroutine update_fbcy_cc_flow_halo(fl, dm)  ! for cylindrical only
    type(t_domain), intent(inout) :: dm
    type(t_flow), intent(in)      :: fl

    if(dm%icase /= ICASE_PIPE) return
    if(dm%icoordinate /= ICYLINDRICAL) return

!----------------------------------------------------------------------------------------------------------
!   ! qx bc in y - direction, interior
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy_qx(1) /= IBC_INTERIOR) call Print_error_msg('Error in ibcy_qx for the centre of the pipe.')
    call get_fbcy_circle_centre(fl%qx, dm%fbcy_qx, dm%knc_sym, dm%dpcc)
!----------------------------------------------------------------------------------------------------------
!   ! qy bc in y - direction, interior
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy_qy(1) /= IBC_INTERIOR) call Print_error_msg('Error in ibcy_qy for the centre of the pipe.') ! check, axial-symmetric at y=2?
    dm%fbcy_qy (:, 1, :) = ZERO
    dm%fbcy_qyr(:, 1, :) = ZERO
!----------------------------------------------------------------------------------------------------------
!   ! qz, gz bc in y - direction, interior
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy_qz(1) /= IBC_INTERIOR) call Print_error_msg('Error in ibcy_qz for the centre of the pipe.') ! 
    call get_fbcy_circle_centre(fl%qz, dm%fbcy_qz, dm%knc_sym, dm%dcpc)
    dm%fbcy_qzr(:, 1, :) = dm%fbcy_qz(:, 1, :) * dm%rci(1)
    dm%fbcy_qzr(:, 3, :) = dm%fbcy_qz(:, 3, :) * dm%rci(2)
!----------------------------------------------------------------------------------------------------------
!   ! pressure bc in y - direction, interior
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy_pr(1) /= IBC_INTERIOR) call Print_error_msg('Error in ibcy_pr for the centre of the pipe.') ! 
    call get_fbcy_circle_centre(fl%pres, dm%fbcy_pr, dm%knc_sym, dm%dccc)

    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine update_fbcy_cc_thermo_halo(fl, tm, dm)  ! for cylindrical only
    use thermo_info_mod
    type(t_domain), intent(inout) :: dm
    type(t_thermo), intent(in)    :: tm
    type(t_flow),   intent(in)    :: fl
    real(WP) :: fbcy(dm%dccc%ysz(1), 4, dm%dccc%ysz(3))
    integer :: i, j, k

    if(.not. dm%is_thermo) return
    if(dm%icase /= ICASE_PIPE) return
    if(dm%icoordinate /= ICYLINDRICAL) return
    
!----------------------------------------------------------------------------------------------------------
!   ! gx bc in y - direction, interior
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy_qx(1) /= IBC_INTERIOR) call Print_error_msg('Error in ibcy_gx for the centre of the pipe.')
    call get_fbcy_circle_centre(fl%gx, dm%fbcy_gx, dm%knc_sym, dm%dpcc)
!----------------------------------------------------------------------------------------------------------
!   ! qy, gy bc in y - direction, interior
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy_qy(1) /= IBC_DIRICHLET) call Print_error_msg('Error in ibcy_gy for the centre of the pipe.') 
    dm%fbcy_gy (:, 1, :) = ZERO
    dm%fbcy_gyr(:, 1, :) = ZERO
!----------------------------------------------------------------------------------------------------------
!   ! gz bc in y - direction, interior
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy_qz(1) /= IBC_INTERIOR) call Print_error_msg('Error in ibcy_qz for the centre of the pipe.') ! 
    call get_fbcy_circle_centre(fl%gz, dm%fbcy_gz, dm%knc_sym, dm%dcpc)
    dm%fbcy_gzr(:, 1, :) = dm%fbcy_gz(:, 1, :) * dm%rci(1)
    dm%fbcy_gzr(:, 3, :) = dm%fbcy_gz(:, 3, :) * dm%rci(2)
!----------------------------------------------------------------------------------------------------------
!   ! thermo bc in y - direction, interior
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy_Th(1) /= IBC_INTERIOR) call Print_error_msg('Error in ibcy_Th for the centre of the pipe.') !
    fbcy = dm%fbcy_ftp%t
    call get_fbcy_circle_centre(tm%tTemp, fbcy, dm%knc_sym, dm%dccc)
    dm%fbcy_ftp%t = fbcy
    call ftp_refresh_thermal_properties_from_T_undim_3D(dm%fbcy_ftp)

    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine apply_fbcx_2dm_halo(fbcx, var, iside, dtmp)
    type(DECOMP_INFO), intent(in) :: dtmp
    integer, intent(in)     :: iside
    real(WP), intent(in)    :: var(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3))
    real(WP), intent(inout) :: fbcx(:, :, :)
    
    if(iside == IFBC(1)) then !122
      fbcx(1, :, :) = var(dtmp%xsz(1),     :, :) ! interior there is a repeated shared points.
      fbcx(3, :, :) = var(dtmp%xsz(1) - 1, :, :)
    else if(iside == IFBC(2)) then !221
      fbcx(2, :, :) = var(1, :, :)
      fbcx(4, :, :) = var(2, :, :)
    else
      call Print_error_msg('Error input for apply_fbcx_2dm_halo')
    end if

    return
  end subroutine
!==========================================================================================================
  subroutine apply_fbcy_2dm_halo(fbcy, var, iside, dtmp)
    type(DECOMP_INFO), intent(in) :: dtmp
    integer, intent(in)     :: iside
    real(WP), intent(in)    :: var(dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3))
    real(WP), intent(inout) :: fbcy(:, :, :)

    real(WP), dimension( dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3) ) :: var_ypencil

    call transpose_x_to_y(var, var_ypencil, dtmp)
    if(iside == IFBC(1)) then
      fbcy(:, 1, :) = var_ypencil(:, dtmp%ysz(2),     :) ! interior there is a repeated shared points.
      fbcy(:, 3, :) = var_ypencil(:, dtmp%ysz(2) - 1, :)
    else if(iside == IFBC(2)) then
      fbcy(:, 2, :) = var_ypencil(:, 1, :)
      fbcy(:, 4, :) = var_ypencil(:, 2, :)
    else
      call Print_error_msg('Error input for apply_fbcy_2dm_halo')
    end if

    return
  end subroutine
!==========================================================================================================
  subroutine apply_fbcz_2dm_halo(fbcz, var, iside, dtmp)
    type(DECOMP_INFO), intent(in) :: dtmp
    integer, intent(in)     :: iside
    real(WP), intent(in)    :: var(dtmp%zsz(1), dtmp%zsz(2), dtmp%zsz(3))
    real(WP), intent(inout) :: fbcz(:, :, :)

    real(WP), dimension( dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3) ) :: var_ypencil
    real(WP), dimension( dtmp%zsz(1), dtmp%zsz(2), dtmp%zsz(3) ) :: var_zpencil

    call transpose_x_to_y(var,         var_ypencil, dtmp)
    call transpose_y_to_z(var_ypencil, var_zpencil, dtmp)
    if(iside == IFBC(1)) then
      fbcz(:, :, 1) = var_zpencil(:, :, dtmp%ysz(2)    ) ! interior there is a repeated shared points.
      fbcz(:, :, 3) = var_zpencil(:, :, dtmp%ysz(2) - 1)
    else if(iside == IFBC(2)) then
      fbcz(:, 2, :) = var_zpencil(:, :, 1)
      fbcz(:, 4, :) = var_zpencil(:, :, 2)
    else
      call Print_error_msg('Error input for apply_fbcz_2dm_halo')
    end if

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine apply_fbc_2dm_flow_halo(dm, fl, iside, bc_type)
    type(t_domain), intent(inout) :: dm
    type(t_flow),   intent(in)    :: fl
    integer, intent(in)           :: iside
    character(len=*), intent(in)  :: bc_type

    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: ac4c_ypencil
    real(WP), dimension( dm%dcpp%ysz(1), 4, dm%dcpp%ysz(3) ) :: ac4p_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), 4 ) :: acc4_zpencil
    real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), 4 ) :: acp4_zpencil
!----------------------------------------------------------------------------------------------------------
!   one repeat on nodes: 1'-2'-3'-4'-5'          no repeat of cells 1--2--3--4--5-
!                                    1'-2'-3'-4'-5'                               -1--2--3--4--5     
!                                 f2 f1 f2                                   f2 f1 f1 f2
!   iside = 122, 1-->2  : |-domain-1----3-1|-domain-2-----|
!   iside = 221, 1<--2  : |-domain-1---|-2-4-domain-2-----|
!----------------------------------------------------------------------------------------------------------
    select case (bc_type)
    case ('x')
        if (dm%ibcx_qx(iside) == IBC_INTERIOR) then
            call apply_fbcx_2dm_halo(dm%fbcx_qx, fl%qx, iside, dm%dpcc)
            if (dm%is_thermo) &
            call apply_fbcx_2dm_halo(dm%fbcx_gx, fl%gx, iside, dm%dpcc)
        end if
        if (dm%ibcx_qy(iside) == IBC_INTERIOR) then
            call apply_fbcx_2dm_halo(dm%fbcx_qy, fl%qy, iside, dm%dcpc)
            if (dm%is_thermo) &
            call apply_fbcx_2dm_halo(dm%fbcx_gy, fl%gy, iside, dm%dcpc)
        end if
        if (dm%ibcx_qz(iside) == IBC_INTERIOR) then
            call apply_fbcx_2dm_halo(dm%fbcx_qz, fl%qz, iside, dm%dccp)
            if (dm%is_thermo) &
            call apply_fbcx_2dm_halo(dm%fbcx_gz, fl%gz, iside, dm%dccp)
        end if
        if (dm%ibcx_pr(iside) == IBC_INTERIOR) then
            call apply_fbcx_2dm_halo(dm%fbcx_pr, fl%pres, iside, dm%dccc)
        end if
    case ('y')
        if (dm%ibcy_qx(iside) == IBC_INTERIOR) then
            call apply_fbcy_2dm_halo(dm%fbcy_qx, fl%qx, iside, dm%dpcc)
            if (dm%is_thermo) &
            call apply_fbcy_2dm_halo(dm%fbcy_gx, fl%gx, iside, dm%dpcc)
        end if
        if (dm%ibcy_qy(iside) == IBC_INTERIOR) then
            call apply_fbcy_2dm_halo(dm%fbcy_qy, fl%qy, iside, dm%dcpc)
            if (dm%is_thermo) &
            call apply_fbcy_2dm_halo(dm%fbcy_gy, fl%gy, iside, dm%dcpc)
            if (dm%icoordinate == ICYLINDRICAL) then
                ac4c_ypencil = dm%fbcy_qy
                call multiple_cylindrical_rn_x4x(ac4c_ypencil, dm%dcpc, dm%rpi, 1, IPENCIL(2))
                dm%fbcy_qyr = ac4c_ypencil
                if (dm%is_thermo) then
                    ac4c_ypencil = dm%fbcy_gy
                    call multiple_cylindrical_rn_x4x(ac4c_ypencil, dm%dcpc, dm%rpi, 1, IPENCIL(2))
                    dm%fbcy_gyr = ac4c_ypencil
                end if
            end if
        end if
        if (dm%ibcy_qz(iside) == IBC_INTERIOR) then
            call apply_fbcy_2dm_halo(dm%fbcy_qz, fl%qz, iside, dm%dccp)
            if (dm%is_thermo) &
            call apply_fbcy_2dm_halo(dm%fbcy_gz, fl%gz, iside, dm%dccp)
            if (dm%icoordinate == ICYLINDRICAL) then
                ac4p_ypencil = dm%fbcy_qz
                call multiple_cylindrical_rn_x4x(ac4p_ypencil, dm%dccp, dm%rci, 1, IPENCIL(2))
                dm%fbcy_qzr = ac4p_ypencil
                if (dm%is_thermo) then
                  ac4p_ypencil = dm%fbcy_gz
                  call multiple_cylindrical_rn_x4x(ac4p_ypencil, dm%dccp, dm%rci, 1, IPENCIL(2))
                  dm%fbcy_gzr = ac4p_ypencil
                end if
            end if
        end if
        if (dm%ibcy_pr(iside) == IBC_INTERIOR) then
            call apply_fbcy_2dm_halo(dm%fbcy_pr, fl%pres, iside, dm%dccc)
        end if
    case ('z')
        if (dm%ibcz_qx(iside) == IBC_INTERIOR) then
            call apply_fbcz_2dm_halo(dm%fbcz_qx, fl%qx, iside, dm%dpcc)
            if (dm%is_thermo) &
            call apply_fbcz_2dm_halo(dm%fbcz_gx, fl%gx, iside, dm%dpcc)
        end if
        if (dm%ibcz_qy(iside) == IBC_INTERIOR) then
            call apply_fbcz_2dm_halo(dm%fbcz_qy, fl%qy, iside, dm%dcpc)
            if (dm%is_thermo) &
            call apply_fbcz_2dm_halo(dm%fbcz_gy, fl%gy, iside, dm%dcpc)
            if (dm%icoordinate == ICYLINDRICAL) then
                acp4_zpencil = dm%fbcz_qy
                call multiple_cylindrical_rn_xx4(acp4_zpencil, dm%dcpc, dm%rpi, 1, IPENCIL(3))
                dm%fbcz_qyr = acp4_zpencil
                if (dm%is_thermo) then
                  acp4_zpencil = dm%fbcz_gy
                  call multiple_cylindrical_rn_xx4(acp4_zpencil, dm%dcpc, dm%rpi, 1, IPENCIL(3))
                  dm%fbcz_gyr = acp4_zpencil
                end if
            end if
        end if
        if (dm%ibcz_qz(iside) == IBC_INTERIOR) then
            call apply_fbcz_2dm_halo(dm%fbcz_qz, fl%qz, iside, dm%dccp)
            if (dm%is_thermo) &
            call apply_fbcz_2dm_halo(dm%fbcz_gz, fl%gz, iside, dm%dccp)
            if (dm%icoordinate == ICYLINDRICAL) then
                acc4_zpencil = dm%fbcz_qz
                call multiple_cylindrical_rn_xx4(acc4_zpencil, dm%dccp, dm%rci, 1, IPENCIL(3))
                dm%fbcz_qzr = acc4_zpencil
                if (dm%is_thermo) then
                  acc4_zpencil = dm%fbcz_gz
                  call multiple_cylindrical_rn_xx4(acc4_zpencil, dm%dccp, dm%rci, 1, IPENCIL(3))
                  dm%fbcz_gzr = acc4_zpencil
                end if
            end if
        end if
        if (dm%ibcz_pr(iside) == IBC_INTERIOR) then
            call apply_fbcz_2dm_halo(dm%fbcz_pr, fl%pres, iside, dm%dccc)
        end if
    end select
  end subroutine apply_fbc_2dm_flow_halo
!==========================================================================================================
!==========================================================================================================
  subroutine update_fbc_2dm_flow_halo(dm1, fl1, dm2, fl2)
    type(t_domain), intent(inout) :: dm1, dm2
    type(t_flow),   intent(in)    :: fl1, fl2

    integer :: n
    ! x-boundary conditions
    call apply_fbc_2dm_flow_halo(dm2, fl1, IFBC(1), 'x')
    call apply_fbc_2dm_flow_halo(dm1, fl2, IFBC(2), 'x')
    ! y-boundary conditions
    call apply_fbc_2dm_flow_halo(dm2, fl1, IFBC(1), 'y')
    call apply_fbc_2dm_flow_halo(dm1, fl2, IFBC(2), 'y')
    ! z-boundary conditions
    call apply_fbc_2dm_flow_halo(dm2, fl1, IFBC(1), 'z')
    call apply_fbc_2dm_flow_halo(dm1, fl2, IFBC(2), 'z')

    ! for turb inlet only, no thermal inlet
    if (dm2%is_thermo .and. (.not.dm1%is_thermo)) then
      do n = 1, 3, 2
        ! x-boundary conditions
        if(dm2%ibcx_qx(1) == IBC_INTERIOR) &
        dm2%fbcx_gx(n, :, :) = dm2%fbcx_qx(n, :, :) * dm2%fbcx_ftp(1, :, :)%d
        if(dm2%ibcx_qy(1) == IBC_INTERIOR) &
        dm2%fbcx_gy(n, :, :) = dm2%fbcx_qy(n, :, :) * dm2%fbcx_ftp(1, :, :)%d
        if(dm2%ibcx_qz(1) == IBC_INTERIOR) &
        dm2%fbcx_gz(n, :, :) = dm2%fbcx_qz(n, :, :) * dm2%fbcx_ftp(1, :, :)%d
        ! x-boundary conditions
        if(dm2%ibcy_qx(1) == IBC_INTERIOR) &
        dm2%fbcy_gx(:, n, :) = dm2%fbcy_qx(:, n, :) * dm2%fbcy_ftp(:, 1, :)%d
        if(dm2%ibcy_qy(1) == IBC_INTERIOR) then
          dm2%fbcy_gy(:, n, :) = dm2%fbcy_qy(:, n, :) * dm2%fbcy_ftp(:, 1, :)%d
          if(dm2%icoordinate == ICYLINDRICAL) &
          dm2%fbcy_gyr(:, n, :)= dm2%fbcy_qyr(:, n, :) * dm2%fbcy_ftp(:, 1, :)%d
        end if
        if(dm2%ibcy_qz(1) == IBC_INTERIOR) then
          dm2%fbcy_gz(:, n, :) = dm2%fbcy_qz(:, n, :) * dm2%fbcy_ftp(:, 1, :)%d
          if(dm2%icoordinate == ICYLINDRICAL) &
          dm2%fbcy_gzr(:, n, :)= dm2%fbcy_qzr(:, n, :) * dm2%fbcy_ftp(:, 1, :)%d
        end if
        ! z-boundary conditions
        if(dm2%ibcz_qx(1) == IBC_INTERIOR) &
        dm2%fbcz_gx(:, :, n) = dm2%fbcz_qx(:, :, n) * dm2%fbcz_ftp(:, :, 1)%d
        if(dm2%ibcz_qy(1) == IBC_INTERIOR) then
          dm2%fbcz_gy(:, :, n) = dm2%fbcz_qy(:, :, n) * dm2%fbcz_ftp(:, :, 1)%d
          if(dm2%icoordinate == ICYLINDRICAL) &
          dm2%fbcz_gyr(:, :, n) = dm2%fbcz_qyr(:, :, n) * dm2%fbcz_ftp(:, :, 1)%d
        end if
        if(dm2%ibcz_qz(1) == IBC_INTERIOR) then
          dm2%fbcz_gz(:, :, n) = dm2%fbcz_qz(:, :, n) * dm2%fbcz_ftp(:, :, 1)%d
          if(dm2%icoordinate == ICYLINDRICAL) &
          dm2%fbcz_gzr(:, :, n) = dm2%fbcz_qzr(:, :, n) * dm2%fbcz_ftp(:, :, 1)%d
        end if
      end do
    end if
    
    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine update_fbc_2dm_thermo_halo(dm1, tm1, dm2, tm2)
    use thermo_info_mod
    type(t_domain), intent(inout) :: dm1, dm2
    type(t_thermo), intent(in)    :: tm1, tm2
    
    integer :: i, j, k
    real(WP), dimension( dm1%dccc%ysz(1), dm1%dccc%ysz(2), dm1%dccc%ysz(3) ) :: accc0_ypencil
    real(WP), dimension( dm1%dccc%zsz(1), dm1%dccc%zsz(2), dm1%dccc%zsz(3) ) :: accc0_zpencil
    real(WP), dimension( dm2%dccc%ysz(1), dm2%dccc%ysz(2), dm2%dccc%ysz(3) ) :: accc1_ypencil
    real(WP), dimension( dm2%dccc%zsz(1), dm2%dccc%zsz(2), dm2%dccc%zsz(3) ) :: accc1_zpencil
!----------------------------------------------------------------------------------------------------------
!   bc in x - direction
!----------------------------------------------------------------------------------------------------------
    ! thermal field, dm1-dm2
    if(dm2%ibcx_Th(1) == IBC_INTERIOR) then
      dm2%fbcx_ftp(1, :, :)%t = tm1%tTemp(dm1%nc(1),     :, :)
      dm2%fbcx_ftp(3, :, :)%t = tm1%tTemp(dm1%nc(1) - 1, :, :)
      call ftp_refresh_thermal_properties_from_T_undim_3D(dm2%fbcx_ftp)
    end if

    if(dm1%ibcx_Th(2) == IBC_INTERIOR) then
      dm1%fbcx_ftp(2, :, :)%t = tm2%tTemp(1, :, :)
      dm1%fbcx_ftp(4, :, :)%t = tm2%tTemp(2, :, :)
      call ftp_refresh_thermal_properties_from_T_undim_3D(dm1%fbcx_ftp)
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in y - direction
!----------------------------------------------------------------------------------------------------------
    ! thermal field, dm1-dm2
    if(dm2%ibcy_Th(1) == IBC_INTERIOR) then
      call transpose_x_to_y(tm1%tTemp, accc0_ypencil, dm1%dccc)
      dm2%fbcy_ftp(:, 1, :)%t = accc0_ypencil(:, dm1%nc(1),     :)
      dm2%fbcy_ftp(:, 3, :)%t = accc0_ypencil(:, dm1%nc(1) - 1, :)
      call ftp_refresh_thermal_properties_from_T_undim_3D(dm2%fbcy_ftp)
    end if

    if(dm1%ibcy_Th(2) == IBC_INTERIOR) then
      call transpose_x_to_y(tm2%tTemp, accc1_ypencil, dm2%dccc)
      dm1%fbcy_ftp(:, 2, :)%t = accc1_ypencil(:, 1, :)
      dm1%fbcy_ftp(:, 4, :)%t = accc1_ypencil(:, 2, :)
      call ftp_refresh_thermal_properties_from_T_undim_3D(dm1%fbcy_ftp)
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in z - direction
!----------------------------------------------------------------------------------------------------------
    ! thermal field, dm1-dm2
    if(dm2%ibcz_Th(1) == IBC_INTERIOR) then
      call transpose_x_to_y(tm1%tTemp,     accc0_ypencil, dm1%dccc)
      call transpose_y_to_z(accc0_ypencil, accc0_zpencil, dm1%dccc)
      dm2%fbcz_ftp(:, :, 1)%t = accc0_zpencil(:, :, dm1%nc(1)    )
      dm2%fbcz_ftp(:, :, 3)%t = accc0_zpencil(:, :, dm1%nc(1) - 1)
      call ftp_refresh_thermal_properties_from_T_undim_3D(dm2%fbcz_ftp)
    end if

    if(dm1%ibcz_Th(2) == IBC_INTERIOR) then
      call transpose_x_to_y(tm2%tTemp,     accc1_ypencil, dm2%dccc)
      call transpose_y_to_z(accc1_ypencil, accc1_zpencil, dm2%dccc)
      dm1%fbcz_ftp(:, :, 2)%t = accc1_zpencil(:, :, 1)
      dm1%fbcz_ftp(:, :, 4)%t = accc1_zpencil(:, :, 2)
      call ftp_refresh_thermal_properties_from_T_undim_3D(dm1%fbcz_ftp)
    end if

    return
  end subroutine

!==========================================================================================================
! to calculate boundary during calculation from primary boundary
  subroutine reconstruct_symmetry_ibc(ibc, mbc, jbc)
    integer, intent(in)  :: ibc(2)
    integer, intent(out) :: mbc(2, 3)
    integer, intent(in), optional :: jbc(2)
    
    integer :: i
    
    mbc(:, JBC_SELF) = ibc(:)
    mbc(:, JBC_GRAD) = ibc(:)
    mbc(:, JBC_PROD) = ibc(:)

    do i = 1, 2
      if(present(jbc)) then
        
        if(ibc(i)==IBC_SYMMETRIC .and. jbc(i)==IBC_SYMMETRIC) then
          mbc(i, JBC_PROD) = IBC_SYMMETRIC
        else if (ibc(i)==IBC_SYMMETRIC .and. jbc(i)==IBC_ASYMMETRIC) then
          mbc(i, JBC_PROD) = IBC_ASYMMETRIC
        else if (ibc(i)==IBC_ASYMMETRIC .and. jbc(i)==IBC_SYMMETRIC) then
          mbc(i, JBC_PROD) = IBC_ASYMMETRIC
        else if (ibc(i)==IBC_ASYMMETRIC .and. jbc(i)==IBC_ASYMMETRIC) then
          mbc(i, JBC_PROD) = IBC_SYMMETRIC
        else 
          if(ibc(i)/=jbc(i)) then
            if(nrank==0) write(*, *) "BCs for the side ", i, " are ", ibc(i), jbc(i) 
            call Print_warning_msg("The two operational variables have different boundary conditions.")
          else
            mbc(i, :) = ibc(i)
          end if
        end if

      else

        if(ibc(i)==IBC_SYMMETRIC) then
          mbc(i, JBC_SELF) = ibc(i)               ! variable itself
          mbc(i, JBC_GRAD) = IBC_ASYMMETRIC       ! d(var)/dn, 
          mbc(i, JBC_PROD) = ibc(i)               ! var * var
        else if(ibc(i)==IBC_ASYMMETRIC) then
          mbc(i, JBC_SELF) = ibc(i)              ! variable itself
          mbc(i, JBC_GRAD) = IBC_SYMMETRIC       ! d(var)/dn, 
          mbc(i, JBC_PROD) = IBC_SYMMETRIC       ! var * var
        else
          mbc(i, :) = ibc(i)
        end if

      end if 

    end do


    return
  end subroutine 
!==========================================================================================================
  subroutine config_calc_eqs_ibc(dm)
    use wtformat_mod
    type(t_domain), intent(inout)   :: dm
    
    integer :: mbc(2, 3), mbc0(2, 3)
    integer :: bc(2)
!----------------------------------------------------------------------------------------------------------
!   x-mom
!----------------------------------------------------------------------------------------------------------
    call reconstruct_symmetry_ibc(dm%ibcx_qx, mbc, dm%ibcx_qx)
    mbcx_cov1(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom x-convection is ", mbcx_cov1 

    call reconstruct_symmetry_ibc(dm%ibcy_qy, mbc, dm%ibcy_qx)
    mbcy_cov1(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom y-convection is ", mbcy_cov1

    call reconstruct_symmetry_ibc(dm%ibcz_qz, mbc, dm%ibcz_qx)
    mbcz_cov1(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom z-convection is ", mbcz_cov1

    call reconstruct_symmetry_ibc(dm%ibcx_qx, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call reconstruct_symmetry_ibc(dm%ibcx_Th, mbc, bc)
    mbcx_tau1(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom x-diffusion  is ", mbcx_tau1

    call reconstruct_symmetry_ibc(dm%ibcy_qx, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call reconstruct_symmetry_ibc(dm%ibcy_Th, mbc, bc)
    call reconstruct_symmetry_ibc(dm%ibcy_Th, mbc0, dm%ibcy_qy)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcy_tau1 is wrong.")
    mbcy_tau1(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom y-diffusion  is ", mbcy_tau1

    call reconstruct_symmetry_ibc(dm%ibcz_qx, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call reconstruct_symmetry_ibc(dm%ibcz_Th, mbc, bc)
    call reconstruct_symmetry_ibc(dm%ibcz_Th, mbc0, dm%ibcz_qz)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcy_tau1 is wrong.")
    mbcz_tau1(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom z-diffusion  is ", mbcz_tau1
!----------------------------------------------------------------------------------------------------------
!   y-mom
!----------------------------------------------------------------------------------------------------------
    call reconstruct_symmetry_ibc(dm%ibcx_qx, mbc, dm%ibcx_qy)
    mbcx_cov2(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom x-convection is ", mbcx_cov2

    call reconstruct_symmetry_ibc(dm%ibcy_qy, mbc, dm%ibcy_qy)
    mbcy_cov2(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom y-convection is ", mbcy_cov2

    call reconstruct_symmetry_ibc(dm%ibcz_qz, mbc, dm%ibcz_qy)
    mbcz_cov2(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom z-convection is ", mbcz_cov2

    if(dm%icoordinate == ICYLINDRICAL) then
      call reconstruct_symmetry_ibc(dm%ibcy_qz, mbc, dm%ibcy_qz)
      mbcr_cov2(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom r-convection is ", mbcr_cov2
    end if

    call reconstruct_symmetry_ibc(dm%ibcx_qy, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call reconstruct_symmetry_ibc(dm%ibcx_Th, mbc, bc)
    call reconstruct_symmetry_ibc(dm%ibcx_Th, mbc0, dm%ibcx_qx)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcx_tau2 is wrong.")
    mbcx_tau2(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom x-diffusion  is ", mbcx_tau2

    call reconstruct_symmetry_ibc(dm%ibcy_qy, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call reconstruct_symmetry_ibc(dm%ibcy_Th, mbc, bc)
    mbcy_tau2(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom y-diffusion  is ", mbcy_tau2

    call reconstruct_symmetry_ibc(dm%ibcz_qy, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call reconstruct_symmetry_ibc(dm%ibcz_Th, mbc, bc)
    call reconstruct_symmetry_ibc(dm%ibcz_Th, mbc0, dm%ibcz_qz)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcz_tau2 is wrong.")
    mbcz_tau2(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom z-diffusion  is ", mbcz_tau2

    if(dm%icoordinate == ICYLINDRICAL) then
      call reconstruct_symmetry_ibc(dm%ibcy_qz, mbc, dm%ibcy_Th)
      mbcr_tau2(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom r-diffusion  is ", mbcr_tau2
    end if
!----------------------------------------------------------------------------------------------------------
!   z-mom
!----------------------------------------------------------------------------------------------------------
    call reconstruct_symmetry_ibc(dm%ibcx_qx, mbc, dm%ibcx_qz)
    mbcx_cov3(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom x-convection is ", mbcx_cov3

    call reconstruct_symmetry_ibc(dm%ibcy_qy, mbc, dm%ibcy_qz)
    mbcy_cov3(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom y-convection is ", mbcy_cov3

    call reconstruct_symmetry_ibc(dm%ibcz_qz, mbc, dm%ibcz_qz)
    mbcz_cov3(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom z-convection is ", mbcz_cov3

    if(dm%icoordinate == ICYLINDRICAL) then
      call reconstruct_symmetry_ibc(dm%ibcy_qy, mbc, dm%ibcy_qz)
      mbcr_cov3(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom r-convection is ", mbcr_cov3
    end if

    call reconstruct_symmetry_ibc(dm%ibcx_qz, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call reconstruct_symmetry_ibc(dm%ibcx_Th, mbc, bc)
    call reconstruct_symmetry_ibc(dm%ibcx_Th, mbc0, dm%ibcx_qx)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcx_tau3 is wrong.")
    mbcx_tau3(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom x-diffusion  is ", mbcx_tau3

    call reconstruct_symmetry_ibc(dm%ibcy_qz, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call reconstruct_symmetry_ibc(dm%ibcy_Th, mbc, bc)
    call reconstruct_symmetry_ibc(dm%ibcy_Th, mbc0, dm%ibcy_qy)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcy_tau3 is wrong.")
    mbcy_tau3(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom y-diffusion  is ", mbcy_tau3

    call reconstruct_symmetry_ibc(dm%ibcz_qz, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call reconstruct_symmetry_ibc(dm%ibcz_Th, mbc, bc)
    mbcz_tau3 = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom z-diffusion  is ", mbcz_tau3

    if(dm%icoordinate == ICYLINDRICAL) then
      call reconstruct_symmetry_ibc(dm%ibcy_Th, mbc, dm%ibcy_qz)
      if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcy_tau3 is wrong.")
      mbcr_tau3(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom r-diffusion  is ", mbcr_tau3
    end if

!----------------------------------------------------------------------------------------------------------
!   energy-eqs
!----------------------------------------------------------------------------------------------------------
    if(dm%is_thermo) then
      call reconstruct_symmetry_ibc(dm%ibcx_qx, mbc, dm%ibcx_Th)
      ebcx_conv(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for energy x-convection is ", ebcx_conv

      call reconstruct_symmetry_ibc(dm%ibcy_qy, mbc, dm%ibcy_Th)
      ebcy_conv(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for energy y-convection is ", ebcy_conv

      call reconstruct_symmetry_ibc(dm%ibcz_qz, mbc, dm%ibcz_Th)
      ebcz_conv(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for energy z-convection is ", ebcz_conv

      call reconstruct_symmetry_ibc(dm%ibcx_Th, mbc)
      bc(1:2) = mbc(1:2, JBC_GRAD)
      call reconstruct_symmetry_ibc(dm%ibcx_Th, mbc, bc)
      ebcx_difu = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for energy x-diffusion is ", ebcx_difu

      call reconstruct_symmetry_ibc(dm%ibcy_Th, mbc)
      bc(1:2) = mbc(1:2, JBC_GRAD)
      call reconstruct_symmetry_ibc(dm%ibcy_Th, mbc, bc)
      ebcy_difu(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for energy y-diffusion is ", ebcy_difu

      call reconstruct_symmetry_ibc(dm%ibcz_Th, mbc)
      bc(1:2) = mbc(1:2, JBC_GRAD)
      call reconstruct_symmetry_ibc(dm%ibcz_Th, mbc, bc)
      ebcz_difu(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for energy z-diffusion is ", ebcz_difu

    end if

!----------------------------------------------------------------------------------------------------------
! preparation for b.c.
!----------------------------------------------------------------------------------------------------------
    is_fbcx_velo_required = .false.
    if(dm%ibcx_qx(1) == IBC_DIRICHLET .or. &
       dm%ibcx_qx(2) == IBC_DIRICHLET .or. &
       dm%ibcx_qy(1) == IBC_DIRICHLET .or. &
       dm%ibcx_qy(2) == IBC_DIRICHLET .or. &
       dm%ibcx_qz(1) == IBC_DIRICHLET .or. &
       dm%ibcx_qz(2) == IBC_DIRICHLET ) then
       is_fbcx_velo_required = .true.
      ! to add neumann later, check
    end if
    is_fbcy_velo_required = .false.
    if(dm%ibcy_qx(1) == IBC_DIRICHLET .or. &
       dm%ibcy_qx(2) == IBC_DIRICHLET .or. &
       dm%ibcy_qy(1) == IBC_DIRICHLET .or. &
       dm%ibcy_qy(2) == IBC_DIRICHLET .or. &
       dm%ibcy_qz(1) == IBC_DIRICHLET .or. &
       dm%ibcy_qz(2) == IBC_DIRICHLET ) then
       is_fbcy_velo_required = .true.
      ! to add neumann later, check
    end if
    is_fbcz_velo_required = .false.
    if(dm%ibcz_qx(1) == IBC_DIRICHLET .or. &
       dm%ibcz_qx(2) == IBC_DIRICHLET .or. &
       dm%ibcz_qy(1) == IBC_DIRICHLET .or. &
       dm%ibcz_qy(2) == IBC_DIRICHLET .or. &
       dm%ibcz_qz(1) == IBC_DIRICHLET .or. &
       dm%ibcz_qz(2) == IBC_DIRICHLET ) then
       is_fbcy_velo_required = .true.
      ! to add neumann later, check
    end if

    return 
  end subroutine


end module
