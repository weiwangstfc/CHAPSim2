module bc_dirichlet_mod
  use udf_type_mod
  use parameters_constant_mod
  use print_msg_mod
  implicit none
  character(18) :: filename(5)


  private :: map_bc_1d_uprofile     
  public  :: initialise_fbcx_given_profile 

  private :: initialise_fbcx_given_const
  private :: initialise_fbcy_given_const
  private :: initialise_fbcz_given_const
  public  :: initialise_fbc_flow_given   ! applied once only, for bc of constant velocity
  public  :: initialise_fbc_thermo_given ! applied once only, for bc of constant temperature
  public  :: enforce_var_from_const
  
  

contains

!==========================================================================================================
!==========================================================================================================
  subroutine  map_bc_1d_uprofile(filename, n, y, u)
    use io_files_mod
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
    integer :: pf_unit, j
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

    if(nrank == 0) then
      open ( newunit = pf_unit,     &
              file    = trim(dir_chkp)//'/check_given_ux_profile.dat', &
              status  = 'replace',         &
              action  = 'write')
      write(pf_unit, '(A)') "#j, y, u - original"
      do j = 1, nn
        write(pf_unit, '(1I3.1, 5ES15.7)') j, yy(j), uprofile(j)
      end do
      write(pf_unit, '(A)') "#j, y, u - interpolation"
      do j = 1, n
        write(pf_unit, '(1I3.1, 5ES15.7)') j, y(j), u(j)
      end do
      close(pf_unit)
    end if

    deallocate(uprofile)
    deallocate(yy)

    return
  end subroutine

  !==========================================================================================================
  !==========================================================================================================
  subroutine initialise_fbcx_given_profile(fbcx, var1y, jst, str)
    use io_files_mod
    implicit none
    real(WP), intent(inout) :: fbcx(:, :, :)
    real(WP), intent(in)    :: var1y(:)
    integer,  intent(in)    :: jst
    character(2), intent(in)   :: str

    integer :: k, j, jj
    integer :: pf_unit

    do k = 1, size(fbcx, 3)
      do j = 1, size(fbcx, 2)
        jj = jst + j - 1
        fbcx(1, j, k) = var1y(jj)
        fbcx(3, j, k) = fbcx(1, j, k)
      end do
    end do

    open ( newunit = pf_unit,     &
            file    = trim(dir_chkp)//'/check_given_'//trim(str)//'_profile.dat', &
            position= 'append',         &
            action  = 'write')
    write(pf_unit, '(A)') "#fbcx"
    do j = 1, size(fbcx, 2)
      write(pf_unit, '(1I3.1, 1ES15.7)') j, fbcx(1, j, 1)
    end do
    close(pf_unit)

    return
  end subroutine
  !==========================================================================================================
  subroutine initialise_fbcx_given_const(fbcx, fbcx_const)
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
  subroutine initialise_fbcy_given_const(fbcy, fbcy_const, ri)
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
  subroutine initialise_fbcz_given_const(fbcz, fbcz_const, ri, jst)
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
!==========================================================================================================
  subroutine initialise_fbc_flow_given (dm) ! apply once only
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
    call initialise_fbcx_given_const(dm%fbcx_qx, dm%fbcx_const(:, 1))
    call initialise_fbcx_given_const(dm%fbcx_qy, dm%fbcx_const(:, 2))
    call initialise_fbcx_given_const(dm%fbcx_qz, dm%fbcx_const(:, 3))
    call initialise_fbcx_given_const(dm%fbcx_pr, dm%fbcx_const(:, 4))
!----------------------------------------------------------------------------------------------------------
! y-bc in y-pencil, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
    call initialise_fbcy_given_const(dm%fbcy_qx, dm%fbcy_const(:, 1))
    call initialise_fbcy_given_const(dm%fbcy_qy, dm%fbcy_const(:, 2))
    call initialise_fbcy_given_const(dm%fbcy_qz, dm%fbcy_const(:, 3)) ! geo_bc, rpi, not rci
    call initialise_fbcy_given_const(dm%fbcy_pr, dm%fbcy_const(:, 4))
!----------------------------------------------------------------------------------------------------------
! z-bc in z-pencil, qx, qy, qz, pr
!----------------------------------------------------------------------------------------------------------
    call initialise_fbcz_given_const(dm%fbcz_qx, dm%fbcz_const(:, 1))
    call initialise_fbcz_given_const(dm%fbcz_qy, dm%fbcz_const(:, 2))
    call initialise_fbcz_given_const(dm%fbcz_qz, dm%fbcz_const(:, 3))
    call initialise_fbcz_given_const(dm%fbcz_pr, dm%fbcz_const(:, 4))
!==========================================================================================================
! y bc for primary variables, uy, uz
!==========================================================================================================
    if(dm%icoordinate == ICYLINDRICAL) then
!----------------------------------------------------------------------------------------------------------
! y-bc in y-pencil, qyr = qy/r = uy; qzr = qz/r = uz
!----------------------------------------------------------------------------------------------------------
      call initialise_fbcy_given_const(dm%fbcy_qyr, dm%fbcy_const(:, 2), dm%rpi)
      call initialise_fbcy_given_const(dm%fbcy_qzr, dm%fbcy_const(:, 3), dm%rci)
!----------------------------------------------------------------------------------------------------------
! z-bc in z-pencil, qyr = qy/r = uy; qzr = qz/r = uz
!----------------------------------------------------------------------------------------------------------
      call initialise_fbcz_given_const(dm%fbcz_qyr, dm%fbcz_const(:, 2), dm%rpi, dm%dcpc%zst(2))
      call initialise_fbcz_given_const(dm%fbcz_qzr, dm%fbcz_const(:, 3), dm%rci, dm%dccp%zst(2))
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
      filename(1) = trim('PF1D_U1Y.DAT') !(undim)
      ! call map_bc_1d_uprofile( filename(1), ny, dm%yc, var1y(1:ny) )
      ! call initialise_fbcx_given_profile(dm%fbcx_qx, var1y, dm%dpcc%xst(2), 'qx')
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, qy(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 2) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%np(2)
      filename(2) = trim('PF1D_V1Y.DAT') !(undim)
      call map_bc_1d_uprofile( filename(2), ny, dm%yp, var1y(1:ny) )
      if(dm%icoordinate == ICARTESIAN) var1y(1:ny) =  var1y(1:ny) / dm%rpi(1:ny)
      call initialise_fbcx_given_profile(dm%fbcx_qy, var1y, dm%dcpc%xst(2), 'qy')
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, qz(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 3) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%nc(2)
      filename(3) = trim('PF1D_W1Y.DAT') !(undim)
      call map_bc_1d_uprofile( filename(3), ny, dm%yc, var1y(1:ny) )
      if(dm%icoordinate == ICARTESIAN) var1y(1:ny) =  var1y(1:ny) / dm%rci(1:ny)
      call initialise_fbcx_given_profile(dm%fbcx_qz, var1y, dm%dccp%xst(2), 'qz')
    end if
!----------------------------------------------------------------------------------------------------------
! x-bc1, pr(x_c, y, z)
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_nominal(1, 4) == IBC_PROFILE1D) then
      var1y = ZERO
      ny = dm%nc(2)
      filename(4) = trim('PF1D_P1Y.DAT') !(undim)
      call map_bc_1d_uprofile( filename(4), ny, dm%yc, var1y(1:ny) )
      call initialise_fbcx_given_profile(dm%fbcx_pr, var1y, dm%dccc%xst(2), 'pr')
    end if

    return
  end subroutine 

!==========================================================================================================
!==========================================================================================================
  subroutine initialise_fbc_thermo_given(tm, dm) ! call this after scaling the fbc_ftp values
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
      call initialise_fbcx_given_profile(dm%fbcx_ftp(:,:,:)%t, var1y, dm%dccc%xst(2),'Ty')
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
        !write(*,*) 'test, bc-T', dm%fbcy_const(n, 5), dm%fbcy_ftp(4, n, 4)%t
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

  subroutine enforce_var_from_const(dm, ux, uy, uz, fbc0)
    use udf_type_mod
    use parameters_constant_mod
    use print_msg_mod
    implicit none 
    type(t_domain), intent(in) :: dm
    real(WP), intent(in), optional :: fbc0(6)
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent (inout) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent (inout) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent (inout) :: uz

    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil
    real(WP) :: fbc(6)

    if(.not. present(fbc0)) then
      fbc = ZERO
    else 
      fbc = fbc0
    end if

    ! -mx_rhs-
    if(dm%ibcx_nominal(1, 1) == IBC_DIRICHLET) ux(1,              :, :) = fbc(1)
    if(dm%ibcx_nominal(2, 1) == IBC_DIRICHLET) ux(dm%dpcc%xsz(1), :, :) = fbc(2)
    !-my_rhs-
    if(dm%ibcy_nominal(1, 2) == IBC_DIRICHLET .or. &
       dm%ibcy_nominal(2, 2) == IBC_DIRICHLET) then
      call transpose_x_to_y(uy, acpc_ypencil, dm%dcpc)
      if(dm%ibcy_nominal(1, 2) == IBC_DIRICHLET) acpc_ypencil(:, 1,              :) = fbc(3)
      if(dm%ibcy_nominal(2, 2) == IBC_DIRICHLET) acpc_ypencil(:, dm%dcpc%ysz(2), :) = fbc(4)
      call transpose_y_to_x(acpc_ypencil, uy, dm%dcpc)
    end if
    !-mz_rhs-
    if(dm%ibcz_nominal(1, 3)  == IBC_DIRICHLET .or. &
       dm%ibcz_nominal(2, 3)  == IBC_DIRICHLET) then
      call transpose_x_to_y(uz, accp_ypencil, dm%dccp)
      call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
      if(dm%ibcz_nominal(1, 3) == IBC_DIRICHLET) accp_zpencil(:, :, 1             ) = fbc(5)
      if(dm%ibcz_nominal(2, 3) == IBC_DIRICHLET) accp_zpencil(:, :, dm%dccp%zsz(3)) = fbc(6)
      call transpose_z_to_y(accp_zpencil, accp_ypencil, dm%dccp)
      call transpose_y_to_x(accp_ypencil, uz, dm%dccp)
    end if

    return
  end subroutine


end module