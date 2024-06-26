module boundary_conditions_mod
  use parameters_constant_mod
  character(12) :: filename(5)
  
  private :: map_bc_1d_uprofile
  private :: map_1d_profile_to_case
  private :: reassign_bc_type
  public  :: update_symmetric_ibc
  public  :: configure_bc_type
  public  :: configure_bc_vars_flow
  public  :: configure_bc_vars_thermo
  public  :: update_bc_interface_flow
  public  :: update_bc_interface_thermo
  public  :: update_flow_bc_1dm_halo
  
  public :: buildup_symmetric_for_eqs

  public :: get_dirichlet_geo_bcx
  public :: get_dirichlet_geo_bcy
  public :: get_dirichlet_geo_bcz

  !public  :: apply_bc_const
  !public  :: apply_convective_outlet

  
contains

  subroutine reassign_bc_type(bc_nominal)
    implicit none 
    integer, intent(inout) :: bc_nominal(2, NBC)
    integer :: n, m

    do n = 1, 2
      do m = 1, NBC
        if (bc_nominal(n, m) == IBC_PROFILE1D)   then
          bc_nominal(n, m) = IBC_DIRICHLET
        else if (bc_nominal(n, m) == IBC_TURBGEN  .or. &
                 bc_nominal(n, m) == IBC_DATABASE )   then
          if(m == 5) then
            bc_nominal(n, m) = IBC_DIRICHLET ! for temperature, default is no incoming thermal flow, check
          else 
            bc_nominal(n, m) = IBC_INTERIOR  ! for u, v, w, p
          end if
        else if (bc_nominal(n, m) == IBC_CONVECTIVE)   then ! check for convetive outlet
          bc_nominal(n, m) = IBC_INTRPL
        else
          bc_nominal(n, m) = bc_nominal(n, m)   
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
! to exclude non-resonable input
!----------------------------------------------------------------------------------------------------------
    do m = 1, NBC
      if(dm%ibcx_nominal(2, m) == IBC_PROFILE1D) call Print_error_msg(" This BC IBC_PROFILE1D is not supported.")
      do n = 1, 2
        if(dm%ibcx_nominal(n, m) >  IBC_OTHERS   ) dm%ibcx_nominal(n, m) = IBC_INTRPL
        if(dm%ibcy_nominal(n, m) >  IBC_OTHERS   ) dm%ibcy_nominal(n, m) = IBC_INTRPL
        if(dm%ibcz_nominal(n, m) >  IBC_OTHERS   ) dm%ibcz_nominal(n, m) = IBC_INTRPL
        if(dm%ibcy_nominal(n, m) == IBC_PROFILE1D) call Print_error_msg(" This BC IBC_PROFILE1D is not supported.")
        if(dm%ibcz_nominal(n, m) == IBC_PROFILE1D) call Print_error_msg(" This BC IBC_PROFILE1D is not supported.")
      end do
    end do
!----------------------------------------------------------------------------------------------------------
! to check velocity symmetric and asymmetric
!----------------------------------------------------------------------------------------------------------
    do n = 1, 2
      if(dm%ibcx_nominal(n, 1) == IBC_SYMMETRIC) &
         dm%ibcx_nominal(n, 1) = IBC_ASYMMETRIC
      if(dm%ibcy_nominal(n, 2) == IBC_SYMMETRIC) &
         dm%ibcy_nominal(n, 2) = IBC_ASYMMETRIC
      if(dm%ibcz_nominal(n, 3) == IBC_SYMMETRIC) &
         dm%ibcz_nominal(n, 3) = IBC_ASYMMETRIC
    end do
!----------------------------------------------------------------------------------------------------------
! to set up real bc for calculation from given nominal b.c.
!----------------------------------------------------------------------------------------------------------
    call reassign_bc_type(dm%ibcx_nominal)
    call reassign_bc_type(dm%ibcy_nominal)
    call reassign_bc_type(dm%ibcz_nominal)
!----------------------------------------------------------------------------------------------------------
! allocate bc to variables
!----------------------------------------------------------------------------------------------------------
    do n = 1, 2
      dm%ibcx_qx(n) = dm%ibcx_nominal(n, 1)
      dm%ibcx_qy(n) = dm%ibcx_nominal(n, 2)
      dm%ibcx_qz(n) = dm%ibcx_nominal(n, 3)
      dm%ibcx_pr(n) = dm%ibcx_nominal(n, 4)
      dm%ibcx_Th(n) = dm%ibcx_nominal(n, 5)

      dm%ibcy_qx(n) = dm%ibcy_nominal(n, 1)
      dm%ibcy_qy(n) = dm%ibcy_nominal(n, 2)
      dm%ibcy_qz(n) = dm%ibcy_nominal(n, 3)
      dm%ibcy_pr(n) = dm%ibcy_nominal(n, 4)
      dm%ibcy_Th(n) = dm%ibcy_nominal(n, 5)

      dm%ibcz_qx(n) = dm%ibcz_nominal(n, 1)
      dm%ibcz_qy(n) = dm%ibcz_nominal(n, 2)
      dm%ibcz_qz(n) = dm%ibcz_nominal(n, 3)
      dm%ibcz_pr(n) = dm%ibcz_nominal(n, 4)
      dm%ibcz_Th(n) = dm%ibcz_nominal(n, 5)
    end do 
!----------------------------------------------------------------------------------------------------------
! correct bc for future calculation. only here uses the IBC_CCC, etc...
!----------------------------------------------------------------------------------------------------------
    ! do n = 1, 2
    !   ! - x - 
    !   if(dm%ibcx_nominal(n, 1) == IBC_DIRICHLET) then ! ux at x-dir
    !     dm%ibcx_qx(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcx_qx(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcx_qx(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcx_qx(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !   end if 
    !   if(dm%ibcx_nominal(n, 2) == IBC_DIRICHLET) then ! uy at x-dir
    !     dm%ibcx_qy(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcx_qy(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcx_qy(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcx_qy(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !   end if 
    !   if(dm%ibcx_nominal(n, 3) == IBC_DIRICHLET) then ! uz at x-dir
    !     dm%ibcx_qz(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcx_qz(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcx_qz(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcx_qz(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !   end if 
    !   ! - y - 
    !   if(dm%ibcy_nominal(n, 1) == IBC_DIRICHLET) then
    !     dm%ibcy_qx(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcy_qx(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcy_qx(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcy_qx(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !   end if 
    !   if(dm%ibcy_nominal(n, 2) == IBC_DIRICHLET) then
    !     dm%ibcy_qy(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcy_qy(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcy_qy(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcy_qy(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !   end if 
    !   if(dm%ibcy_nominal(n, 3) == IBC_DIRICHLET) then
    !     dm%ibcy_qz(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcy_qz(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcy_qz(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcy_qz(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !   end if 
    !   ! - z - 
    !   if(dm%ibcz_nominal(n, 1) == IBC_DIRICHLET) then
    !     dm%ibcz_qx(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcz_qx(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcz_qx(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcz_qx(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !   end if 
    !   if(dm%ibcz_nominal(n, 2) == IBC_DIRICHLET) then
    !     dm%ibcz_qy(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcz_qy(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcz_qy(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcz_qy(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !   end if 
    !   if(dm%ibcz_nominal(n, 3)== IBC_DIRICHLET) then
    !     dm%ibcz_qz(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcz_qz(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcz_qz(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !     dm%ibcz_qz(n) = IBC_ASYMMETRIC ! default, non-slip wall, check
    !   end if 

    ! end do

    return
  end subroutine 
!==========================================================================================================
!==========================================================================================================
  subroutine map_1d_profile_to_case(nin, yin, uin, nout, ycase, ucase)
    use cubic_spline_interpolation
    use precision_mod
    implicit none
    integer, intent(in) :: nin
    real(WP), dimension(nin), intent(in) :: yin
    real(WP), dimension(nin), intent(in) :: uin
    integer, intent(in) :: nout
    real(WP), dimension(nout), intent(in)  :: ycase
    real(WP), dimension(nout), intent(out) :: ucase

    integer :: i
    real(WP), allocatable :: cs_b(:), cs_c(:), cs_d(:)

    allocate(cs_b(nin))
    allocate(cs_c(nin))
    allocate(cs_d(nin))

    call cubic_spline (nin, yin, uin, cs_b, cs_c, cs_d)

    do i = 1, nout
      ucase(i) = spline_interpolation(nin, uin, cs_b, cs_c, cs_d, ycase(i))
    end do 

    deallocate(cs_b)
    deallocate(cs_c)
    deallocate(cs_d)

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
  subroutine configure_bc_vars_flow(dm)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout)  :: dm

    real(WP) :: var1y(1:dm%np(2))
    integer  :: m, n, k, ny
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

    if(dm%is_thermo) then
      allocate( dm%fbcx_gx(             4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)) )! default x pencil
      allocate( dm%fbcy_gx(dm%dpcc%ysz(1),              4, dm%dpcc%ysz(3)) )! default y pencil
      allocate( dm%fbcz_gx(dm%dpcc%zsz(1), dm%dpcc%zsz(2),              4) )! default z pencil

      allocate( dm%fbcx_gy(             4, dm%dcpc%xsz(2), dm%dcpc%xsz(3)) )! default x pencil
      allocate( dm%fbcy_gy(dm%dcpc%ysz(1),              4, dm%dcpc%ysz(3)) )! default y pencil
      allocate( dm%fbcz_gy(dm%dcpc%zsz(1), dm%dcpc%zsz(2),              4) )! default z pencil

      allocate( dm%fbcx_gz(             4, dm%dccp%xsz(2), dm%dccp%xsz(3)) )! default x pencil
      allocate( dm%fbcy_gz(dm%dccp%ysz(1),              4, dm%dccp%ysz(3)) )! default y pencil
      allocate( dm%fbcz_gz(dm%dccp%zsz(1), dm%dccp%zsz(2),              4) )! default z pencil
    end if
    ! -3-1-||||-2-4
    !do m = 1, NBC ! u, v, w, p, T(dim)

! ux - 
    do n = 1, 2
      dm%fbcx_qx(n, :, :) = dm%fbcx_const(n, 1)
      dm%fbcy_qx(:, n, :) = dm%fbcy_const(n, 1)
      dm%fbcz_qx(:, :, n) = dm%fbcz_const(n, 1)
    end do
    do n = 3, 4
      dm%fbcx_qx(n, :, :) = dm%fbcx_const(n - 2, 1)
      dm%fbcy_qx(:, n, :) = dm%fbcy_const(n - 2, 1)
      dm%fbcz_qx(:, :, n) = dm%fbcz_const(n - 2, 1)
    end do
! uy- 
    do n = 1, 2
      dm%fbcx_qy(n, :, :) = dm%fbcx_const(n, 2)
      dm%fbcy_qy(:, n, :) = dm%fbcy_const(n, 2)
      dm%fbcz_qy(:, :, n) = dm%fbcz_const(n, 2)
    end do
    do n = 3, 4
      dm%fbcx_qy(n, :, :) = dm%fbcx_const(n - 2, 2)
      dm%fbcy_qy(:, n, :) = dm%fbcy_const(n - 2, 2)
      dm%fbcz_qy(:, :, n) = dm%fbcz_const(n - 2, 2)
    end do
! uz- 
    do n = 1, 2
      dm%fbcx_qz(n, :, :) = dm%fbcx_const(n, 3)
      dm%fbcy_qz(:, n, :) = dm%fbcy_const(n, 3)
      dm%fbcz_qz(:, :, n) = dm%fbcz_const(n, 3)
    end do
    do n = 3, 4
      dm%fbcx_qz(n, :, :) = dm%fbcx_const(n - 2, 3)
      dm%fbcy_qz(:, n, :) = dm%fbcy_const(n - 2, 3)
      dm%fbcz_qz(:, :, n) = dm%fbcz_const(n - 2, 3)
    end do
! pr- 
    do n = 1, 2
      dm%fbcx_pr(n, :, :) = dm%fbcx_const(n, 4)
      dm%fbcy_pr(:, n, :) = dm%fbcy_const(n, 4)
      dm%fbcz_pr(:, :, n) = dm%fbcz_const(n, 4)
    end do
    do n = 3, 4
      dm%fbcx_pr(n, :, :) = dm%fbcx_const(n - 2, 4)
      dm%fbcy_pr(:, n, :) = dm%fbcy_const(n - 2, 4)
      dm%fbcz_pr(:, :, n) = dm%fbcz_const(n - 2, 4)
    end do

!----------------------------------------------------------------------------------------------------------
! to build up bc for var(x_const, y, z)
!----------------------------------------------------------------------------------------------------------
    ! ! ux - 
    ! if(dm%ibcx_nominal(1, 1) == IBC_PROFILE1D) then
    !   filename = 'pf1d_u1y.dat' !(undim)
    !   call map_bc_1d_uprofile( filename, dm%nc(2), dm%yc, var1y(1:dm%nc(2)) )
    !   call transpose_y_to_x
    !   dm%fbcx_qx(1, :, :) = 
    ! end if
    ! filename(2) = 'pf1d_v1y.dat' !(undim)
    ! filename(3) = 'pf1d_w1y.dat' !(undim)
    ! filename(4) = 'pf1d_p1y.dat' !(undim)
    ! filename(5) = 'pf1d_T1y.dat' !(dim  )
    ! do m = 1, NBC
    !   if(dm%ibcx_nominal(1, m) == IBC_PROFILE1D) then

    !     if(m /= 2) then
    !       ny = dm%nc(2)
    !       call map_bc_1d_uprofile( filename(m), ny, dm%yc, var1y(1:ny) )
    !     else
    !       ny = dm%np(2)
    !       call map_bc_1d_uprofile( filename(m), ny, dm%yp, var1y(1:ny) )
    !     end if

    !     do k = 1, size(dm%fbcx_var, 3) 
    !       do j = 1, 
    !       dm%fbcx_var(1, 1:ny, k, m) = var1y(1:ny) ! wrong
    !     end do

    !   end if
    ! end do

    return
  end subroutine 
!==========================================================================================================
!==========================================================================================================
subroutine configure_bc_vars_thermo(dm)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm

    if( .not. dm%is_thermo) return

    allocate( dm%ftpbcx_var(             4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)) )! default x pencil
    allocate( dm%ftpbcy_var(dm%dcpc%ysz(1),              4, dm%dcpc%ysz(3)) )! default y pencil
    allocate( dm%ftpbcz_var(dm%dccp%zsz(1), dm%dccp%zsz(2),             4)  )! default z pencil

    return
  end subroutine 

!==========================================================================================================
!==========================================================================================================
  subroutine update_flow_bc_1dm_halo(dm, fl) ! similar to asymmetric
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    type(t_flow), intent(inout)    :: fl

    real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: apcc_ypencil
    real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: apcc_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil !
    real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: acpc_zpencil !
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil

    !------------ux at x-bc-----------
    if(dm%ibcx_qx(1) == IBC_DIRICHLET ) then
      dm%fbcx_qx(1, :, :) = dm%fbcx_const(1, 1)
      dm%fbcx_qx(3, :, :) = dm%fbcx_const(1, 1)
      fl%qx(1, :, :) = dm%fbcx_const(1, 1)
    end if
    if(dm%ibcx_qx(2) == IBC_DIRICHLET ) then
      dm%fbcx_qx(2, :, :) = dm%fbcx_const(2, 1)
      dm%fbcx_qx(4, :, :) = dm%fbcx_const(2, 1)
      fl%qx(dm%dpcc%xsz(1), :, :) = dm%fbcx_const(2, 1)
    end if
    !------------ux at y-bc-----------
    if(dm%ibcy_qx(1) == IBC_INTERIOR .or. &
       dm%ibcy_qx(2) == IBC_INTERIOR) then
      call transpose_x_to_y(fl%qx, apcc_ypencil, dm%dpcc)
    end if
    ! ux at y-bc-y1
    if(dm%ibcy_qx(1) == IBC_INTERIOR) then
      if( dm%dpcc%yst(2)    == 1) dm%fbcy_qx(:, 1, :) = - apcc_ypencil(:, 1, :) 
      if((dm%dpcc%yst(2)+1) == 2) dm%fbcy_qx(:, 3, :) = - apcc_ypencil(:, 2, :)
    end if
    ! ux at y-bc-yn
    if(dm%ibcy_qx(2) == IBC_INTERIOR) then
      if( dm%dpcc%yen(2)    ==  dm%nc(2)   ) dm%fbcy_qx(:, 2, :) = - apcc_ypencil(:, dm%nc(2),     :)
      if((dm%dpcc%yen(2)-1) == (dm%nc(2)-1)) dm%fbcy_qx(:, 4, :) = - apcc_ypencil(:, dm%nc(2) - 1, :)
    end if
    !-----------ux at z-bc-----------
    if(dm%ibcz_qx(1) == IBC_INTERIOR .or. &
       dm%ibcz_qx(2) == IBC_INTERIOR) then
      call transpose_x_to_y(fl%qx,        apcc_ypencil, dm%dpcc)
      call transpose_y_to_z(apcc_ypencil, apcc_zpencil, dm%dpcc)
    end if
    ! ux at z-bc
    if(dm%ibcz_qx(1) == IBC_INTERIOR) then
      if( dm%dpcc%zst(3)    == 1) dm%fbcz_qx(:, :, 1) = - apcc_zpencil(:, :, 1) 
      if((dm%dpcc%zst(3)+1) == 2) dm%fbcz_qx(:, :, 3) = - apcc_zpencil(:, :, 2)
    end if
    if(dm%ibcz_qx(2) == IBC_INTERIOR) then
      if( dm%dpcc%zen(3)    ==  dm%nc(3)   ) dm%fbcz_qx(:, :, 2) = - apcc_zpencil(:, :, dm%nc(3)    )
      if((dm%dpcc%zen(3)-1) == (dm%nc(3)-1)) dm%fbcz_qx(:, :, 4) = - apcc_zpencil(:, :, dm%nc(3) - 1)
    end if

    !------------uy at y-bc-----------
    if(dm%ibcy_qy(1) == IBC_DIRICHLET ) then
      dm%fbcy_qy(:, 1, :) = dm%fbcy_const(1, 2)
      dm%fbcy_qy(:, 3, :) = dm%fbcy_const(1, 2)
      call transpose_x_to_y(fl%qy, acpc_ypencil, dm%dcpc)
      acpc_ypencil(:, 1, :) = dm%fbcy_const(1, 2)
      call transpose_y_to_x(acpc_ypencil, fl%qy, dm%dcpc)
    end if
    if(dm%ibcy_qy(2) == IBC_DIRICHLET ) then
      dm%fbcy_qy(:, 2, :) = dm%fbcy_const(2, 2)
      dm%fbcy_qy(:, 4, :) = dm%fbcy_const(2, 2)
      call transpose_x_to_y(fl%qy, acpc_ypencil, dm%dcpc)
      acpc_ypencil(:, dm%dcpc%ysz(2), :) = dm%fbcy_const(2, 2)
      call transpose_y_to_x(acpc_ypencil, fl%qy, dm%dcpc)
    end if
    !-----------uy at x-bc-----------
    if(dm%ibcx_qy(1) == IBC_INTERIOR) then
      if( dm%dcpc%xst(1)    == 1) dm%fbcx_qy(1, :, :) = - fl%qy(1, :, :) 
      if((dm%dcpc%xst(1)+1) == 2) dm%fbcx_qy(3, :, :) = - fl%qy(2, :, :)
    end if
    if(dm%ibcx_qy(2) == IBC_INTERIOR) then
      if( dm%dcpc%xen(1)    ==  dm%nc(1)   ) dm%fbcx_qy(2, :, :) = - fl%qy(dm%nc(1),     :, :)
      if((dm%dcpc%xen(1)-1) == (dm%nc(1)-1)) dm%fbcx_qy(4, :, :) = - fl%qy(dm%nc(1) - 1, :, :)
    end if
    !-----------uy at z-bc-----------
    if(dm%ibcz_qy(1) == IBC_INTERIOR .or. &
       dm%ibcz_qy(2) == IBC_INTERIOR) then
      call transpose_x_to_y(fl%qy,         acpc_ypencil, dm%dcpc)
      call transpose_y_to_z(acpc_ypencil, acpc_zpencil, dm%dcpc)
    end if
    if(dm%ibcz_qy(1) == IBC_INTERIOR) then
      if( dm%dcpc%zst(3)    == 1) dm%fbcz_qy(:, :, 1) = - acpc_zpencil(:, :, 1) 
      if((dm%dcpc%zst(3)+1) == 2) dm%fbcz_qy(:, :, 3) = - acpc_zpencil(:, :, 2)
    end if
    if(dm%ibcz_qy(2) == IBC_INTERIOR) then
      if( dm%dpcc%zen(3)    ==  dm%nc(3)   ) dm%fbcz_qy(:, :, 2) = - acpc_zpencil(:, :, dm%nc(3)    )
      if((dm%dpcc%zen(3)-1) == (dm%nc(3)-1)) dm%fbcz_qy(:, :, 4) = - acpc_zpencil(:, :, dm%nc(3) - 1)
    end if

    !------------uz at z-bc-----------
    ! check, to add transport
    if(dm%ibcz_qz(1) == IBC_DIRICHLET ) then
      dm%fbcz_qz(:, 1, :) = dm%fbcz_const(1, 3)
      dm%fbcz_qz(:, 3, :) = dm%fbcz_const(1, 3)
    end if
    if(dm%ibcz_qz(2) == IBC_DIRICHLET ) then
      dm%fbcz_qz(:, 2, :) = dm%fbcz_const(2, 3)
      dm%fbcz_qz(:, 4, :) = dm%fbcz_const(2, 3)
    end if
    !-----------uz at x-bc-----------
    if(dm%ibcx_qz(1) == IBC_INTERIOR) then
      if( dm%dccp%xst(1)    == 1) dm%fbcx_qz(1, :, :) = - fl%qz(1, :, :) 
      if((dm%dccp%xst(1)+1) == 2) dm%fbcx_qz(3, :, :) = - fl%qz(2, :, :)
    end if
    if(dm%ibcx_qz(2) == IBC_INTERIOR) then
      if( dm%dccp%xen(1)    ==  dm%nc(1)   ) dm%fbcx_qz(2, :, :) = - fl%qz(dm%nc(1),     :, :)
      if((dm%dccp%xen(1)-1) == (dm%nc(1)-1)) dm%fbcx_qz(4, :, :) = - fl%qz(dm%nc(1) - 1, :, :)
    end if
    !-----------uz at y-bc-----------
    if(dm%ibcy_qz(1) == IBC_INTERIOR .or. &
       dm%ibcy_qz(2) == IBC_INTERIOR) then
      call transpose_x_to_y(fl%qz, accp_ypencil, dm%dpcc)
    end if
    if(dm%ibcy_qz(1) == IBC_INTERIOR) then
      if( dm%dccp%yst(2)    == 1) dm%fbcy_qz(:, 1, :) = - accp_ypencil(:, 1, :) 
      if((dm%dccp%yst(2)+1) == 2) dm%fbcy_qz(:, 3, :) = - accp_ypencil(:, 2, :)
    end if
    if(dm%ibcy_qz(2) == IBC_INTERIOR) then
      if( dm%dccp%yen(2)    ==  dm%nc(2)   ) dm%fbcy_qz(:, 2, :) = - accp_ypencil(:, dm%nc(2),     :)
      if((dm%dccp%yen(2)-1) == (dm%nc(2)-1)) dm%fbcy_qz(:, 4, :) = - accp_ypencil(:, dm%nc(2) - 1, :)
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
!----------------------------------------------------------------------------------------------------------
!   all in x-pencil
!   no overlap of values
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
!   bc in x - direction
!----------------------------------------------------------------------------------------------------------
    ! ux, dm0-dm1, np
    if(dm1%ibcx_qx(1) == IBC_INTERIOR) then
      dm1%fbcx_qx(1, :, :) = fl0%qx(dm0%np_geo(1) - 1, :, :)
      dm1%fbcx_qx(3, :, :) = fl0%qx(dm0%np_geo(1) - 2, :, :)
    end if
    if(dm0%ibcx_qx(2) == IBC_INTERIOR) then
      dm0%fbcx_qx(2, :, :) = fl1%qx(2, :, :)
      dm0%fbcx_qx(4, :, :) = fl1%qx(3, :, :)
    end if

    ! uy, dm0-dm1, nc
    if(dm1%ibcx_qy(1) == IBC_INTERIOR) then
      dm1%fbcx_qy(1, :, :) = fl0%qy(dm0%nc(1),     :, :)
      dm1%fbcx_qy(3, :, :) = fl0%qy(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx_qy(2) == IBC_INTERIOR) then
      dm0%fbcx_qy(2, :, :) = fl1%qy(1, :, :)
      dm0%fbcx_qy(4, :, :) = fl1%qy(2, :, :)
    end if

    ! uz, dm0-dm1, nc
    if(dm1%ibcx_qz(1) == IBC_INTERIOR) then
      dm1%fbcx_qz(1, :, :) = fl0%qz(dm0%nc(1),     :, :)
      dm1%fbcx_qz(3, :, :) = fl0%qz(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx_qz(2) == IBC_INTERIOR) then
      dm0%fbcx_qz(2, :, :) = fl1%qz(1, :, :)
      dm0%fbcx_qz(4, :, :) = fl1%qz(2, :, :)
    end if

    ! p, dm0-dm1, nc
    if(dm1%ibcx_pr(1) == IBC_INTERIOR) then
      dm1%fbcx_pr(1, :, :) = fl0%pres(dm0%nc(1),     :, :)
      dm1%fbcx_pr(3, :, :) = fl0%pres(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx_pr(2) == IBC_INTERIOR) then
      dm0%fbcx_pr(2, :, :) = fl1%pres(1, :, :)
      dm0%fbcx_pr(4, :, :) = fl1%pres(2, :, :)
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in y - direction
!----------------------------------------------------------------------------------------------------------
    ! ux, dm0-dm1, nc
    if(dm1%ibcy_qx(1) == IBC_INTERIOR) then
      dm1%fbcy_qx(:, 1, :) = fl0%qx(:, dm0%nc(2),     :)
      dm1%fbcy_qx(:, 3, :) = fl0%qx(:, dm0%nc(2) - 1, :)
    end if
    if(dm0%ibcy_qx(2) == IBC_INTERIOR) then
      dm0%fbcy_qx(:, 2, :) = fl1%qx(:, 1, :)
      dm0%fbcy_qx(:, 4, :) = fl1%qx(:, 2, :)
    end if

    ! uy, dm0-dm1, np
    if(dm1%ibcy_qy(1) == IBC_INTERIOR) then
      dm1%fbcy_qy(:, 1, :) = fl0%qy(:, dm0%np_geo(2) - 1, :)
      dm1%fbcy_qy(:, 3, :) = fl0%qy(:, dm0%np_geo(2) - 2, :)
    end if
    if(dm0%ibcy_qy(2) == IBC_INTERIOR) then
      dm0%fbcy_qy(:, 2, :) = fl1%qy(:, 2, :)
      dm0%fbcy_qy(:, 4, :) = fl1%qy(:, 3, :)
    end if

    ! uz, dm0-dm1, nc
    if(dm1%ibcy_qz(1) == IBC_INTERIOR) then
      dm1%fbcy_qz(:, 1, :) = fl0%qz(:, dm0%nc(2),     :)
      dm1%fbcy_qz(:, 3, :) = fl0%qz(:, dm0%nc(2) - 1, :)
    end if
    if(dm0%ibcy_qz(2) == IBC_INTERIOR) then
      dm0%fbcy_qz(:, 2, :) = fl1%qz(:, 1, :)
      dm0%fbcy_qz(:, 4, :) = fl1%qz(:, 2, :)
    end if

    ! p, dm0-dm1, nc
    if(dm1%ibcy_pr(1) == IBC_INTERIOR) then
      dm1%fbcy_pr(:, 1, :) = fl0%pres(:, dm0%nc(2),     :)
      dm1%fbcy_pr(:, 3, :) = fl0%pres(:, dm0%nc(2) - 1, :)
    end if
    if(dm0%ibcy_pr(2) == IBC_INTERIOR) then
      dm0%fbcy_pr(:, 2, :) = fl1%pres(:, 1, :)
      dm0%fbcy_pr(:, 4, :) = fl1%pres(:, 2, :)
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in z - direction
!----------------------------------------------------------------------------------------------------------
    ! ux, dm0-dm1, nc
    if(dm1%ibcz_qx(1) == IBC_INTERIOR) then
      dm1%fbcz_qx(:, :, 1) = fl0%qx(:, :, dm0%nc(2)    )
      dm1%fbcz_qx(:, :, 3) = fl0%qx(:, :, dm0%nc(2) - 1)
    end if
    if(dm0%ibcz_qx(2) == IBC_INTERIOR) then
      dm0%fbcz_qx(:, :, 2) = fl1%qx(:, :, 1)
      dm0%fbcz_qx(:, :, 4) = fl1%qx(:, :, 2)
    end if

    ! uy, dm0-dm1, nc
    if(dm1%ibcz_qy(1) == IBC_INTERIOR) then
      dm1%fbcz_qy(:, :, 1) = fl0%qy(:, :, dm0%nc(2)    )
      dm1%fbcz_qy(:, :, 3) = fl0%qy(:, :, dm0%nc(2) - 1)
    end if
    if(dm0%ibcz_qy(2) == IBC_INTERIOR) then
      dm0%fbcz_qy(:, :, 2) = fl1%qy(:, :, 1)
      dm0%fbcz_qy(:, :, 4) = fl1%qy(:, :, 2)
    end if

    ! uz, dm0-dm1, np
    if(dm1%ibcz_qz(1) == IBC_INTERIOR) then
      dm1%fbcz_qz(:, :, 1) = fl0%qz(:, :, dm0%np_geo(2) - 1)
      dm1%fbcz_qz(:, :, 3) = fl0%qz(:, :, dm0%np_geo(2) - 2)
    end if
    if(dm0%ibcz_qz(2) == IBC_INTERIOR) then
      dm0%fbcz_qz(:, :, 2) = fl1%qz(:, :, 2)
      dm0%fbcz_qz(:, :, 4) = fl1%qz(:, :, 3)
    end if

    ! p, dm0-dm1, nc
    if(dm1%ibcz_pr(1) == IBC_INTERIOR) then
      dm1%fbcz_pr(:, :, 1) = fl0%pres(:, :, dm0%nc(2)    )
      dm1%fbcz_pr(:, :, 3) = fl0%pres(:, :, dm0%nc(2) - 1)
    end if
    if(dm0%ibcz_pr(2) == IBC_INTERIOR) then
      dm0%fbcz_pr(:, :, 2) = fl1%pres(:, :, 1)
      dm0%fbcz_pr(:, :, 4) = fl1%pres(:, :, 2)
    end if
    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine update_bc_interface_thermo(dm0, fl0, tm0, dm1, fl1, tm1)
    use parameters_constant_mod
    use thermo_info_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm0, dm1
    type(t_flow), intent(in)      :: fl0, fl1
    type(t_thermo), intent(in)    :: tm0, tm1
    
    integer :: i, j, k
!----------------------------------------------------------------------------------------------------------
!   all in x-pencil
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
!   bc in x - direction
!----------------------------------------------------------------------------------------------------------
    ! gx, dm0-dm1, np
    if(dm1%ibcx_qx(1) == IBC_INTERIOR) then
      dm1%fbcx_gx(1, :, :) = fl0%gx(dm0%np_geo(1) - 1,     :, :)
      dm1%fbcx_gx(3, :, :) = fl0%gx(dm0%np_geo(1) - 2, :, :)
    end if
    if(dm0%ibcx_qx(2) == IBC_INTERIOR) then
      dm0%fbcx_gx(2, :, :) = fl1%gx(2, :, :)
      dm0%fbcx_gx(4, :, :) = fl1%gx(3, :, :)
    end if

    ! gy, dm0-dm1, nc
    if(dm1%ibcx_qy(1) == IBC_INTERIOR) then
      dm1%fbcx_gy(1, :, :) = fl0%gy(dm0%nc(1),     :, :)
      dm1%fbcx_gy(3, :, :) = fl0%gy(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx_qy(2) == IBC_INTERIOR) then
      dm0%fbcx_gy(2, :, :) = fl1%gy(1, :, :)
      dm0%fbcx_gy(4, :, :) = fl1%gy(2, :, :)
    end if

    ! gz, dm0-dm1, nc
    if(dm1%ibcx_qz(1) == IBC_INTERIOR) then
      dm1%fbcx_gz(1, :, :) = fl0%gz(dm0%nc(1),     :, :)
      dm1%fbcx_gz(3, :, :) = fl0%gz(dm0%nc(1) - 1, :, :)
    end if
    if(dm0%ibcx_qz(2) == IBC_INTERIOR) then
      dm0%fbcx_gz(2, :, :) = fl1%gz(1, :, :)
      dm0%fbcx_gz(4, :, :) = fl1%gz(2, :, :)
    end if

    ! thermal field, dm0-dm1
    if(dm1%ibcx_Th(1) == IBC_INTERIOR) then
      dm1%ftpbcx_var(1, :, :)%t = tm0%tTemp(dm0%nc(1),     :, :)
      dm1%ftpbcx_var(3, :, :)%t = tm0%tTemp(dm0%nc(1) - 1, :, :)
      do k = 1, dm1%np(3)
        do j = 1, dm1%np(2)
          call ftp_refresh_thermal_properties_from_T_undim(dm1%ftpbcx_var(1, j, k))
          call ftp_refresh_thermal_properties_from_T_undim(dm1%ftpbcx_var(3, j, k))
        end do
      end do
    end if
    if(dm0%ibcx_Th(2) == IBC_INTERIOR) then
      dm0%ftpbcx_var(2, :, :)%t = tm1%tTemp(1, :, :)
      dm0%ftpbcx_var(4, :, :)%t = tm1%tTemp(2, :, :)
      do k = 1, dm0%np(3)
        do j = 1, dm0%np(2)
          call ftp_refresh_thermal_properties_from_T_undim(dm0%ftpbcx_var(2, j, k))
          call ftp_refresh_thermal_properties_from_T_undim(dm0%ftpbcx_var(4, j, k))
        end do
      end do
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in y - direction
!----------------------------------------------------------------------------------------------------------
    ! gx, dm0-dm1, nc
    if(dm1%ibcy_qx(1) == IBC_INTERIOR) then
      dm1%fbcy_gx(:, 1, :) = fl0%gx(:, dm0%nc(1),     :)
      dm1%fbcy_gx(:, 3, :) = fl0%gx(:, dm0%nc(1) - 1, :)
    end if
    if(dm0%ibcy_qx(2) == IBC_INTERIOR) then
      dm0%fbcy_gx(:, 2, :) = fl1%gx(:, 1, :)
      dm0%fbcy_gx(:, 4, :) = fl1%gx(:, 2, :)
    end if

    ! gy, dm0-dm1, nc
    if(dm1%ibcy_qy(1) == IBC_INTERIOR) then
      dm1%fbcy_gy(:, 1, :) = fl0%gy(:, dm0%np_geo(1) - 1, :)
      dm1%fbcy_gy(:, 3, :) = fl0%gy(:, dm0%np_geo(1) - 2, :)
    end if
    if(dm0%ibcy_qy(2) == IBC_INTERIOR) then
      dm0%fbcy_gy(:, 2, :) = fl1%gy(:, 2, :)
      dm0%fbcy_gy(:, 4, :) = fl1%gy(:, 3, :)
    end if

    ! gz, dm0-dm1, nc
    if(dm1%ibcy_qz(1) == IBC_INTERIOR) then
      dm1%fbcy_gz(:, 1, :) = fl0%gz(:, dm0%nc(1),     :)
      dm1%fbcy_gz(:, 3, :) = fl0%gz(:, dm0%nc(1) - 1, :)
    end if
    if(dm0%ibcy_qz(2) == IBC_INTERIOR) then
      dm0%fbcy_gz(:, 2, :) = fl1%gz(:, 1, :)
      dm0%fbcy_gz(:, 4, :) = fl1%gz(:, 2, :)
    end if

    ! thermal field, dm0-dm1
    if(dm1%ibcy_Th(1) == IBC_INTERIOR) then
      dm1%ftpbcy_var(:, 1, :)%t = tm0%tTemp(:, dm0%nc(1),     :)
      dm1%ftpbcy_var(:, 3, :)%t = tm0%tTemp(:, dm0%nc(1) - 1, :)
      do k = 1, dm1%np(3)
        do i = 1, dm1%np(1)
          call ftp_refresh_thermal_properties_from_T_undim(dm1%ftpbcy_var(i, 1, k))
          call ftp_refresh_thermal_properties_from_T_undim(dm1%ftpbcy_var(i, 3, k))
        end do
      end do
    end if
    if(dm0%ibcy_Th(2) == IBC_INTERIOR) then
      dm0%ftpbcy_var(:, 2, :)%t = tm1%tTemp(:, 1, :)
      dm0%ftpbcy_var(:, 4, :)%t = tm1%tTemp(:, 2, :)
      do k = 1, dm0%np(3)
        do i = 1, dm0%np(1)
          call ftp_refresh_thermal_properties_from_T_undim(dm0%ftpbcy_var(i, 2, k))
          call ftp_refresh_thermal_properties_from_T_undim(dm0%ftpbcy_var(i, 4, k))
        end do
      end do
    end if
!----------------------------------------------------------------------------------------------------------
!   bc in z - direction
!----------------------------------------------------------------------------------------------------------
    ! gx, dm0-dm1, nc
    if(dm1%ibcz_qx(1) == IBC_INTERIOR) then
      dm1%fbcz_gx(:, :, 1) = fl0%gx(:, :, dm0%nc(1)    )
      dm1%fbcz_gx(:, :, 3) = fl0%gx(:, :, dm0%nc(1) - 1)
    end if
    if(dm0%ibcz_qx(2) == IBC_INTERIOR) then
      dm0%fbcz_gx(:, :, 2) = fl1%gx(:, :, 1)
      dm0%fbcz_gx(:, :, 4) = fl1%gx(:, :, 2)
    end if

    ! gy, dm0-dm1, nc
    if(dm1%ibcz_qy(1) == IBC_INTERIOR) then
      dm1%fbcz_gy(:, :, 1) = fl0%gy(:, :, dm0%nc(1)    )
      dm1%fbcz_gy(:, :, 3) = fl0%gy(:, :, dm0%nc(1) - 1)
    end if
    if(dm0%ibcz_qy(2) == IBC_INTERIOR) then
      dm0%fbcz_gy(:, :, 2) = fl1%gy(:, :, 1)
      dm0%fbcz_gy(:, :, 4) = fl1%gy(:, :, 2)
    end if

    ! gz, dm0-dm1, np
    if(dm1%ibcz_qz(1) == IBC_INTERIOR) then
      dm1%fbcz_gz(:, :, 1) = fl0%gz(:, :, dm0%np_geo(1) - 1)
      dm1%fbcz_gz(:, :, 3) = fl0%gz(:, :, dm0%np_geo(1) - 2)
    end if
    if(dm0%ibcz_qz(2) == IBC_INTERIOR) then
      dm0%fbcz_gz(:, :, 2) = fl1%gz(:, :, 2)
      dm0%fbcz_gz(:, :, 4) = fl1%gz(:, :, 3)
    end if

    ! thermal field, dm0-dm1
    if(dm1%ibcz_Th(1) == IBC_INTERIOR) then
      dm1%ftpbcz_var(:, :, 1)%t = tm0%tTemp(:, :, dm0%nc(1)    )
      dm1%ftpbcz_var(:, :, 3)%t = tm0%tTemp(:, :, dm0%nc(1) - 1)
      do j = 1, dm1%np(2)
        do i = 1, dm1%np(1)
          call ftp_refresh_thermal_properties_from_T_undim(dm1%ftpbcz_var(i, j, 1))
          call ftp_refresh_thermal_properties_from_T_undim(dm1%ftpbcz_var(i, j, 3))
        end do
      end do
    end if
    if(dm0%ibcz_Th(2) == IBC_INTERIOR) then
      dm0%ftpbcz_var(:, :, 2)%t = tm1%tTemp(:, :, 1)
      dm0%ftpbcz_var(:, :, 4)%t = tm1%tTemp(:, :, 2)
      do j = 1, dm0%np(2)
        do i = 1, dm0%np(1)
          call ftp_refresh_thermal_properties_from_T_undim(dm0%ftpbcz_var(i, j, 2))
          call ftp_refresh_thermal_properties_from_T_undim(dm0%ftpbcz_var(i, j, 4))
        end do
      end do
    end if
    return
  end subroutine
!==========================================================================================================
! to calculate boundary during calculation from primary boundary
  subroutine update_symmetric_ibc(ibc, mbc, jbc)
    use parameters_constant_mod
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
            if(nrank==0) write(*, *) "BCs for size ", i, " are ", ibc(i), jbc(i) 
            call Print_warning_msg("The two operational variables have different boundary conditions.")
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
  subroutine buildup_symmetric_for_eqs(dm)
    use parameters_constant_mod
    use udf_type_mod
    use wtformat_mod
    implicit none
    type(t_domain), intent( inout)   :: dm
    
    integer :: mbc(2, 3), mbc0(2, 3)
    integer :: bc(2)
!----------------------------------------------------------------------------------------------------------
!   x-mom
!----------------------------------------------------------------------------------------------------------
    call update_symmetric_ibc(dm%ibcx_qx, mbc, dm%ibcx_qx)
    mbcx_cov1(:) = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom x-convection is ", mbcx_cov1 

    call update_symmetric_ibc(dm%ibcy_qy, mbc, dm%ibcy_qx)
    mbcy_cov1(:) = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom y-convection is ", mbcy_cov1

    call update_symmetric_ibc(dm%ibcz_qz, mbc, dm%ibcz_qx)
    mbcz_cov1(:) = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom z-convection is ", mbcz_cov1

    call update_symmetric_ibc(dm%ibcx_qx, mbc)
    bc(:) = mbc(:, JBC_GRAD)
    call update_symmetric_ibc(dm%ibcx_Th, mbc, bc)
    mbcx_tau1 = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom x-viscous is ", mbcx_tau1

    call update_symmetric_ibc(dm%ibcy_qx, mbc)
    bc(:) = mbc(:, JBC_GRAD)
    call update_symmetric_ibc(dm%ibcy_Th, mbc, bc)
    call update_symmetric_ibc(dm%ibcy_Th, mbc0, dm%ibcy_qy)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcy_tau1 is wrong.")
    mbcy_tau1 = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom y-viscous is ", mbcy_tau1

    call update_symmetric_ibc(dm%ibcz_qx, mbc)
    bc(:) = mbc(:, JBC_GRAD)
    call update_symmetric_ibc(dm%ibcy_Th, mbc, bc)
    call update_symmetric_ibc(dm%ibcy_Th, mbc0, dm%ibcz_qz)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcy_tau1 is wrong.")
    mbcz_tau1 = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom z-viscous is ", mbcz_tau1
!----------------------------------------------------------------------------------------------------------
!   y-mom
!----------------------------------------------------------------------------------------------------------
    call update_symmetric_ibc(dm%ibcx_qx, mbc, dm%ibcx_qy)
    mbcx_cov2(:) = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom x-convection is ", mbcx_cov2

    call update_symmetric_ibc(dm%ibcy_qy, mbc, dm%ibcy_qy)
    mbcy_cov2(:) = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom y-convection is ", mbcy_cov2

    call update_symmetric_ibc(dm%ibcz_qz, mbc, dm%ibcz_qy)
    mbcz_cov2(:) = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom z-convection is ", mbcz_cov2

    if(dm%icoordinate == ICYLINDRICAL) then
      call update_symmetric_ibc(dm%ibcy_qz, mbc, dm%ibcy_qz)
      mbcr_cov2(:) = mbc(:, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom r-convection is ", mbcr_cov2
    end if

    call update_symmetric_ibc(dm%ibcy_qy, mbc)
    bc(:) = mbc(:, JBC_GRAD)
    call update_symmetric_ibc(dm%ibcy_Th, mbc, bc)
    mbcy_tau2 = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom y-viscous is ", mbcy_tau2

    call update_symmetric_ibc(dm%ibcx_qy, mbc)
    bc(:) = mbc(:, JBC_GRAD)
    call update_symmetric_ibc(dm%ibcx_Th, mbc, bc)
    call update_symmetric_ibc(dm%ibcx_Th, mbc0, dm%ibcx_qx)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcx_tau2 is wrong.")
    mbcx_tau2 = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom x-viscous is ", mbcx_tau2

    call update_symmetric_ibc(dm%ibcz_qy, mbc)
    bc(:) = mbc(:, JBC_GRAD)
    call update_symmetric_ibc(dm%ibcz_Th, mbc, bc)
    call update_symmetric_ibc(dm%ibcz_Th, mbc0, dm%ibcz_qz)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcz_tau2 is wrong.")
    mbcz_tau2 = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom z-viscous is ", mbcz_tau2

    if(dm%icoordinate == ICYLINDRICAL) then
      call update_symmetric_ibc(dm%ibcy_qz, mbc, dm%ibcy_Th)
      mbcr_tau2 = mbc(:, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom r-vicous is ", mbcr_tau2
    end if
!----------------------------------------------------------------------------------------------------------
!   z-mom
!----------------------------------------------------------------------------------------------------------
    call update_symmetric_ibc(dm%ibcx_qx, mbc, dm%ibcx_qz)
    mbcx_cov3(:) = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom x-convection is ", mbcx_cov3

    call update_symmetric_ibc(dm%ibcy_qy, mbc, dm%ibcy_qz)
    mbcy_cov3(:) = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom y-convection is ", mbcy_cov3

    call update_symmetric_ibc(dm%ibcz_qz, mbc, dm%ibcz_qz)
    mbcz_cov3(:) = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom z-convection is ", mbcz_cov3

    if(dm%icoordinate == ICYLINDRICAL) then
      call update_symmetric_ibc(dm%ibcy_qy, mbc, dm%ibcy_qz)
      mbcr_cov3(:) = mbc(:, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom r-convection is ", mbcr_cov3
    end if

    call update_symmetric_ibc(dm%ibcz_qz, mbc)
    bc(:) = mbc(:, JBC_GRAD)
    call update_symmetric_ibc(dm%ibcz_Th, mbc, bc)
    mbcz_tau3 = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom z-viscous is ", mbcz_tau3

    call update_symmetric_ibc(dm%ibcx_qz, mbc)
    bc(:) = mbc(:, JBC_GRAD)
    call update_symmetric_ibc(dm%ibcx_Th, mbc, bc)
    call update_symmetric_ibc(dm%ibcx_Th, mbc0, dm%ibcx_qx)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcx_tau3 is wrong.")
    mbcx_tau3 = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom x-viscous is ", mbcx_tau3

    call update_symmetric_ibc(dm%ibcy_qz, mbc)
    bc(:) = mbc(:, JBC_GRAD)
    call update_symmetric_ibc(dm%ibcy_Th, mbc, bc)
    call update_symmetric_ibc(dm%ibcy_Th, mbc0, dm%ibcz_qy)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcy_tau3 is wrong.")
    mbcy_tau3 = mbc(:, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom y-viscous is ", mbcy_tau3

    if(dm%icoordinate == ICYLINDRICAL) then
      call update_symmetric_ibc(dm%ibcy_Th, mbc, dm%ibcy_qz)
      if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcy_tau3 is wrong.")
      mbcr_tau3 = mbc(:, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom r-viscous is ", mbcr_tau3
    end if

!----------------------------------------------------------------------------------------------------------
!   energy-eqs
!----------------------------------------------------------------------------------------------------------
    if(dm%is_thermo) then
      call update_symmetric_ibc(dm%ibcx_qx, mbc, dm%ibcx_Th)
      ebcx_conv(:) = mbc(:, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for energy x-convection is ", ebcx_conv

      call update_symmetric_ibc(dm%ibcy_qy, mbc, dm%ibcy_Th)
      ebcy_conv(:) = mbc(:, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for energy y-convection is ", ebcy_conv

      call update_symmetric_ibc(dm%ibcz_qz, mbc, dm%ibcz_Th)
      ebcz_conv(:) = mbc(:, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for energy z-convection is ", ebcz_conv

      call update_symmetric_ibc(dm%ibcx_Th, mbc)
      bc(:) = mbc(:, JBC_GRAD)
      call update_symmetric_ibc(dm%ibcx_Th, mbc, bc)
      ebcx_difu = mbc(:, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for energy x-diffusion is ", ebcx_difu

      call update_symmetric_ibc(dm%ibcy_Th, mbc)
      bc(:) = mbc(:, JBC_GRAD)
      call update_symmetric_ibc(dm%ibcy_Th, mbc, bc)
      ebcy_difu = mbc(:, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for energy y-diffusion is ", ebcy_difu

      call update_symmetric_ibc(dm%ibcz_Th, mbc)
      bc(:) = mbc(:, JBC_GRAD)
      call update_symmetric_ibc(dm%ibcz_Th, mbc, bc)
      ebcz_difu = mbc(:, JBC_PROD)
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
!   subroutine apply_bc_const (dm, fl) ! check, necessary?
!     use parameters_constant_mod
!     use udf_type_mod
!     implicit none
!     type(t_domain), intent( in    )   :: dm
!     type(t_flow),   intent( inout )   :: fl

!     integer :: m, s
!     type(DECOMP_INFO) :: dtmp

!     real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil !
!     real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil ! 
!     real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil !


!     !  if(dm%ibcx(1, 1) /= IBC_DIRICHLET .and. &
!     !     dm%ibcx(2, 1) /= IBC_DIRICHLET .and. &
!     !     dm%ibcy(1, 2) /= IBC_DIRICHLET .and. &
!     !     dm%ibcy(2, 2) /= IBC_DIRICHLET .and. &
!     !     dm%ibcz(1, 3) /= IBC_DIRICHLET .and. &
!     !     dm%ibcz(2, 3) /= IBC_DIRICHLET ) return
! !----------------------------------------------------------------------------------------------------------
! !   all in x-pencil
! !----------------------------------------------------------------------------------------------------------

! !----------------------------------------------------------------------------------------------------------
! !   x-pencil, ux stored at nodes, bc is nodes.
! !----------------------------------------------------------------------------------------------------------
!     dtmp = dm%dpcc
!     if(dm%ibcx(1, 1) == IBC_DIRICHLET .and. dtmp%xst(1) == 1) then ! ux at x-begin-xpencil
!       fl%qx     (1, 1:dtmp%xsz(2), 1:dtmp%xsz(3)) = &
!       dm%fbcx_qx(1, 1:dtmp%xsz(2), 1:dtmp%xsz(3))
!     end if
!     if(dm%ibcx(2, 1) == IBC_DIRICHLET .and. dtmp%xen(1) == dm%np(1)) then ! ux at x-end-xpencil
!       fl%qx     (dtmp%xsz(m), 1:dtmp%xsz(2), 1:dtmp%xsz(3)) = &
!       dm%fbcx_qx(2,           1:dtmp%xsz(2), 1:dtmp%xsz(3))
!     end if
! !----------------------------------------------------------------------------------------------------------
! !   x-pencil, uy stored at cells, bc is nodes. symmetric to get bc.
! !----------------------------------------------------------------------------------------------------------
!     dtmp = dm%dcpc
!     if(dm%ibcx_nominal(1, 2) == IBC_DIRICHLET .and. dtmp%xst(1) == 1) then ! ux at x-begin-xpencil
!       dm%ibcx(1, 2) = IBC_INTERIOR
!       dm%fbcx_qy(1, 1:dtmp%xsz(2), 1:dtmp%xsz(3)) = two * dm%fbcx_const(1, 2) - 
!       fl%qy     (1, 1:dtmp%xsz(2), 1:dtmp%xsz(3)) 
!     end if
!     if(dm%ibcx(2, 1) == IBC_DIRICHLET .and. dtmp%xen(1) == dm%np(1)) then ! ux at x-end-xpencil
!       fl%qx     (dtmp%xsz(m), 1:dtmp%xsz(2), 1:dtmp%xsz(3)) = &
!       dm%fbcx_qx(2,           1:dtmp%xsz(2), 1:dtmp%xsz(3))
!     end if

! !----------------------------------------------------------------------------------------------------------
! !   uy at y-direction. BC of others at y-direction are given in operations directly.
! !----------------------------------------------------------------------------------------------------------
!     m = 2
!     dtmp = dm%dcpc
!     do s = 1, 2
!       if(dm%ibcy_nominal(s, m) == IBC_DIRICHLET) then
!         call transpose_x_to_y(fl%qy, acpc_ypencil, dm%dcpc)
!         if(dtmp%yst(m) == 1) then
!           acpc_ypencil(1:dtmp%ysz(1), 1, 1:dtmp%ysz(3)) = &
!             dm%fbcy_qy(1:dtmp%ysz(1), 1, 1:dtmp%ysz(3))
!           !if(dm%is_thermo) fl%gy(:, 1, :) = fl%qy(:, 1, :) * dm%fbc_dend(s, m)
!         end if
!         if(dtmp%yen(m) == dm%np(m)) then
!           acpc_ypencil(1:dtmp%ysz(1), dtmp%ysz(m), 1:dtmp%ysz(3)) = &
!             dm%fbcy_qy(1:dtmp%ysz(1), 2,           1:dtmp%ysz(3))
!           !if(dm%is_thermo) fl%gy(:, dtmp%xsz(m), :) = fl%qy(:, dtmp%xsz(m), :)  * dm%fbc_dend(s, m)
!         end if
!         call transpose_y_to_x(acpc_ypencil, fl%qy, dm%dcpc)
!       end if
!     end do

! !----------------------------------------------------------------------------------------------------------
! !   uz at z-direction. BC of others at z-direction are given in oeprations directly.
! !----------------------------------------------------------------------------------------------------------
!     m = 3
!     dtmp = dm%dccp
!     do s = 1, 2
!       if(dm%ibcz_nominal(s, m) == IBC_DIRICHLET) then
!         call transpose_x_to_y(fl%qz,        accp_ypencil, dm%dccp)
!         call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
!         if(dtmp%zst(m) == 1) then
!           accp_zpencil(1:dtmp%zsz(1), 1:dtmp%zsz(2), 1) = &
!             dm%fbcz_qz(1:dtmp%zsz(1), 1:dtmp%zsz(2), 1)
!           !if(dm%is_thermo) fl%gz(:, :, 1) = fl%qz(:, :, 1) * dm%fbc_dend(s, m)
!         end if
!         if(dtmp%zen(m) == dm%np(m)) then
!           accp_zpencil(1:dtmp%zsz(1), 1:dtmp%zsz(2), dtmp%zsz(m)) = &
!             dm%fbcz_qz(1:dtmp%zsz(1), 1:dtmp%zsz(2), 2)
!           !if(dm%is_thermo)  fl%gz(:, :, dtmp%xsz(m)) = fl%qz(:, :, dtmp%xsz(m)) *  dm%fbc_dend(s, m)
!         end if
!         call transpose_z_to_y(accp_zpencil, accp_ypencil, dm%dccp )
!         call transpose_y_to_x(accp_ypencil, fl%qz, dm%dccp)
!       end if
!     end do

!     return
!   end subroutine


!==========================================================================================================
!==========================================================================================================
!   subroutine Apply_x_convective_outlet ? check necessary?
!     use solver_tools_mod

! !----------------------------------------------------------------------------------------------------------
! !  to get convective outlet bc, ux_cbc
! !----------------------------------------------------------------------------------------------------------
!     call Find_max_min_3d(fl%qx, umax, umin)
!     ux_cbc = HALF * (umax + umin)

!     return
!   end subroutine
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


end module
