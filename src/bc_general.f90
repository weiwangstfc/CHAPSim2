module boundary_conditions_mod
  use bc_dirichlet_mod
  use bc_convective_outlet_mod
  use bc_ndomain_interior_mod
  use udf_type_mod
  use parameters_constant_mod
  use print_msg_mod
  implicit none

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

  private :: reassign_ibc           ! re-assign calcuation ibc and keep the nominal bc
  public  :: config_calc_basic_ibc  ! applied once only, just before calculation

  public  :: allocate_fbc_flow   ! applied once only
  public  :: allocate_fbc_thermo ! applied once only

  private :: get_fbcy_circle_centre
  public  :: update_fbcy_cc_flow_halo   ! for pipe only, applied every NS, cc for circle central point and var stored in xcx
  public  :: update_fbcy_cc_thermo_halo ! for pipe only, applied every NS, cc for circle central point and var stored in xcx

  public  :: build_bc_symm_operation    ! applied if necessary
  public  :: config_calc_eqs_ibc

  public  :: get_fbcx_iTh
  public  :: get_fbcy_iTh
  public  :: get_fbcz_iTh

contains
!==========================================================================================================
!==========================================================================================================
  subroutine reassign_ibc(bc_nominal, ibc)
    integer, intent(in) :: bc_nominal(2, 5)
    integer, intent(out) :: ibc(2, 5)
    integer :: n, m

    do n = 1, 2
      do m = 1, 5
        if (bc_nominal(n, m) == IBC_PROFILE1D)   then
          ibc(n, m) = IBC_DIRICHLET
        else if (bc_nominal(n, m) == IBC_TURBGEN )   then
          if(m == 5) then
            ibc(n, m) = IBC_DIRICHLET ! for temperature, default is no incoming thermal flow, it is initilazed temperature
          else if(m == 4) then
            ibc(n, m) = IBC_NEUMANN    ! for p
          else 
            ibc(n, m) = IBC_DIRICHLET  ! for u, v, w
          end if
        else if (bc_nominal(n, m) == IBC_DATABASE )   then
          if(m == 5) then
            ibc(n, m) = IBC_DIRICHLET ! for temperature, default is no incoming thermal flow, it is initilazed temperature
          else if(m == 4) then
            ibc(n, m) = IBC_NEUMANN    ! for p
          else 
            ibc(n, m) = IBC_DIRICHLET  ! for u, v, w, p, check!!
          end if
        else if (bc_nominal(n, m) == IBC_CONVECTIVE)   then ! check for convetive outlet
          ibc(n, m) = IBC_DIRICHLET
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
      dm%ibcx_ftp(n) = dm%ibcx_Th(n)
      dm%ibcy_ftp(n) = dm%ibcy_Th(n)
      dm%ibcz_ftp(n) = dm%ibcz_Th(n)
      if(dm%ibcx_Th(n) == IBC_NEUMANN) dm%ibcx_ftp(n) = IBC_DIRICHLET
      if(dm%ibcy_Th(n) == IBC_NEUMANN) dm%ibcy_ftp(n) = IBC_DIRICHLET
      if(dm%ibcz_Th(n) == IBC_NEUMANN) dm%ibcz_ftp(n) = IBC_DIRICHLET
    end do 

    if(nrank == 0) then

      write (*, wrtfmt1s) '  Boundary type options : '
      write (*, wrtfmt1s) '          0  = IBC_INTERIOR'
      write (*, wrtfmt1s) '          1  = IBC_PERIODIC'
      write (*, wrtfmt1s) '          2  = IBC_SYMMETRIC'
      write (*, wrtfmt1s) '          3  = IBC_ASYMMETRIC'
      write (*, wrtfmt1s) '          4  = IBC_DIRICHLET'
      write (*, wrtfmt1s) '          5  = IBC_NEUMANN'
      write (*, wrtfmt1s) '          6  = IBC_INTRPL'
      write (*, wrtfmt1s) '          7  = IBC_CONVECTIVE'
      write (*, wrtfmt1s) '          8  = IBC_TURBGEN'
      write (*, wrtfmt1s) '          9  = IBC_PROFILE1D'
      write (*, wrtfmt1s) '          10 = IBC_DATABASE'

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

    if(dm%is_record_xoutlet) then
      allocate (dm%fbcx_qx_outl1(dm%dxcc%xsz(1), dm%dxcc%xsz(2), dm%dxcc%xsz(3)) )
      allocate (dm%fbcx_qx_outl2(dm%dxcc%xsz(1), dm%dxcc%xsz(2), dm%dxcc%xsz(3)) )
      allocate (dm%fbcx_qy_outl1(dm%dxpc%xsz(1), dm%dxpc%xsz(2), dm%dxpc%xsz(3)) )
      allocate (dm%fbcx_qy_outl2(dm%dxpc%xsz(1), dm%dxpc%xsz(2), dm%dxpc%xsz(3)) )
      allocate (dm%fbcx_qz_outl1(dm%dxcp%xsz(1), dm%dxcp%xsz(2), dm%dxcp%xsz(3)) )
      allocate (dm%fbcx_qz_outl2(dm%dxcp%xsz(1), dm%dxcp%xsz(2), dm%dxcp%xsz(3)) )
      allocate (dm%fbcx_pr_outl1(dm%dxcc%xsz(1), dm%dxcc%xsz(2), dm%dxcc%xsz(3)) )
      allocate (dm%fbcx_pr_outl2(dm%dxcc%xsz(1), dm%dxcc%xsz(2), dm%dxcc%xsz(3)) )
    end if

    if(dm%is_read_xinlet) then
      allocate (dm%fbcx_qx_inl1(dm%dxcc%xsz(1), dm%dxcc%xsz(2), dm%dxcc%xsz(3)) )
      allocate (dm%fbcx_qx_inl2(dm%dxcc%xsz(1), dm%dxcc%xsz(2), dm%dxcc%xsz(3)) )
      allocate (dm%fbcx_qy_inl1(dm%dxpc%xsz(1), dm%dxpc%xsz(2), dm%dxpc%xsz(3)) )
      allocate (dm%fbcx_qy_inl2(dm%dxpc%xsz(1), dm%dxpc%xsz(2), dm%dxpc%xsz(3)) )
      allocate (dm%fbcx_qz_inl1(dm%dxcp%xsz(1), dm%dxcp%xsz(2), dm%dxcp%xsz(3)) )
      allocate (dm%fbcx_qz_inl2(dm%dxcp%xsz(1), dm%dxcp%xsz(2), dm%dxcp%xsz(3)) )
      allocate (dm%fbcx_pr_inl1(dm%dxcc%xsz(1), dm%dxcc%xsz(2), dm%dxcc%xsz(3)) )
      allocate (dm%fbcx_pr_inl2(dm%dxcc%xsz(1), dm%dxcc%xsz(2), dm%dxcc%xsz(3)) )
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
! to calculate boundary during calculation from primary boundary
  subroutine build_bc_symm_operation(ibc, mbc, jbc)
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
            if(ibc(i) == IBC_DIRICHLET) mbc(i, :) = ibc(i)
            if(jbc(i) == IBC_DIRICHLET) mbc(i, :) = jbc(i)
            if(ibc(i) == IBC_PERIODIC .or. jbc(i) == IBC_PERIODIC) then
              if(nrank==0) write(*, '(A20, I2.1, A5, I2.1)') "BCs for the side ", i, " are ", ibc(i), jbc(i) 
              call Print_warning_msg("The two operational variables have different boundary conditions.")
            end if
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
    integer :: bc(2), n
!----------------------------------------------------------------------------------------------------------
!   x-mom
!----------------------------------------------------------------------------------------------------------
    call build_bc_symm_operation(dm%ibcx_qx, mbc, dm%ibcx_qx)
    mbcx_cov1(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom x-convection is ", mbcx_cov1 

    call build_bc_symm_operation(dm%ibcy_qy, mbc, dm%ibcy_qx)
    mbcy_cov1(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom y-convection is ", mbcy_cov1

    call build_bc_symm_operation(dm%ibcz_qz, mbc, dm%ibcz_qx)
    mbcz_cov1(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom z-convection is ", mbcz_cov1

    call build_bc_symm_operation(dm%ibcx_qx, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call build_bc_symm_operation(dm%ibcx_ftp, mbc, bc)
    mbcx_tau1(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom x-diffusion  is ", mbcx_tau1

    call build_bc_symm_operation(dm%ibcy_qx, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call build_bc_symm_operation(dm%ibcy_ftp, mbc, bc)
    call build_bc_symm_operation(dm%ibcy_ftp, mbc0, dm%ibcy_qy)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcy_tau1 is wrong.")
    mbcy_tau1(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom y-diffusion  is ", mbcy_tau1

    call build_bc_symm_operation(dm%ibcz_qx, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call build_bc_symm_operation(dm%ibcz_ftp, mbc, bc)
    call build_bc_symm_operation(dm%ibcz_ftp, mbc0, dm%ibcz_qz)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcy_tau1 is wrong.")
    mbcz_tau1(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for x-mom z-diffusion  is ", mbcz_tau1
!----------------------------------------------------------------------------------------------------------
!   y-mom
!----------------------------------------------------------------------------------------------------------
    call build_bc_symm_operation(dm%ibcx_qx, mbc, dm%ibcx_qy)
    mbcx_cov2(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom x-convection is ", mbcx_cov2

    call build_bc_symm_operation(dm%ibcy_qy, mbc, dm%ibcy_qy)
    mbcy_cov2(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom y-convection is ", mbcy_cov2

    call build_bc_symm_operation(dm%ibcz_qz, mbc, dm%ibcz_qy)
    mbcz_cov2(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom z-convection is ", mbcz_cov2

    if(dm%icoordinate == ICYLINDRICAL) then
      call build_bc_symm_operation(dm%ibcy_qz, mbc, dm%ibcy_qz)
      mbcr_cov2(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom r-convection is ", mbcr_cov2
    end if

    call build_bc_symm_operation(dm%ibcx_qy, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call build_bc_symm_operation(dm%ibcx_ftp, mbc, bc)
    call build_bc_symm_operation(dm%ibcx_ftp, mbc0, dm%ibcx_qx)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcx_tau2 is wrong.")
    mbcx_tau2(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom x-diffusion  is ", mbcx_tau2

    call build_bc_symm_operation(dm%ibcy_qy, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call build_bc_symm_operation(dm%ibcy_ftp, mbc, bc)
    mbcy_tau2(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom y-diffusion  is ", mbcy_tau2

    call build_bc_symm_operation(dm%ibcz_qy, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call build_bc_symm_operation(dm%ibcz_ftp, mbc, bc)
    call build_bc_symm_operation(dm%ibcz_ftp, mbc0, dm%ibcz_qz)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcz_tau2 is wrong.")
    mbcz_tau2(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom z-diffusion  is ", mbcz_tau2

    if(dm%icoordinate == ICYLINDRICAL) then
      call build_bc_symm_operation(dm%ibcy_qz, mbc, dm%ibcy_ftp)
      mbcr_tau2(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for y-mom r-diffusion  is ", mbcr_tau2
    end if
!----------------------------------------------------------------------------------------------------------
!   z-mom
!----------------------------------------------------------------------------------------------------------
    call build_bc_symm_operation(dm%ibcx_qx, mbc, dm%ibcx_qz)
    mbcx_cov3(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom x-convection is ", mbcx_cov3

    call build_bc_symm_operation(dm%ibcy_qy, mbc, dm%ibcy_qz)
    mbcy_cov3(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom y-convection is ", mbcy_cov3

    call build_bc_symm_operation(dm%ibcz_qz, mbc, dm%ibcz_qz)
    mbcz_cov3(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom z-convection is ", mbcz_cov3

    if(dm%icoordinate == ICYLINDRICAL) then
      call build_bc_symm_operation(dm%ibcy_qy, mbc, dm%ibcy_qz)
      mbcr_cov3(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom r-convection is ", mbcr_cov3
    end if

    call build_bc_symm_operation(dm%ibcx_qz, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call build_bc_symm_operation(dm%ibcx_ftp, mbc, bc)
    call build_bc_symm_operation(dm%ibcx_ftp, mbc0, dm%ibcx_qx)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcx_tau3 is wrong.")
    mbcx_tau3(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom x-diffusion  is ", mbcx_tau3

    call build_bc_symm_operation(dm%ibcy_qz, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call build_bc_symm_operation(dm%ibcy_ftp, mbc, bc)
    call build_bc_symm_operation(dm%ibcy_ftp, mbc0, dm%ibcy_qy)
    if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcy_tau3 is wrong.")
    mbcy_tau3(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom y-diffusion  is ", mbcy_tau3

    call build_bc_symm_operation(dm%ibcz_qz, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call build_bc_symm_operation(dm%ibcz_ftp, mbc, bc)
    mbcz_tau3 = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom z-diffusion  is ", mbcz_tau3

    if(dm%icoordinate == ICYLINDRICAL) then
      call build_bc_symm_operation(dm%ibcy_ftp, mbc, dm%ibcy_qz)
      if(mbc0(1, JBC_PROD)/= mbc(1, JBC_PROD)) call Print_error_msg("BC in mbcy_tau3 is wrong.")
      mbcr_tau3(1:2) = mbc(1:2, JBC_PROD)
      if(nrank==0) write(*, wrtfmt2i) "The bc for z-mom r-diffusion  is ", mbcr_tau3
    end if
!----------------------------------------------------------------------------------------------------------
!   energy-eqs
!----------------------------------------------------------------------------------------------------------
    if(dm%is_thermo)  then
    call build_bc_symm_operation(dm%ibcx_qx, mbc, dm%ibcx_ftp)
    ebcx_conv(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for energy x-convection is ", ebcx_conv

    call build_bc_symm_operation(dm%ibcy_qy, mbc, dm%ibcy_ftp)
    ebcy_conv(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for energy y-convection is ", ebcy_conv

    call build_bc_symm_operation(dm%ibcz_qz, mbc, dm%ibcz_ftp)
    ebcz_conv(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for energy z-convection is ", ebcz_conv

    call build_bc_symm_operation(dm%ibcx_Th, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call build_bc_symm_operation(dm%ibcx_ftp, mbc, bc)
    ebcx_difu = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for energy x-diffusion  is ", ebcx_difu

    call build_bc_symm_operation(dm%ibcy_Th, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call build_bc_symm_operation(dm%ibcy_ftp, mbc, bc)
    ebcy_difu(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for energy y-diffusion  is ", ebcy_difu

    call build_bc_symm_operation(dm%ibcz_Th, mbc)
    bc(1:2) = mbc(1:2, JBC_GRAD)
    call build_bc_symm_operation(dm%ibcz_ftp, mbc, bc)
    ebcz_difu(1:2) = mbc(1:2, JBC_PROD)
    if(nrank==0) write(*, wrtfmt2i) "The bc for energy z-diffusion  is ", ebcz_difu
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
!==========================================================================================================
  subroutine get_fbcx_iTh(ibc, dm, fbc)
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    integer, intent(in) :: ibc(2)
    type(t_domain), intent(in) :: dm
    real(WP), intent(out) :: fbc(4, dm%dpcc%xsz(2), dm%dpcc%xsz(3))

    integer :: n

    fbc = ZERO
    do n = 1, 2
      if(ibc(n) == IBC_DIRICHLET) then    
        fbc(n, :, :) = dm%fbcx_ftp(n, :, :)%t
      else if(ibc(n) == IBC_NEUMANN) then
        fbc(n, :, :) = dm%fbcx_qw(n, :, :)
      else
        fbc(n, :, :) = ZERO
      end if
    end do 
    return
  end subroutine 
!==========================================================================================================
  subroutine get_fbcy_iTh(ibc, dm, fbc)
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    integer, intent(in) :: ibc(2)
    type(t_domain), intent(in) :: dm
    real(WP), intent(out) :: fbc(dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3))

    integer :: n

    fbc = ZERO
    do n = 1, 2
      if(ibc(n) == IBC_DIRICHLET) then    
        fbc(:, n, :) = dm%fbcy_ftp(:, n, :)%t
      else if(ibc(n) == IBC_NEUMANN) then
        fbc(:, n, :) = dm%fbcy_qw(:, n, :)
      else
        fbc(:, n, :) = ZERO
      end if
    end do 
    return
  end subroutine 
!==========================================================================================================
  subroutine get_fbcz_iTh(ibc, dm, fbc)
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    integer, intent(in) :: ibc(2)
    type(t_domain), intent(in) :: dm
    real(WP), intent(out) :: fbc(dm%dccp%zsz(1), dm%dccp%zsz(2), 4)

    integer :: n

    fbc = ZERO
    do n = 1, 2
      if(ibc(n) == IBC_DIRICHLET) then    
        fbc(:, :, n) = dm%fbcz_ftp(:, :, n)%t
      else if(ibc(n) == IBC_NEUMANN) then
        fbc(:, :, n) = dm%fbcz_qw(:, :, n)
      else
        fbc(:, :, n) = ZERO
      end if
    end do 
    return
  end subroutine 

!==========================================================================================================
!==========================================================================================================
  subroutine extract_dirichlet_fbcx(fbc, var, dtmp)
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), intent(out) :: fbc(4,           dtmp%xsz(2), dtmp%xsz(3))
    real(WP), intent(in)  :: var(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3))

    if(dtmp%xsz(1) /= dtmp%xen(1)) call Print_error_msg("Error. This is not x-pencil.")
    fbc(1,   :, :) = var(1,           :, :)
    fbc(2,   :, :) = var(dtmp%xsz(1), :, :)
    fbc(3:4, :, :) = fbc(1:2,         :, :)

    return
  end subroutine 
!==========================================================================================================
  subroutine extract_dirichlet_fbcy(fbc, var, dtmp)
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), intent(out) :: fbc(dtmp%ysz(1), 4,           dtmp%ysz(3))
    real(WP), intent(in)  :: var(dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3))

    if(dtmp%ysz(2) /= dtmp%yen(2)) call Print_error_msg("Error. This is not y-pencil.")
    fbc(:, 1,   :) = var(:, 1,           :)
    fbc(:, 2,   :) = var(:, dtmp%ysz(2), :)
    fbc(:, 3:4, :) = fbc(:, 1:2,         :)

    return
  end subroutine 
!==========================================================================================================
  subroutine extract_dirichlet_fbcz(fbc, var, dtmp)
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), intent(in)  :: var(dtmp%zsz(1), dtmp%zsz(2), dtmp%zsz(3))
    real(WP), intent(out) :: fbc(dtmp%zsz(1), dtmp%zsz(2), 4          )
    

    if(dtmp%zsz(3) /= dtmp%zen(3)) call Print_error_msg("Error. This is not z-pencil.")
    fbc(:, :, 1, ) = var(:, :, 1,         )
    fbc(:, :, 2, ) = var(:, :, dtmp%zsz(3))
    fbc(:, :, 3:4) = fbc(:, :, 1:2,       )

    return
  end subroutine 

!==========================================================================================================


!==========================================================================================================
end module
