module test_operations_mod
  use operations

  public  :: Test_interpolation
  public  :: Test_1st_derivative
  public  :: Test_2nd_derivative
  
contains
!===============================================================================
!===============================================================================
!> \brief To test this subroutine for mid-point interpolation.
!>
!> This subroutine is called in \ref Test_algorithms. Define the logicals to choose
!> which test section is required. 
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!_______________________________________________________________________________
  subroutine Test_interpolation(dm)
    use operations
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    integer :: i, j, k
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    real(WP) :: ref
    real(WP) :: err, err_Linf, err_L2
    integer  :: wrt_unit

    real(WP) :: fxc (dm%nc(1)), fyc (dm%nc(2)), fzc (dm%nc(3))
    real(WP) :: fxp (dm%np(1)), fyp (dm%np(2)), fzp (dm%np(3))
    real(WP) :: fgxc(dm%nc(1)), fgyc(dm%nc(2)), fgzc(dm%nc(3))
    real(WP) :: fgxp(dm%np(1)), fgyp(dm%np(2)), fgzp(dm%np(3))

    real(WP) :: scale, shift

    open (newunit = wrt_unit, file = 'check_test_algorithms.dat', position="append")

    if (dm%ibcx(1, 5) == IBC_PERIODIC) then
      scale = ONE
      shift = ZERO
    else if (dm%ibcx(1, 5) == IBC_SYMMETRIC) then
      scale = ONE
      shift = PI / TWO
    else if (dm%ibcx(1, 5) == IBC_ASYMMETRIC) then
      scale = TWO
      shift = ZERO
    else if (dm%ibcx(1, 5) == IBC_DIRICHLET) then
      scale = THREE
      shift = ZERO
      dm%fbcx(1, 5) = ZERO
      dm%fbcx(2, 5) = sin_wp(TWO * PI / THREE)
      dm%fbcy(:, 5) = dm%fbcx(:, 5)
      dm%fbcz(:, 5) = dm%fbcx(:, 5)
    else if (dm%ibcx(1, 5) == IBC_NEUMANN) then
      scale = THREE
      shift = ZERO
      dm%fbcx(1, 5) = ONE / THREE * cos_wp(ZERO / THREE)
      dm%fbcx(2, 5) = ONE / THREE * cos_wp(TWO * PI / THREE)
      dm%fbcy(:, 5) = dm%fbcx(:, 5)
      dm%fbcz(:, 5) = dm%fbcx(:, 5)
    else 
      scale = THREE
      shift = ZERO
      dm%fbcx(1, 5) = ZERO
      dm%fbcx(2, 5) = sin_wp(TWO * PI / THREE)
      dm%fbcy(:, 5) = dm%fbcx(:, 5)
      dm%fbcz(:, 5) = dm%fbcx(:, 5)
    end if


! x direction, xc
    do i = 1, dm%nc(1)
      xc = dm%h(1) * (real(i - 1, WP) + HALF)
      fxc(i) = sin_wp ( xc / scale + shift)
    end do
! x: c2p
    call Get_x_midp_C2P_1D (fxc, fgxp, dm, dm%ibcx(:, 5), dm%fbcx(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do i = 1, dm%np(1)
      xp = dm%h(1) * real(i - 1, WP)
      ref = sin_wp(xp / scale + shift)
      err = dabs(fgxp(i) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'x-interp-c2p ', i, xp, ref, fgxp(i), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(1)) 
    write(wrt_unit, *) 'x-interp-c2p ', dm%np(1), err_Linf, err_L2

! y direction, yc
    do j = 1, dm%nc(2)
      yc = dm%yc(j)
      fyc(j) = sin_wp ( yc / scale + shift)
    end do
! y: c2p
    call Get_y_midp_C2P_1D (fyc, fgyp, dm, dm%ibcy(:, 5), dm%fbcy(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do j = 1, dm%np(2)
      yp = dm%yp(j)
      ref = sin_wp(yp / scale + shift)
      err = dabs(fgyp(j) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'y-interp-c2p ', j, yp, ref, fgyp(j), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(2)) 
    write(wrt_unit, *) 'y-interp-c2p ', dm%np(2), err_Linf, err_L2

 ! z direction, zc
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      fzc(k) = sin_wp ( zc / scale + shift)
    end do
! z: c2p
    call Get_z_midp_C2P_1D (fzc, fgzp, dm, dm%ibcz(:, 5), dm%fbcz(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      ref = sin_wp(zp / scale + shift)
      err = dabs(fgzp(k) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'z-interp-c2p ', k, zp, ref, fgzp(k), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(3)) 
    write(wrt_unit, *) 'z-interp-c2p ', dm%np(3), err_Linf, err_L2

! x direction, xp
    do i = 1, dm%np(1)
      xp = dm%h(1) * real(i - 1, WP)
      fxp(i) = sin_wp ( xp / scale + shift)
    end do
! x: p2c
    call Get_x_midp_P2C_1D (fxp, fgxc, dm, dm%ibcx(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do i = 1, dm%nc(1)
      xc = dm%h(1) * (real(i - 1, WP) + HALF)
      ref = sin_wp(xc / scale + shift)
      err = dabs(fgxc(i) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'x-interp-p2c ', i, xc, ref, fgxc(i), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(1)) 
    write(wrt_unit, *) 'x-interp-p2c ', dm%nc(1), err_Linf, err_L2

! y direction, yp
    do j = 1, dm%np(2)
      yp = dm%yp(j)
      fyp(j) = sin_wp ( yp / scale + shift)
    end do
! y: p2c
    call Get_y_midp_P2C_1D (fyp, fgyc, dm, dm%ibcy(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do j = 1, dm%nc(2)
      yc = dm%yc(j)
      ref = sin_wp(yc / scale + shift)
      err = dabs(fgyc(j) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'y-interp-p2c ', j, yc, ref, fgyc(j), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(2)) 
    write(wrt_unit, *) 'y-interp-p2c ', dm%nc(2), err_Linf, err_L2


! z direction, zp
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      fzp(k) = sin_wp ( zp / scale + shift)
    end do
! z: p2c
    call Get_z_midp_P2C_1D (fzp, fgzc, dm, dm%ibcz(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      ref = sin_wp(zc / scale + shift)
      err = dabs(fgzc(k) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'z-interp-p2c ', k, zc, ref, fgzc(k), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(3)) 
    write(wrt_unit, *) 'z-interp-p2c ', dm%nc(3), err_Linf, err_L2
    close(wrt_unit)

    return 
  end subroutine

!===============================================================================
!===============================================================================
!> \brief To test this subroutine for mid-point interpolation.
!>
!> This subroutine is called in \ref Test_algorithms. Define the logicals to choose
!> which test section is required. 
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!_______________________________________________________________________________
  subroutine Test_1st_derivative(dm)
    use operations
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    integer :: i, j, k
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    real(WP) :: ref
    real(WP) :: err, err_Linf, err_L2
    integer  :: wrt_unit

    real(WP) :: fxc (dm%nc(1)), fyc (dm%nc(2)), fzc (dm%nc(3))
    real(WP) :: fxp (dm%np(1)), fyp (dm%np(2)), fzp (dm%np(3))
    real(WP) :: fgxc(dm%nc(1)), fgyc(dm%nc(2)), fgzc(dm%nc(3))
    real(WP) :: fgxp(dm%np(1)), fgyp(dm%np(2)), fgzp(dm%np(3))
    
    real(WP) :: scale, shift

    if (dm%ibcx(1, 5) == IBC_PERIODIC) then
      scale = ONE
      shift = ZERO
    else if (dm%ibcx(1, 5) == IBC_SYMMETRIC) then
      scale = ONE
      shift = PI / TWO
    else if (dm%ibcx(1, 5) == IBC_ASYMMETRIC) then
      scale = TWO
      shift = ZERO
    else if (dm%ibcx(1, 5) == IBC_DIRICHLET) then
      scale = THREE
      shift = ZERO
      dm%fbcx(1, 5) = ZERO
      dm%fbcx(2, 5) = sin_wp(TWO * PI / THREE)
      dm%fbcy(:, 5) = dm%fbcx(:, 5)
      dm%fbcz(:, 5) = dm%fbcx(:, 5)
    else if (dm%ibcx(1, 5) == IBC_NEUMANN) then
      scale = THREE
      shift = ZERO
      dm%fbcx(1, 5) = ONE / THREE * cos_wp(ZERO / THREE + shift)
      dm%fbcx(2, 5) = ONE / THREE * cos_wp(TWO * PI / THREE + shift)
      dm%fbcy(:, 5) = dm%fbcx(:, 5)
      dm%fbcz(:, 5) = dm%fbcx(:, 5)
    else 
      scale = THREE
      shift = ZERO
      dm%fbcx(1, 5) = ZERO
      dm%fbcx(2, 5) = sin_wp(TWO * PI / THREE)
      dm%fbcy(:, 5) = dm%fbcx(:, 5)
      dm%fbcz(:, 5) = dm%fbcx(:, 5)
    end if

    open (newunit = wrt_unit, file = 'check_test_algorithms.dat', position="append")

! x direction, xc
    do i = 1, dm%nc(1)
      xc = dm%h(1) * (real(i - 1, WP) + HALF)
      fxc(i) = sin_wp ( xc / scale + shift)
    end do
! x: c2c
    call Get_x_1st_derivative_C2C_1D (fxc, fgxc, dm, dm%ibcx(:, 5), dm%fbcx(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do i = 1, dm%nc(1)
      xc = dm%h(1) * (real(i - 1, WP) + HALF)
      ref = ONE/scale * cos_wp(xc / scale + shift)
      err = dabs(fgxc(i) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'x-1stder-c2c ', i, xc, ref, fgxc(i), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(1)) 
    write(wrt_unit, *) 'x-1stder-c2c ', dm%nc(1), err_Linf, err_L2

! y direction, yc
    do j = 1, dm%nc(2)
      yc = dm%yc(j)
      fyc(j) = sin_wp ( yc / scale + shift)
    end do
! y: c2c
    call Get_y_1st_derivative_C2C_1D (fyc, fgyc, dm, dm%ibcy(:, 5), dm%fbcy(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do j = 1, dm%nc(2)
      yc = dm%yc(j)
      ref = ONE/scale * cos_wp(yc / scale + shift)
      err = dabs(fgyc(j) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'y-1stder-c2c ', j, yc, ref, fgyc(j), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(2)) 
    write(wrt_unit, *) 'y-1stder-c2c ', dm%nc(2), err_Linf, err_L2

! z direction, zc
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      fzc(k) = sin_wp ( zc / scale + shift)
    end do
! z: c2c
    call Get_z_1st_derivative_C2C_1D (fzc, fgzc, dm, dm%ibcz(:, 5), dm%fbcz(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      ref = ONE/scale * cos_wp(zc / scale + shift)
      err = dabs(fgzc(k) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'z-1stder-c2c ', k, zc, ref, fgzc(k), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(3)) 
    write(wrt_unit, *) 'z-1stder-c2c ', dm%nc(3), err_Linf, err_L2

 ! x direction, xp
    do i = 1, dm%np(1)
      xp = dm%h(1) * real(i - 1, WP)
      fxp(i) = sin_wp ( xp / scale + shift)
    end do
! x: p2p
    call Get_x_1st_derivative_P2P_1D (fxp, fgxp, dm, dm%ibcx(:, 5), dm%fbcx(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do i = 1, dm%np(1)
      xp = dm%h(1) * real(i - 1, WP)
      ref = ONE/scale * cos_wp(xp / scale + shift)
      err = dabs(fgxp(i) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'x-1stder-p2p', i, xp, ref, fgxp(i), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(1)) 
    write(wrt_unit, *) 'x-1stder-p2p ', dm%np(1), err_Linf, err_L2
! y direction, yp
    do j = 1, dm%np(2)
      yp = dm%yp(j)
      fyp(j) = sin_wp ( yp / scale + shift)
    end do
! y: p2p
    call Get_y_1st_derivative_P2P_1D (fyp, fgyp, dm, dm%ibcy(:, 5), dm%fbcy(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do j = 1, dm%np(2)
      yp = dm%yp(j)
      ref = ONE/scale * cos_wp(yp / scale + shift)
      err = dabs(fgyp(j) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'y-1stder-p2p ', j, yp, ref, fgyp(j), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(2)) 
    write(wrt_unit, *) 'y-1stder-p2p ', dm%np(2), err_Linf, err_L2
! z direction, zp
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      fzp(k) = sin_wp ( zp / scale + shift)
    end do
! z: p2p
    call Get_z_1st_derivative_P2P_1D (fzp, fgzp, dm, dm%ibcz(:, 5), dm%fbcz(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      ref = ONE/scale * cos_wp(zp / scale + shift)
      err = dabs(fgzp(k) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'z-1stder-p2p ', k, zp, ref, fgzp(k), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(3)) 
    write(wrt_unit, *) 'z-1stder-p2p ', dm%np(3), err_Linf, err_L2

! x: c2p
    call Get_x_1st_derivative_C2P_1D (fxc, fgxp, dm, dm%ibcx(:, 5), dm%fbcx(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do i = 1, dm%np(1)
      xp = dm%h(1) * real(i - 1, WP)
      ref = ONE/scale * cos_wp(xp / scale + shift)
      err = dabs(fgxp(i) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'x-1stder-c2p ', i, xp, ref, fgxp(i), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(1)) 
    write(wrt_unit, *) 'x-1stder-c2p ', dm%np(1), err_Linf, err_L2

! y: c2p
    call Get_y_1st_derivative_C2P_1D (fyc, fgyp, dm, dm%ibcy(:, 5), dm%fbcy(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do j = 1, dm%np(2)
      yp = dm%yp(j)
      ref = ONE/scale * cos_wp(yp / scale + shift)
      err = dabs(fgyp(j) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'y-1stder-c2p ', j, yp, ref, fgyp(j), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(2)) 
    write(wrt_unit, *) 'y-1stder-c2p ', dm%np(2), err_Linf, err_L2

! z: c2p
    call Get_z_1st_derivative_C2P_1D (fzc, fgzp, dm, dm%ibcz(:, 5), dm%fbcz(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      ref = ONE/scale * cos_wp(zp / scale + shift)
      err = dabs(fgzp(k) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'z-1stder-c2p ', k, zp, ref, fgzp(k), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(3)) 
    write(wrt_unit, *) 'z-1stder-c2p ', dm%np(3), err_Linf, err_L2

! x: p2c
    call Get_x_1st_derivative_P2C_1D (fxp, fgxc, dm, dm%ibcx(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do i = 1, dm%nc(1)
      xc = dm%h(1) * (real(i - 1, WP) + HALF)
      ref = ONE/scale * cos_wp(xc / scale + shift)
      err = dabs(fgxc(i) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'x-1stder-p2c ', i, xc, ref, fgxc(i), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(1)) 
    write(wrt_unit, *) 'x-1stder-p2c ', dm%nc(1), err_Linf, err_L2

! y: p2c
    call Get_y_1st_derivative_P2C_1D (fyp, fgyc, dm, dm%ibcy(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do j = 1, dm%nc(2)
      yc = dm%yc(j)
      ref = ONE/scale * cos_wp(yc / scale + shift)
      err = dabs(fgyc(j) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'y-1stder-p2c ', j, yc, ref, fgyc(j), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(2)) 
    write(wrt_unit, *) 'y-1stder-p2c ', dm%nc(2), err_Linf, err_L2

! z: p2c
    call Get_z_1st_derivative_P2C_1D (fzp, fgzc, dm, dm%ibcz(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      ref = ONE/scale * cos_wp(zc / scale + shift)
      err = dabs(fgzc(k) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'z-1stder-p2c ', k, zc, ref, fgzc(k), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(3)) 
    write(wrt_unit, *) 'z-1stder-p2c ', dm%nc(3), err_Linf, err_L2

    close(wrt_unit)

    return 
  end subroutine

!===============================================================================
!===============================================================================
!> \brief To test this subroutine for mid-point interpolation.
!>
!> This subroutine is called in \ref Test_algorithms. Define the logicals to choose
!> which test section is required. 
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!_______________________________________________________________________________
  subroutine Test_2nd_derivative(dm)
    use operations
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    integer :: i, j, k
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    real(WP) :: ref
    real(WP) :: err, err_Linf, err_L2
    integer  :: wrt_unit

    real(WP) :: fxc (dm%nc(1)), fyc (dm%nc(2)), fzc (dm%nc(3))
    real(WP) :: fxp (dm%np(1)), fyp (dm%np(2)), fzp (dm%np(3))
    real(WP) :: fgxc(dm%nc(1)), fgyc(dm%nc(2)), fgzc(dm%nc(3))
    real(WP) :: fgxp(dm%np(1)), fgyp(dm%np(2)), fgzp(dm%np(3))
    
    real(WP) :: scale, shift

    if (dm%ibcx(1, 5) == IBC_PERIODIC) then
      scale = ONE
      shift = ZERO
    else if (dm%ibcx(1, 5) == IBC_SYMMETRIC) then
      scale = ONE
      shift = PI / TWO
    else if (dm%ibcx(1, 5) == IBC_ASYMMETRIC) then
      scale = TWO
      shift = ZERO
    else if (dm%ibcx(1, 5) == IBC_DIRICHLET) then
      scale = THREE
      shift = ZERO
      dm%fbcx(1, 5) = ZERO
      dm%fbcx(2, 5) = sin_wp(TWO * PI / THREE)
      dm%fbcy(:, 5) = dm%fbcx(:, 5)
      dm%fbcz(:, 5) = dm%fbcx(:, 5)
    else if (dm%ibcx(1, 5) == IBC_NEUMANN) then
      scale = THREE
      shift = ZERO
      dm%fbcx(1, 5) = ONE / THREE * cos_wp(ZERO / THREE)
      dm%fbcx(2, 5) = ONE / THREE * cos_wp(TWO * PI / THREE)
      dm%fbcy(:, 5) = dm%fbcx(:, 5)
      dm%fbcz(:, 5) = dm%fbcx(:, 5)
    else 
      scale = THREE
      shift = ZERO
      dm%fbcx(1, 5) = ZERO
      dm%fbcx(2, 5) = sin_wp(TWO * PI / THREE)
      dm%fbcy(:, 5) = dm%fbcx(:, 5)
      dm%fbcz(:, 5) = dm%fbcx(:, 5)
    end if

    open (newunit = wrt_unit, file = 'check_test_algorithms.dat', position="append")

! x direction, xc
    do i = 1, dm%nc(1)
      xc = dm%h(1) * (real(i - 1, WP) + HALF)
      fxc(i) = sin_wp ( xc / scale + shift)
    end do
! x: c2c
    call Get_x_2nd_derivative_C2C_1D (fxc, fgxc, dm, dm%ibcx(:, 5), dm%fbcx(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do i = 1, dm%nc(1)
      xc = dm%h(1) * (real(i - 1, WP) + HALF)
      ref = - (ONE/scale)**2 * sin_wp(xc / scale + shift)
      err = dabs(fgxc(i) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'x-2ndder-c2c ', i, xc, ref, fgxc(i), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(1)) 
    write(wrt_unit, *) 'x-2ndder-c2c ', dm%nc(1), err_Linf, err_L2

! y direction, yc
    do j = 1, dm%nc(2)
      yc = dm%yc(j)
      fyc(j) = sin_wp ( yc / scale + shift)
    end do
! y: c2c
    call Get_y_2nd_derivative_C2C_1D (fyc, fgyc, dm, dm%ibcy(:, 5), dm%fbcy(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do j = 1, dm%nc(2)
      yc = dm%yc(j)
      ref = - (ONE/scale)**2 * sin_wp(yc / scale + shift)
      err = dabs(fgyc(j) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
     ! !write(*,'(A,1I5.1,4ES13.5)') 'y-2ndder-c2c ', j, yc, ref, fgyc(j), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(2)) 
    write(wrt_unit, *) 'y-2ndder-c2c ', dm%nc(2), err_Linf, err_L2

! z direction, zc
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      fzc(k) = sin_wp ( zc / scale + shift)
    end do
! z: c2c
    call Get_z_2nd_derivative_C2C_1D (fzc, fgzc, dm, dm%ibcz(:, 5), dm%fbcz(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      ref = - (ONE/scale)**2 * sin_wp(zc / scale + shift)
      err = dabs(fgzc(k) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'z-2ndder-c2c ', k, zc, ref, fgzc(k), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(3)) 
    write(wrt_unit, *) 'z-2ndder-c2c ', dm%nc(3), err_Linf, err_L2

! x direction, xp
    do i = 1, dm%np(1)
      xp = dm%h(1) * real(i - 1, WP)
      fxp(i) = sin_wp ( xp / scale + shift)
    end do
! x: p2p
    call Get_x_2nd_derivative_P2P_1D (fxp, fgxp, dm, dm%ibcx(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do i = 1, dm%np(1)
      xp = dm%h(1) * real(i - 1, WP)
      ref = - (ONE/scale)**2 * sin_wp(xp / scale + shift)
      err = dabs(fgxp(i) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'x-2ndder-p2p ', i, xp, ref, fgxp(i), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(1)) 
    write(wrt_unit, *) 'x-2ndder-p2p ', dm%np(1), err_Linf, err_L2


! y direction, yp
    do j = 1, dm%np(2)
      yp = dm%yp(j)
      fyp(j) = sin_wp ( yp / scale + shift)
    end do
! y: p2p
    call Get_y_2nd_derivative_P2P_1D (fyp, fgyp, dm, dm%ibcy(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do j = 1, dm%np(2)
      yp = dm%yp(j)
      ref = - (ONE/scale)**2 * sin_wp(yp / scale + shift)
      err = dabs(fgyp(j) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'y-2ndder-p2p ', j, yp, ref, fgyp(j), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(2)) 
    write(wrt_unit, *) 'y-2ndder-p2p ', dm%np(2), err_Linf, err_L2


! z direction, zp
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      fzp(k) = sin_wp ( zp / scale + shift)
    end do
! z: p2p
    call Get_z_2nd_derivative_P2P_1D (fzp, fgzp, dm, dm%ibcz(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      ref = - (ONE/scale)**2 * sin_wp(zp / scale + shift)
      err = dabs(fgzp(k) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      !write(*,'(A,1I5.1,4ES13.5)') 'z-2ndder-p2p ', k, zp, ref, fgzp(k), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(3)) 
    write(wrt_unit, *) 'z-2ndder-p2p ', dm%np(3), err_Linf, err_L2
    close(wrt_unit)

    return 
  end subroutine
end module

module burgers_eq_mod
  use precision_mod
  use parameters_constant_mod

  integer, parameter :: ICASE_HEATEQ = 12
  integer, parameter :: ICASE_INVSD_BURGERS = 13
  real(WP),parameter :: alpha = ONE, beta = ZERO

  ! udf variables
  integer :: icase = 12 ! which case
  integer :: idir = 1   ! which direction to test!

  private :: Compute_burgers_rhs
  private :: Validate_burgers_error
  public  :: Initialize_burgers_flow
  public  :: Solve_burgers_eq_iteration
  public  :: Plot_burgers_profile

contains
  subroutine  Initialize_burgers_flow(dm, ux, uy, uz, p)
    use udf_type_mod, only : t_domain, t_flow
    use math_mod, only : sin_wp
    use parameters_constant_mod!, only : HALF, ZERO, SIXTEEN, TWO
    use input_general_mod
    implicit none

    type(t_domain), intent(inout)   :: dm
    real(WP),       intent(inout) :: ux(:, :, :), &
                                     uy(:, :, :), &
                                     uz(:, :, :), &
                                     p (:, :, :)

    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    integer(4) :: i, j, k
    
    ux = ZERO
    uy = ZERO
    uz = ZERO
    p  = ZERO

    dm%icase = icase

  !==============================================================================
  ! example 1 : input sin(x)
  ! example 2 : input alpha * x + beta for inviscid Burgers' equation
  !===============================================================================
    if(icase == ICASE_HEATEQ .or. icase == ICASE_BURGERS) then
      ! bc
      dm%ibcx(:,:) = IBC_PERIODIC
      dm%ibcy(:,:) = IBC_PERIODIC
      dm%ibcz(:,:) = IBC_PERIODIC
      do i = 1, dm%np(idir)
        xp = dm%h(idir) * real(i - 1, WP)
        if(idir == 1) ux(i, :, :) =  sin_wp ( PI * xp )
        if(idir == 2) uy(:, i, :) =  sin_wp ( PI * xp )
        if(idir == 3) uz(:, :, i) =  sin_wp ( PI * xp )
      end do 
    else if (icase == ICASE_INVSD_BURGERS) then
      
      do i = 1, dm%np(idir)
        xp = dm%h(idir) * real(i - 1, WP)
        if(idir == 1) ux(i, :, :) =  alpha * xp + beta
        if(idir == 2) uy(:, i, :) =  alpha * xp + beta
        if(idir == 3) uz(:, :, i) =  alpha * xp + beta
      end do 
      if(idir == 1) then
        dm%ibcx(:,:) = IBC_DIRICHLET
        dm%fbcx(1, idir) = beta / (ONE)
        dm%fbcx(2, idir) = (alpha * dm%lxx + beta) / (ONE)
      else if(idir == 2) then
        dm%ibcy(:,:) = IBC_DIRICHLET
        dm%fbcy(1, idir) = beta / (ONE)
        dm%fbcy(2, idir) = (alpha * dm%lyt + beta) / (ONE)
      else if(idir == 3) then
        dm%ibcz(:,:) = IBC_DIRICHLET
        dm%fbcz(1, idir) = beta / (ONE)
        dm%fbcz(2, idir) = (alpha * dm%lzz + beta) / (ONE)
      else
      end if 
    else
    end if
    
    return
  end subroutine Initialize_burgers_flow

  subroutine Compute_burgers_rhs(fl, dm, isub)
    use udf_type_mod
    use parameters_constant_mod
    use operations
    use input_general_mod
    use boundary_conditions_mod
    use decomp_2d
    implicit none

    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(in   ) :: dm
    integer(4),     intent(in   ) :: isub  

    real(WP) :: fbc(2)
    integer :: i


    ! natural position as in staggered storage
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: mx_rhs ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: rhsx_dummy
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: my_rhs ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: rhsy_dummy
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: mz_rhs ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: rhsz_dummy
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: qy_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: my_rhs_ypencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: qz_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: qz_zpencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: mz_rhs_zpencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: mz_rhs_ypencil


    if(idir == 1) then
! xpencil 
      fl%mx_rhs = ZERO

      ! for x-mom convection term : d(qx * qx)/dx at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_INVSD_BURGERS) then
        do i = 1, 2
          fbc(i) = dm%ibcx(i, 1) * dm%ibcx(i, 1)
        end do
        call Get_x_1st_derivative_P2P_3D(-fl%qx * fl%qx * HALF, mx_rhs, dm, dm%ibcx(:, 1), fbc(:))
        fl%mx_rhs = fl%mx_rhs + mx_rhs
      end if

      ! for x-mom diffusion term , \mu * Ljj(ux) at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_HEATEQ) then
        call Get_x_2nd_derivative_P2P_3D( fl%qx, mx_rhs, dm, dm%ibcx(:, 1) )
        fl%mx_rhs = fl%mx_rhs + fl%rre * mx_rhs
      end if

      rhsx_dummy(:, :, :) = fl%mx_rhs(:, :, :)
      fl%mx_rhs(:, :, :) = dm%tGamma(isub) * fl%mx_rhs(:, :, :) + &
                           dm%tZeta (isub) * fl%mx_rhs0(:, :, :)
      fl%mx_rhs0(:, :, :) = rhsx_dummy(:, :, :)

      fl%qx(:, :, :) = fl%qx(:, :, :) + dm%dt * fl%mx_rhs(:, :, :)
      
    else if (idir == 2) then
! y pencil
      fl%my_rhs = ZERO
      my_rhs =  ZERO
      my_rhs_ypencil = ZERO

      call transpose_x_to_y (fl%qy,  qy_ypencil, dm%dcpc)     
      
      ! for y-mom convection term : d(qy * qy)/dy at (i, j', k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_INVSD_BURGERS) then
        do i = 1, 2
          fbc(i) = dm%ibcy(i, 2) * dm%ibcy(i, 2)
        end do
        call Get_y_1st_derivative_P2P_3D(-qy_ypencil * qy_ypencil * HALF, my_rhs_ypencil, dm, dm%ibcy(:, 2), fbc(:))
        call transpose_y_to_x (my_rhs_ypencil,  my_rhs)     
        fl%my_rhs = fl%my_rhs + my_rhs
      end if
      ! for x-mom diffusion term , \mu * Ljj(ux) at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_HEATEQ) then
        call Get_y_2nd_derivative_P2P_3D(qy_ypencil, my_rhs_ypencil, dm, dm%ibcy(:, 2) )
        call transpose_y_to_x (my_rhs_ypencil,  my_rhs)     
        fl%my_rhs = fl%my_rhs + fl%rre * my_rhs
      end if

      rhsy_dummy(:, :, :) = fl%my_rhs(:, :, :)
      fl%my_rhs(:, :, :)  = dm%tGamma(isub) * fl%my_rhs(:, :, :) + &
                            dm%tZeta (isub) * fl%my_rhs0(:, :, :)
      fl%my_rhs0(:, :, :) = rhsy_dummy(:, :, :)

      fl%qy(:, :, :) = fl%qy(:, :, :) + dm%dt * fl%my_rhs(:, :, :)

    else if (idir == 3) then
! z pencil
      call transpose_x_to_y (fl%qz,       qz_ypencil, dm%dccp)     
      call transpose_y_to_z (qz_ypencil,  qz_zpencil, dm%dccp)     

      fl%mz_rhs = ZERO
      mz_rhs =  ZERO
      mz_rhs_zpencil = ZERO

      ! for x-mom convection term : d(qx * qx)/dx at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_INVSD_BURGERS) then
        do i = 1, 2
          fbc(i) = dm%ibcz(i, 3) * dm%ibcz(i, 3)
        end do
        call Get_z_1st_derivative_P2P_3D(-qz_zpencil * qz_zpencil * HALF, mz_rhs_zpencil, dm, dm%ibcz(:, 3), fbc(:))
        call transpose_z_to_y (mz_rhs_zpencil,  mz_rhs_ypencil, dm%dccp)  
        call transpose_y_to_x (mz_rhs_ypencil,  mz_rhs,         dm%dccp)     
        fl%mz_rhs = fl%mz_rhs + mz_rhs
      end if
      ! for x-mom diffusion term , \mu * Ljj(ux) at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_HEATEQ) then
        call Get_z_2nd_derivative_P2P_3D( qz_zpencil, mz_rhs_zpencil, dm, dm%ibcz(:, 3))
        call transpose_z_to_y (mz_rhs_zpencil,  mz_rhs_ypencil, dm%dccp)  
        call transpose_y_to_x (mz_rhs_ypencil,  mz_rhs,         dm%dccp)    
        fl%mz_rhs = fl%mz_rhs + fl%rre * mz_rhs
      end if

      rhsz_dummy(:, :, :) = fl%mz_rhs(:, :, :)
      fl%mz_rhs(:, :, :) = dm%tGamma(isub) * fl%mz_rhs(:, :, :) + &
                           dm%tZeta (isub) * fl%mz_rhs0(:, :, :)
      fl%mz_rhs0(:, :, :) = rhsz_dummy(:, :, :)

      fl%qz(:, :, :) = fl%qz(:, :, :) + dm%dt * fl%mz_rhs(:, :, :)
    else
    end if

    ! apply bc
    call Apply_BC_velocity (dm, fl%qx, fl%qy, fl%qz)

    if(icase == ICASE_INVSD_BURGERS) then 
      if(idir == 1) then
        if (dm%dpcc%xst(1) == 1 )        fl%qx(1,        :, :) = beta / (alpha * fl%time + ONE)
        if (dm%dpcc%xen(1) == dm%np(1) ) fl%qx(dm%np(1), :, :) = (alpha * dm%lxx + beta) / (alpha * fl%time + ONE)
      else if(idir == 2) then
        if (dm%dcpc%xst(2) == 1 )        fl%qy(:, 1,        :) = beta / (alpha * fl%time + ONE)
        if (dm%dcpc%xen(2) == dm%np(2) ) fl%qy(:, dm%np(2), :) = (alpha * dm%lyt + beta) / (alpha * fl%time + ONE)
      else if(idir == 3) then
        if (dm%dccp%xst(3) == 1 )        fl%qz(:, :, 1      ) = beta / (alpha * fl%time + ONE)
        if (dm%dccp%xen(3) == dm%np(3) ) fl%qz(:, :, dm%np(3)) = (alpha * dm%lzz + beta) / (alpha * fl%time + ONE)
      else
      end if 
    end if

    return
  end subroutine Compute_burgers_rhs
!===============================================================================
  subroutine Validate_burgers_error(fl, dm)
    use udf_type_mod,            only : t_flow, t_domain
    use parameters_constant_mod
    use operations
    use math_mod
    use input_general_mod
    use mpi_mod
    implicit none

    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(in   ) :: dm
    integer :: i, j, k
    real(WP) :: xp, ux, uerr, uerr2, uerrmax, wavenum
    real(WP) :: uerr2_work, uerrmax_work
    real(WP) :: dd
    integer :: nx, ny, nz
    integer :: wrt_unit
    character( len = 128) :: filename
    logical :: file_exists = .FALSE.
    integer :: nsz

    
    uerr2 = ZERO
    uerrmax = ZERO

    dd = dm%h(idir)
    if(idir == 1) then
      wavenum = TWO * PI / dm%lxx
      nx = dm%dpcc%xsz(1)
      ny = dm%dpcc%xsz(2)
      nz = dm%dpcc%xsz(3)
    else if (idir == 2) then
      wavenum = TWO * PI / dm%lyt
      nx = dm%dcpc%xsz(1)
      ny = dm%dcpc%xsz(2)
      nz = dm%dcpc%xsz(3)
    else if (idir == 3) then
      wavenum = TWO * PI / dm%lzz
      nx = dm%dccp%xsz(1)
      ny = dm%dccp%xsz(2)
      nz = dm%dccp%xsz(3)
    else
    end if
    nsz =  nx * ny * nz

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if(idir == 1) xp = dd * real(i - 1, WP)
          if(idir == 2) xp = dd * real(j - 1, WP)
          if(idir == 3) xp = dd * real(k - 1, WP)
          if(icase == ICASE_BURGERS) then
            ux = sin_wp ( PI * xp ) * exp(- TWO * fl%rre * fl%time * wavenum * wavenum)
          else if(icase == ICASE_HEATEQ) then
            ux = sin_wp ( PI * xp ) * exp(- TWO * fl%rre * fl%time * wavenum * wavenum) ! check
          else if(icase == ICASE_INVSD_BURGERS) then
            ux = (alpha * xp + beta )/(alpha * fl%time + ONE) ! check
          else
          end if
          if(idir == 1) uerr = fl%qx(i, j, k) - ux
          if(idir == 2) uerr = fl%qy(i, j, k) - ux
          if(idir == 3) uerr = fl%qz(i, j, k) - ux

          uerr2 = uerr2 + uerr**2
          if(dabs(uerr) > uerrmax) uerrmax = dabs(uerr)
          !if(k==d%nc(3)/2 .and. j == d%nc(2)/2) write(*,*) k, j, i, ux, fl%qx(i, j, k), uerr
        end do 
      end do
    end do

    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(uerr2,   uerr2_work,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    call mpi_allreduce(uerrmax, uerrmax_work, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
    uerr2_work = sqrt_wp(uerr2_work / real(nsz, WP) )
    uerr2_work = sqrt_wp(uerr2_work)

    if(nrank == 0 ) then
      filename = 'Validation_Burgers.dat'

      INQUIRE(FILE = trim(filename), exist = file_exists)

      if(.not.file_exists) then
        open(newunit = wrt_unit, file = trim(filename), action = "write", status = "new")
        write(wrt_unit, '(A)') 'Time, SD(uerr), Max(uerr)'
      else
        open(newunit = wrt_unit, file = trim(filename), action = "write", status = "old", position="append")
      end if
      ! data convert to cell centre data...

      write(wrt_unit, '(1F10.4, 2ES15.7)') fl%time, uerr2_work, uerrmax_work
      close(wrt_unit)
    end if

  end subroutine 
  !===============================================================================
  subroutine Plot_burgers_profile(fl, dm, iter)
    use udf_type_mod,            only : t_flow, t_domain
    use parameters_constant_mod
    use input_general_mod
    use operations
    use math_mod
    use typeconvert_mod
    use mpi_mod
    implicit none

    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(in   ) :: dm
    integer, intent(in) :: iter
    integer :: i, j, k
    real(WP) :: xp, ux, uerr, uerr2, uerrmax, wavenum
    real(WP) :: dd
    integer :: nx, ny, nz

    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: qy_ypencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: qz_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: qz_zpencil

    integer :: wrt_unit
    character( len = 128) :: filename
    logical :: file_exists = .FALSE.
    

    if(nrank == 0) then

    filename = 'Plot_Burgers_profile'//trim(int2str(iter))//'.dat'

    INQUIRE(FILE = trim(filename), exist = file_exists)

    if(.not.file_exists) then
      open(newunit = wrt_unit, file = trim(filename), action = "write", status = "new")
      write(wrt_unit, '(A)') 'x qx'
    else
      open(newunit = wrt_unit, file = trim(filename), action = "write", status = "old", position="append")
     end if
    ! data convert to cell centre data...
    end if

    wavenum = ONE
    uerr = ZERO
    uerrmax = ZERO

    dd = dm%h(idir)
    if(idir == 1) then
      do i = 1, dm%np(idir)
        if(nrank == 0) write(wrt_unit, '(1F10.4, 2ES15.7)') dd*real(i, WP), fl%qx(i, dm%nc(2)/2, dm%nc(3)/2)
      end do
    else if (idir == 2) then
      call transpose_x_to_y (fl%qy,       qy_ypencil, dm%dcpc)    
      do i = 1, dm%np(idir)
        if(nrank == 0) write(wrt_unit, '(1F10.4, 2ES15.7)') dd*real(i, WP), qy_ypencil(dm%nc(1)/2, i, dm%nc(3)/2)
      end do
    else if (idir == 3) then
      call transpose_x_to_y (fl%qz,       qz_ypencil, dm%dccp)    
      call transpose_y_to_z (qz_ypencil,  qz_zpencil, dm%dccp)    
      do i = 1, dm%np(idir)
        if(nrank == 0) write(wrt_unit, '(1F10.4, 2ES15.7)') dd*real(i, WP), qz_zpencil(dm%nc(1)/2, dm%nc(2)/2, i)
      end do
    else
    end if

    if(nrank == 0)close(wrt_unit)

  end subroutine 
!===============================================================================
  subroutine Solve_burgers_eq_iteration
    use parameters_constant_mod
    use mpi_mod
    use vars_df_mod
    use solver_tools_mod
    use thermo_info_mod
    use code_performance_mod
    use input_general_mod
    implicit none

    logical :: is_flow   = .false.
    logical :: is_thermo = .false.
    integer :: i
    integer :: iter, isub
    integer :: nrsttckpt
    integer :: niter
    real(wp)   :: rtmp

    
    nrsttckpt = HUGE(0)
    niter     = 0
    do i = 1, nxdomain
      if( flow(i)%nrsttckpt < nrsttckpt) nrsttckpt = flow(i)%nrsttckpt
      if( flow(i)%nIterFlowEnd > niter)  niter     = flow(i)%nIterFlowEnd
      if( is_any_energyeq) then
        if (thermo(i)%nIterThermoEnd > niter) niter = thermo(i)%nIterThermoEnd
      end if
    end do

    do iter = nrsttckpt + 1, niter
      call Call_cpu_time(CPU_TIME_ITER_START, nrsttckpt, niter, iter)
      do i = 1, nxdomain
!===============================================================================
!      setting up 1/re, 1/re/prt, gravity, etc
!===============================================================================
        call Update_Re(iter, flow(i))
        if(domain(i)%ithermo == 1) &
        call Update_PrGr(flow(i), thermo(i))
!===============================================================================
!      setting up flow solver
!===============================================================================
        if ( (iter >= flow(i)%nIterFlowStart) .and. (iter <=flow(i)%nIterFlowEnd)) then
          is_flow = .true.
          flow(i)%time = flow(i)%time + domain(i)%dt
          !call Check_cfl_diffusion (domain(i)%h2r(:), flow(i)%rre, domain(i)%dt)
          !call Check_cfl_convection(flow(i)%qx, flow(i)%qy, flow(i)%qz, domain(i))
        end if
!===============================================================================
!     setting up thermo solver
!===============================================================================
        if(domain(i)%ithermo == 1) then
          if ( (iter >= thermo(i)%nIterThermoStart) .and. (iter <= thermo(i)%nIterThermoEnd)) then
            is_thermo = .true.
            thermo(i)%time = thermo(i)%time  + domain(i)%dt
          end if
        end if
!===============================================================================
!     main solver
!===============================================================================
        do isub = 1, domain(i)%nsubitr
          !if(is_thermo) call Solve_energy_eq  (flow(i), thermo(i), domain(i), isub)
          !if(is_flow)   call Solve_momentum_eq(flow(i), domain(i), isub)
          call Compute_burgers_rhs(flow(i), domain(i), isub)
        end do
  
        call Validate_burgers_error (flow(i), domain(i))
        if( MOD(iter, domain(i)%nvisu) == 0 ) call Plot_burgers_profile(flow(i), domain(i), iter)

      end do
      
      

    end do
  
  
    call Call_cpu_time(CPU_TIME_CODE_END, nrsttckpt, niter)
    call Finalise_mpi()
    stop
    return
  end subroutine Solve_burgers_eq_iteration

end module


!===============================================================================
!===============================================================================
!> \brief In-code independent test code for algorithms and schemes
!>
!> This subroutine is only called in the main program for testing.
!> Please select the test options which you are interested in.
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
subroutine Test_algorithms()
  use vars_df_mod
  use test_operations_mod
  use burgers_eq_mod
  implicit none

  logical :: is_TDMA = .false.
  logical :: is_operations = .false.
  logical :: is_burgers = .true.

  if(is_TDMA) then
    call Test_TDMA_cyclic
    call Test_TDMA_noncyclic
  end if
  
  if(is_operations) then
    call Test_interpolation (domain(1))
    call Test_1st_derivative(domain(1))
    call Test_2nd_derivative(domain(1))
  end if

  if (is_burgers) then
    call Solve_burgers_eq_iteration
  end if

  return 
end subroutine 