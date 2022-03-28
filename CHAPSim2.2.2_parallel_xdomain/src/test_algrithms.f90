module test_algorithms_mod
  use operations

  public   :: Test_algorithms
  private  :: Test_interpolation
  private  :: Test_1st_derivative
  private  :: Test_2nd_derivative
  
contains
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
  implicit none

  logical is_TDMA = .false.
  logical is_operations = .false.
  logical is_burgers = .true.

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
  end if

  call 
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

!===============================================================================
!===============================================================================

module burgers_eq_mod

  private :: Compute_burgers_rhs
  private :: Validate_burgers_error
  public  :: Solve_burgers_eq_iteration
  public  :: Plot_burgers_profile

contains
  subroutine Compute_burgers_rhs(fl, dm, isub)
    use udf_type_mod
    use parameters_constant_mod
    use operations
    use input_general_mod
    use boundary_conditions_mod
    use nvtx
    implicit none

    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(in   ) :: m
    integer(4),     intent(in   ) :: isub  

    real(WP),parameter :: alpha = ONE, beta = ZERO


    ! natural position as in staggered storage
    real(WP), dimension( dm%np(1), dm%nc(2), dm%nc(3) ) :: m1_rhs ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( dm%np(1), dm%nc(2), dm%nc(3) ) :: rhs1_dummy
    real(WP), dimension( dm%nc(1), dm%np(2), dm%nc(3) ) :: m2_rhs ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( dm%nc(1), dm%np(2), dm%nc(3) ) :: rhs2_dummy
    real(WP), dimension( dm%nc(1), dm%nc(2), dm%np(3) ) :: m3_rhs ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( dm%nc(1), dm%nc(2), dm%np(3) ) :: rhs3_dummy

    if(idir == 1) then
      f%m1_rhs = ZERO
      call nvtxStartRange("Get_x_1st_derivative_P2P_3dArray")
      ! for x-mom convection term : d(qx * qx)/dx at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_INVSD_BURGERS) then
        call Get_x_1st_derivative_P2P_3dArray( d, -f%qx * f%qx * HALF, m1_rhs )
        f%m1_rhs = f%m1_rhs + m1_rhs
      end if
      call nvtxEndRange

      call nvtxStartRange("Get_x_2nd_derivative_P2P_3dArray")
      ! for x-mom diffusion term , \mu * Ljj(ux) at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_HEATEQ) then
        call Get_x_2nd_derivative_P2P_3dArray( d, f%qx, m1_rhs )
        f%m1_rhs = f%m1_rhs + f%rre * m1_rhs
      end if
      call nvtxEndRange

      call nvtxStartRange("RK")
      rhs1_dummy(:, :, :) = f%m1_rhs(:, :, :)
      f%m1_rhs(:, :, :) = tGamma(isub) * f%m1_rhs(:, :, :) + &
                          tZeta (isub) * f%m1_rhs0(:, :, :)
      f%m1_rhs0(:, :, :) = rhs1_dummy(:, :, :)

      f%qx(:, :, :) = f%qx(:, :, :) + dt * f%m1_rhs(:, :, :)
      call nvtxEndRange
    else if (idir == 2) then
      f%m2_rhs = ZERO
      ! for x-mom convection term : d(qx * qx)/dx at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_INVSD_BURGERS) then
        call Get_y_1st_derivative_P2P_3dArray( d, -f%qy * f%qy * HALF, m2_rhs )
        f%m2_rhs = f%m2_rhs + m2_rhs
      end if
      ! for x-mom diffusion term , \mu * Ljj(ux) at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_HEATEQ) then
        call Get_y_2nd_derivative_P2P_3dArray( d, f%qy, m2_rhs )
        f%m2_rhs = f%m2_rhs + f%rre * m2_rhs
      end if

      rhs2_dummy(:, :, :) = f%m2_rhs(:, :, :)
      f%m2_rhs(:, :, :) = tGamma(isub) * f%m2_rhs(:, :, :) + &
                          tZeta (isub) * f%m2_rhs0(:, :, :)
      f%m2_rhs0(:, :, :) = rhs2_dummy(:, :, :)

      f%qy(:, :, :) = f%qy(:, :, :) + dt * f%m2_rhs(:, :, :)
    else if (idir == 3) then
      f%m3_rhs = ZERO
      ! for x-mom convection term : d(qx * qx)/dx at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_INVSD_BURGERS) then
        call Get_z_1st_derivative_P2P_3dArray( d, -f%qz * f%qz * HALF, m3_rhs )
        f%m3_rhs = f%m3_rhs + m3_rhs
      end if
      ! for x-mom diffusion term , \mu * Ljj(ux) at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_HEATEQ) then
        call Get_z_2nd_derivative_P2P_3dArray( d, f%qz, m3_rhs )
        f%m3_rhs = f%m3_rhs + f%rre * m3_rhs
      end if

      rhs3_dummy(:, :, :) = f%m3_rhs(:, :, :)
      f%m3_rhs(:, :, :) = tGamma(isub) * f%m3_rhs(:, :, :) + &
                          tZeta (isub) * f%m3_rhs0(:, :, :)
      f%m3_rhs0(:, :, :) = rhs3_dummy(:, :, :)

      f%qz(:, :, :) = f%qz(:, :, :) + dt * f%m3_rhs(:, :, :)
    else
    end if

    call Apply_BC_velocity (f%qx, f%qy, f%qz, d)

    if(icase == ICASE_INVSD_BURGERS) then 
      if(idir == 1) then
        f%qx(1,       :, :) = beta / (alpha * f%time + ONE)
        f%qx(d%np(1), :, :) = (alpha * lxx + beta) / (alpha * f%time + ONE)
      else if(idir == 2) then
        f%qy(:, 1,       :) = beta / (alpha * f%time + ONE)
        f%qy(:, d%np(2), :) = (alpha * lyt + beta) / (alpha * f%time + ONE)
      else if(idir == 3) then
        f%qz(:, :, 1      ) = beta / (alpha * f%time + ONE)
        f%qz(:, :, d%np(3)) = (alpha * lzz + beta) / (alpha * f%time + ONE)
      else
      end if 
    end if

    return
  end subroutine Compute_burgers_rhs
!===============================================================================
  subroutine Validate_burgers_error(f, d)
    use udf_type_mod,            only : t_flow, t_domain
    use parameters_constant_mod
    use operations
    use math_mod
    use input_general_mod
    implicit none

    type(t_flow),   intent(inout) :: f
    type(t_domain), intent(in   ) :: d
    integer :: i, j, k
    real(WP) :: xp, ux, uerr, uerr2, uerrmax, wavenum
    real(WP) :: dd
    integer :: nx, ny, nz

    integer :: output_unit
    character( len = 128) :: filename
    logical :: file_exists = .FALSE.
    real(WP),parameter :: alpha = ONE, beta = ZERO

    

    filename = 'Validation_Burgers.dat'

    INQUIRE(FILE = trim(filename), exist = file_exists)

    if(.not.file_exists) then
      open(newunit = output_unit, file = trim(filename), action = "write", status = "new")
      write(output_unit, '(A)') 'Time, SD(uerr), Max(uerr)'
    else
      open(newunit = output_unit, file = trim(filename), action = "write", status = "old", position="append")
     end if
    ! data convert to cell centre data...


    wavenum = TWO * PI / lxx
    uerr2 = ZERO
    uerrmax = ZERO

    dd = d%h(idir)
    if(idir == 1) then
      nx = d%np(1)
      ny = d%nc(2)
      nz = d%nc(3)
    else if (idir == 2) then
      nx = d%nc(1)
      ny = d%np(2)
      nz = d%nc(3)
    else if (idir == 3) then
      nx = d%nc(1)
      ny = d%np(2)
      nz = d%nc(3)
    else
    end if

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if(idir == 1) xp = dd * real(i - 1, WP)
          if(idir == 2) xp = dd * real(j - 1, WP)
          if(idir == 3) xp = dd * real(k - 1, WP)
          if(icase == ICASE_BURGERS) then
            ux = sin_wp ( PI * xp ) * exp(- TWO * f%rre * f%time * wavenum * wavenum)
          else if(icase == ICASE_HEATEQ) then
            ux = sin_wp ( PI * xp ) * exp(- TWO * f%rre * f%time * wavenum * wavenum) ! check
          else if(icase == ICASE_INVSD_BURGERS) then
            ux = (alpha * xp + beta )/(alpha * f%time + ONE) ! check
          else
          end if
          if(idir == 1) uerr = f%qx(i, j, k) - ux
          if(idir == 2) uerr = f%qy(i, j, k) - ux
          if(idir == 3) uerr = f%qz(i, j, k) - ux

          uerr2 = uerr2 + uerr**2
          if(dabs(uerr) > uerrmax) uerrmax = dabs(uerr)
          !if(k==d%nc(3)/2 .and. j == d%nc(2)/2) write(*,*) k, j, i, ux, f%qx(i, j, k), uerr
        end do 
      end do
    end do
    uerr2 = uerr2 / real(d%np(1), wp) / real(d%nc(2), wp) / real(d%nc(3), wp)
    uerr2 = sqrt_wp(uerr2)

    write(output_unit, '(1F10.4, 2ES15.7)') f%time, uerr2, uerrmax
    close(output_unit)

  end subroutine 
  !===============================================================================
  subroutine Plot_burgers_profile(f, d, iter)
    use udf_type_mod,            only : t_flow, t_domain
    use parameters_constant_mod
    use input_general_mod
    use operations
    use math_mod
    use typeconvert_mod
    implicit none

    type(t_flow),   intent(inout) :: f
    type(t_domain), intent(in   ) :: d
    integer, intent(in) :: iter
    integer :: i, j, k
    real(WP) :: xp, ux, uerr, uerr2, uerrmax, wavenum
    real(WP) :: dd
    integer :: nx, ny, nz

    integer :: output_unit
    character( len = 128) :: filename
    logical :: file_exists = .FALSE.
    

    filename = 'Plot_Burgers_profile'//trim(int2str(iter))//'.dat'

    INQUIRE(FILE = trim(filename), exist = file_exists)

    if(.not.file_exists) then
      open(newunit = output_unit, file = trim(filename), action = "write", status = "new")
      write(output_unit, '(A)') 'Time, SD(uerr), Max(uerr)'
    else
      open(newunit = output_unit, file = trim(filename), action = "write", status = "old", position="append")
     end if
    ! data convert to cell centre data...


    wavenum = ONE
    uerr = ZERO
    uerrmax = ZERO

    dd = d%h(idir)
    if(idir == 1) then
      do i = 1, d%np(idir)
        write(output_unit, '(1F10.4, 2ES15.7)') dd*real(i, WP), f%qx(i, d%nc(2)/2, d%nc(3)/2)
      end do
    else if (idir == 2) then
      do i = 1, d%np(idir)
        write(output_unit, '(1F10.4, 2ES15.7)') dd*real(i, WP), f%qy(d%nc(1)/2, i, d%nc(3)/2)
      end do
    else if (idir == 3) then
      do i = 1, d%np(idir)
        write(output_unit, '(1F10.4, 2ES15.7)') dd*real(i, WP), f%qz(d%nc(1)/2, d%nc(2)/2, i)
      end do
    else
    end if

    close(output_unit)

  end subroutine 
!===============================================================================
  subroutine Solve_burgers_eq_iteration
    use input_general_mod!,  only :  ithermo, nIterFlowEnd, nIterThermoEnd, &
                                    !nIterFlowStart, nIterThermoStart, &
                                    !tThermo, tFlow, nIterFlowEnd, nrsttckpt, &
                                    !dt, nsubitr, niter
    use type_vars_mod!,      only : flow, thermo, domain
    use flow_variables_mod!, only : Calculate_RePrGr, Check_maximum_velocity
    use eq_momentum_mod!,    only : Solve_momentum_eq
    use solver_tools_mod!,   only : Check_cfl_diffusion, Check_cfl_convection
    use continuity_eq_mod
    use poisson_mod
    use code_performance_mod
    use parameters_constant_mod
    use solver_tools_mod
    implicit none

    integer(4) :: iter, isub
    real(wp)   :: rtmp

    niter = nIterFlowEnd
    flow%time = tFlow 
    flow%rre = ONE / ren
    call Check_cfl_diffusion_1d(domain%h2r(IDIR), flow%rre)
    !call Plot_burgers_profile(flow, domain, 0)

    do iter = nrsttckpt + 1, niter
      call Call_cpu_time(CPU_TIME_ITER_START, nrsttckpt, niter, iter)
  !===============================================================================
  !     main solver
  !===============================================================================
      flow%time = flow%time + dt
      do isub = 1, nsubitr
        call Compute_burgers_rhs(flow, domain, isub)
      end do
  !===============================================================================
  !     validation
  !===============================================================================
      !call Validate_burgers_error (flow, domain)

      !if(MOD(iter, nvisu) == 0) call Plot_burgers_profile(flow, domain, iter)

      call Call_cpu_time(CPU_TIME_ITER_END, nrsttckpt, niter, iter)

    end do

    return
  end subroutine Solve_burgers_eq_iteration

end module
