module test_algorithms_mod
  use operations

  public   :: Test_schemes
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
subroutine Test_schemes()
  use vars_df_mod
  implicit none

  ! Please use below information for input file
  ! x = 0, 2pi
  ! y = 0, 2pi
  ! z = 0, 2pi

  !call Test_TDMA_cyclic
  !call Test_TDMA_noncyclic
  call Test_interpolation (domain(1))
  call Test_1st_derivative(domain(1))
  call Test_2nd_derivative(domain(1))
  return 
end subroutine 

!===============================================================================
!===============================================================================
!> \brief To test this subroutine for mid-point interpolation.
!>
!> This subroutine is called in \ref Test_schemes. Define the logicals to choose
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
      write(*,'(A,1I5.1,4ES13.5)') 'x-interp-c2p ', i, xp, ref, fgxp(i), err !test
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
      write(*,'(A,1I5.1,4ES13.5)') 'y-interp-c2p ', j, yp, ref, fgyp(j), err !test
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
      write(*,'(A,1I5.1,4ES13.5)') 'z-interp-c2p ', k, zp, ref, fgzp(k), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(3)) 
    write(wrt_unit, *) 'z-interp-c2p ', dm%np(3), err_Linf, err_L2

! x direction, xp
    do i = 1, dm%np(1)
      xp = dm%h(1) * real(i - 1, WP)
      fxp(i) = sin_wp ( xp / scale + shift)
    end do
! x: p2c
    call Get_x_midp_P2C_1D (fxp, fgxc, dm, dm%ibcx(:, 5), dm%fbcx(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do i = 1, dm%nc(1)
      xc = dm%h(1) * (real(i - 1, WP) + HALF)
      ref = sin_wp(xc / scale + shift)
      err = dabs(fgxc(i) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      write(*,'(A,1I5.1,4ES13.5)') 'x-interp-p2c ', i, xc, ref, fgxc(i), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(1)) 
    write(wrt_unit, *) 'x-interp-p2c ', dm%nc(1), err_Linf, err_L2

! y direction, yp
    do j = 1, dm%np(2)
      yp = dm%yp(j)
      fyp(j) = sin_wp ( yp / scale + shift)
    end do
! y: p2c
    call Get_y_midp_P2C_1D (fyp, fgyc, dm, dm%ibcy(:, 5), dm%fbcy(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do j = 1, dm%nc(2)
      yc = dm%yc(j)
      ref = sin_wp(yc / scale + shift)
      err = dabs(fgyc(j) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      write(*,'(A,1I5.1,4ES13.5)') 'y-interp-p2c ', j, yc, ref, fgyc(j), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(2)) 
    write(wrt_unit, *) 'y-interp-p2c ', dm%nc(2), err_Linf, err_L2


! z direction, zp
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      fzp(k) = sin_wp ( zp / scale + shift)
    end do
! z: p2c
    call Get_z_midp_P2C_1D (fzp, fgzc, dm, dm%ibcz(:, 5), dm%fbcz(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      ref = sin_wp(zc / scale + shift)
      err = dabs(fgzc(k) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      write(*,'(A,1I5.1,4ES13.5)') 'z-interp-p2c ', k, zc, ref, fgzc(k), err !test
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
!> This subroutine is called in \ref Test_schemes. Define the logicals to choose
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
      write(*,'(A,1I5.1,4ES13.5)') 'x-1stder-c2c ', i, xc, ref, fgxc(i), err !test
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
      write(*,'(A,1I5.1,4ES13.5)') 'y-1stder-c2c ', j, yc, ref, fgyc(j), err !test
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
      write(*,'(A,1I5.1,4ES13.5)') 'z-1stder-c2c ', k, zc, ref, fgzc(k), err !test
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
      write(*,'(A,1I5.1,4ES13.5)') 'x-1stder-p2p', i, xp, ref, fgxp(i), err !test
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
      write(*,'(A,1I5.1,4ES13.5)') 'y-1stder-p2p ', j, yp, ref, fgyp(j), err !test
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
      write(*,'(A,1I5.1,4ES13.5)') 'z-1stder-p2p ', k, zp, ref, fgzp(k), err !test
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
      write(*,'(A,1I5.1,4ES13.5)') 'x-1stder-c2p ', i, xp, ref, fgxp(i), err !test
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
      write(*,'(A,1I5.1,4ES13.5)') 'y-1stder-c2p ', j, yp, ref, fgyp(j), err !test
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
      write(*,'(A,1I5.1,4ES13.5)') 'z-1stder-c2p ', k, zp, ref, fgzp(k), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(3)) 
    write(wrt_unit, *) 'z-1stder-c2p ', dm%np(3), err_Linf, err_L2

! x: p2c
    call Get_x_1st_derivative_P2C_1D (fxp, fgxc, dm, dm%ibcx(:, 5), dm%fbcx(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do i = 1, dm%nc(1)
      xc = dm%h(1) * (real(i - 1, WP) + HALF)
      ref = ONE/scale * cos_wp(xc / scale + shift)
      err = dabs(fgxc(i) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      write(*,'(A,1I5.1,4ES13.5)') 'x-1stder-p2c ', i, xc, ref, fgxc(i), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(1)) 
    write(wrt_unit, *) 'x-1stder-p2c ', dm%nc(1), err_Linf, err_L2

! y: p2c
    call Get_y_1st_derivative_P2C_1D (fyp, fgyc, dm, dm%ibcy(:, 5), dm%fbcy(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do j = 1, dm%nc(2)
      yc = dm%yc(j)
      ref = ONE/scale * cos_wp(yc / scale + shift)
      err = dabs(fgyc(j) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      write(*,'(A,1I5.1,4ES13.5)') 'y-1stder-p2c ', j, yc, ref, fgyc(j), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(2)) 
    write(wrt_unit, *) 'y-1stder-p2c ', dm%nc(2), err_Linf, err_L2

! z: p2c
    call Get_z_1st_derivative_P2C_1D (fzp, fgzc, dm, dm%ibcz(:, 5), dm%fbcz(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      ref = ONE/scale * cos_wp(zc / scale + shift)
      err = dabs(fgzc(k) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      write(*,'(A,1I5.1,4ES13.5)') 'z-1stder-p2c ', k, zc, ref, fgzc(k), err !test
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
!> This subroutine is called in \ref Test_schemes. Define the logicals to choose
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
      write(*,'(A,1I5.1,4ES13.5)') 'x-2ndder-c2c ', i, xc, ref, fgxc(i), err !test
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
     ! write(*,'(A,1I5.1,4ES13.5)') 'y-2ndder-c2c ', j, yc, ref, fgyc(j), err !test
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
      write(*,'(A,1I5.1,4ES13.5)') 'z-2ndder-c2c ', k, zc, ref, fgzc(k), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%nc(3)) 
    write(wrt_unit, *) 'z-2ndder-c2c ', dm%nc(3), err_Linf, err_L2

! x direction, xp
    do i = 1, dm%np(1)
      xp = dm%h(1) * real(i - 1, WP)
      fxp(i) = sin_wp ( xp / scale + shift)
    end do
! x: p2p
    call Get_x_2nd_derivative_P2P_1D (fxp, fgxp, dm, dm%ibcx(:, 5), dm%fbcx(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do i = 1, dm%np(1)
      xp = dm%h(1) * real(i - 1, WP)
      ref = - (ONE/scale)**2 * sin_wp(xp / scale + shift)
      err = dabs(fgxp(i) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      write(*,'(A,1I5.1,4ES13.5)') 'x-2ndder-p2p ', i, xp, ref, fgxp(i), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(1)) 
    write(wrt_unit, *) 'x-2ndder-p2p ', dm%np(1), err_Linf, err_L2


! y direction, yp
    do j = 1, dm%np(2)
      yp = dm%yp(j)
      fyp(j) = sin_wp ( yp / scale + shift)
    end do
! y: p2p
    call Get_y_2nd_derivative_P2P_1D (fyp, fgyp, dm, dm%ibcy(:, 5), dm%fbcy(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do j = 1, dm%np(2)
      yp = dm%yp(j)
      ref = - (ONE/scale)**2 * sin_wp(yp / scale + shift)
      err = dabs(fgyp(j) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      write(*,'(A,1I5.1,4ES13.5)') 'y-2ndder-p2p ', j, yp, ref, fgyp(j), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(2)) 
    write(wrt_unit, *) 'y-2ndder-p2p ', dm%np(2), err_Linf, err_L2


! z direction, zp
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      fzp(k) = sin_wp ( zp / scale + shift)
    end do
! z: p2p
    call Get_z_2nd_derivative_P2P_1D (fzp, fgzp, dm, dm%ibcz(:, 5), dm%fbcz(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      ref = - (ONE/scale)**2 * sin_wp(zp / scale + shift)
      err = dabs(fgzp(k) - ref)
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
      write(*,'(A,1I5.1,4ES13.5)') 'z-2ndder-p2p ', k, zp, ref, fgzp(k), err !test
    end do
    err_L2 = sqrt_wp(err_L2 / dm%np(3)) 
    write(wrt_unit, *) 'z-2ndder-p2p ', dm%np(3), err_Linf, err_L2
    close(wrt_unit)

    return 
  end subroutine
end module