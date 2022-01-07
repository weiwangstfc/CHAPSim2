
module burgers_eq_mod

  private :: Compute_burgers_rhs
  private :: Validate_burgers_error
  public  :: Solve_burgers_eq_iteration
  public  :: Plot_burgers_profile

contains
  subroutine Compute_burgers_rhs(f, d, isub)
    use udf_type_mod,            only : t_flow, t_domain
    use parameters_constant_mod
    use operations
    use input_general_mod
    use boundary_conditions_mod
    implicit none

    type(t_flow),   intent(inout) :: f
    type(t_domain), intent(in   ) :: d
    integer(4),     intent(in   ) :: isub  

    real(WP),parameter :: alpha = ONE, beta = ZERO


    ! natural position as in staggered storage
    real(WP), dimension( d%np(1), d%nc(2), d%nc(3) ) :: m1_rhs ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( d%np(1), d%nc(2), d%nc(3) ) :: rhs1_dummy
    real(WP), dimension( d%nc(1), d%np(2), d%nc(3) ) :: m2_rhs ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( d%nc(1), d%np(2), d%nc(3) ) :: rhs2_dummy
    real(WP), dimension( d%nc(1), d%nc(2), d%np(3) ) :: m3_rhs ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( d%nc(1), d%nc(2), d%np(3) ) :: rhs3_dummy

    if(idir == 1) then
      f%m1_rhs = ZERO
      ! for x-mom convection term : d(qx * qx)/dx at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_INVSD_BURGERS) then
        call Get_x_1st_derivative_P2P_3dArray( d, -f%qx * f%qx * HALF, m1_rhs )
        f%m1_rhs = f%m1_rhs + m1_rhs
      end if
      ! for x-mom diffusion term , \mu * Ljj(ux) at (i', j, k)
      if(icase == ICASE_BURGERS .or. icase == ICASE_HEATEQ) then
        call Get_x_2nd_derivative_P2P_3dArray( d, f%qx, m1_rhs )
        f%m1_rhs = f%m1_rhs + f%rre * m1_rhs
      end if

      rhs1_dummy(:, :, :) = f%m1_rhs(:, :, :)
      f%m1_rhs(:, :, :) = tGamma(isub) * f%m1_rhs(:, :, :) + &
                          tZeta (isub) * f%m1_rhs0(:, :, :)
      f%m1_rhs0(:, :, :) = rhs1_dummy(:, :, :)

      f%qx(:, :, :) = f%qx(:, :, :) + dt * f%m1_rhs(:, :, :)
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
    implicit none

    integer(4) :: iter, isub
    real(wp)   :: rtmp

    niter = nIterFlowEnd
    flow%time = tFlow 
    flow%rre = ONE / ren

    call Plot_burgers_profile(flow, domain, 0)

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
      call Validate_burgers_error (flow, domain)

      if(MOD(iter, nvisu) == 0) call Plot_burgers_profile(flow, domain, iter)

      call Call_cpu_time(CPU_TIME_ITER_END, nrsttckpt, niter, iter)

    end do

    return
  end subroutine Solve_burgers_eq_iteration

end module
