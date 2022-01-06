
module burgers_eq_mod

  private :: Compute_burgers_rhs
  private :: Validate_burgers_error
  public  :: Solve_burgers_eq_iteration

contains
  subroutine Compute_burgers_rhs(f, d, isub)
    use udf_type_mod,            only : t_flow, t_domain
    use parameters_constant_mod
    use operations
    use input_general_mod
    implicit none

    type(t_flow),   intent(inout) :: f
    type(t_domain), intent(in   ) :: d
    integer(4),     intent(in   ) :: isub  
  !-------------------------------------------------------------------------------
  ! common vars
  !-------------------------------------------------------------------------------
    ! natural position as in staggered storage
    real(WP), dimension( d%np(1), d%nc(2), d%nc(3) ) :: m1_rhs ! rhs for momentum-x at (xp, yc, zc)
    real(WP), dimension( d%np(1), d%nc(2), d%nc(3) ) :: rhs_dummy


    f%m1_rhs = ZERO
   ! for x-mom convection term (x-c1/3): d(qx * qx)/dx at (i', j, k)
    call Get_x_1st_derivative_P2P_3dArray( d, -f%qx * f%qx * HALF, m1_rhs )
    f%m1_rhs = f%m1_rhs + m1_rhs


    ! for x-mom diffusion term (x-v1/1), \mu * Ljj(ux) at (i', j, k)
    call Get_x_2nd_derivative_P2P_3dArray( d, f%qx, m1_rhs )
    f%m1_rhs = f%m1_rhs + f%rre * m1_rhs

  !-------------------------------------------------------------------------------
  ! to build up rhs in total, in all directions
  !_______________________________________________________________________________ 
    ! x-momentum
    ! add explicit terms
    rhs_dummy(:, :, :) = f%m1_rhs(:, :, :)
    f%m1_rhs(:, :, :) = tGamma(isub) * f%m1_rhs(:, :, :) + &
                        tZeta (isub) * f%m1_rhs0(:, :, :)
    f%m1_rhs0(:, :, :) = rhs_dummy(:, :, :)

    ! times the time step 
    f%m1_rhs(:, :, :) = dt *f%m1_rhs(:, :, :)

    f%qx(:, :, :) = f%qx(:, :, :) + f%m1_rhs(:, :, :)

    return
  end subroutine Compute_burgers_rhs
!===============================================================================
  subroutine Validate_burgers_error(f, d)
    use udf_type_mod,            only : t_flow, t_domain
    use parameters_constant_mod
    use operations
    use math_mod
    implicit none

    type(t_flow),   intent(inout) :: f
    type(t_domain), intent(in   ) :: d
    integer :: i, j, k
    real(WP) :: xp, ux, uerr, uerr2, uerrmax, wavenum

    integer :: output_unit
    character( len = 128) :: filename
    logical :: file_exists = .FALSE.
    

    filename = 'Validation_Burgers.dat'

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
    do k = 1, d%nc(3)
      do j = 1, d%nc(2)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          ux =  sin_wp ( xp ) * exp(- TWO * f%rre * f%time * wavenum * wavenum)
          uerr = f%qx(i, j, k) - ux
          uerr = uerr + uerr**2
          if(dabs(uerr) > uerrmax) uerrmax = dabs(uerr)
          !if(k==d%nc(3)/2 .and. j == d%nc(2)/2) write(*,*) k, j, i, ux, f%qx(i, j, k), uerr
        end do 
      end do
    end do
    uerr = uerr / real(d%np(1), wp) / real(d%nc(2), wp) / real(d%nc(3), wp)
    uerr = sqrt_wp(uerr)

    write(output_unit, '(1F10.4, 2ES15.7)') f%time, uerr, uerrmax
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
    implicit none

    integer(4) :: iter, isub
    real(wp)   :: rtmp

    niter = nIterFlowEnd
    flow%time = tFlow 

    do iter = nrsttckpt + 1, niter
      call Call_cpu_time(CPU_TIME_ITER_START, nrsttckpt, niter, iter)
  !===============================================================================
  !      setting up 1/re, 1/re/prt, gravity, etc
  !===============================================================================
      call Calculate_RePrGr(flow, thermo, iter)
  !===============================================================================
  !     main solver
  !===============================================================================
      flow%time = flow%time + dt
      do isub = 1, nsubitr
        !if(is_thermo) call Solve_energy_eq  (flow, thermo, domain, isub)
        call Compute_burgers_rhs(flow, domain, isub)
      end do
  !
      !comment this part code for testing 
      ! below is for validation
      ! cpu time will be calculated later today 
  !===============================================================================
  !     validation
  !===============================================================================
      call Validate_burgers_error (flow, domain)

      call Call_cpu_time(CPU_TIME_ITER_END, nrsttckpt, niter, iter)

    end do

    return
  end subroutine Solve_burgers_eq_iteration

end module
