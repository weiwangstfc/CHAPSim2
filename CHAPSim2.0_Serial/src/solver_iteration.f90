module solver_iteration_mod


contains

  subroutine Solve_eqs_iteration
    use input_general_mod
    use parameters_constant_mod
    implicit none
    integer(4) :: iter
    logical    :: is_flow   = .false.
    logical    :: is_thermo = .false.
    
    if(ithermo == 1) then
      niter = MAX(nIterFlow1, nIterThermo1)
      thermo%t = tThermo
    else
      niter = nIterFlow1
      flow%t = tFlow
    end if

    do iter = iterchkpt + 1, niter
!===============================================================================
!      setting up flow solver
!===============================================================================
      if ( (iter >= nIterFlow0) .and. (iter <=nIterFlow1)) then
        is_flow = .true.
        flow%t = flow%t + dt
!-------------------------------------------------------------------------------
!       ->set up Reynolds number, which may ramp up/down based on user's setting
!_______________________________________________________________________________
        if(iter < nIterIniRen) then
          flow%rre = ONE / renIni
        else
          flow%rre = ONE / ren
        end if
!-------------------------------------------------------------------------------
!       ->check CFL number from diffusion and convection, for stability
!_______________________________________________________________________________
        call Check_cfl_diffusion (domain%h2r(:), flow%rre)
        call Check_cfl_convection(flow%qx, flow%qy, flow%qz, domain)
      end if
!===============================================================================
!     setting up thermo solver
!===============================================================================
      if ( (iter >= nIterThermo0) .and. (iter <=nIterThermo1)) then
        is_thermo = .true.
        thermo%t = thermo%t + dt
!-------------------------------------------------------------------------------
!       ->set up 1/(Pr*Re)
!_______________________________________________________________________________
        thermo%rPrRen = flow%rre * tpRef0%k / tpRef0%m / tpRef0%cp
      end if
!===============================================================================
!     main solver
!===============================================================================
      do isub = 1, nsubitr
        if(is_thermo) call Solve_energy_eq(isub)
        if(is_flow)   call Solve_momentum_eq(isub)
      end do

    end do



  end subroutine Solve_eqs_iteration


end module solver_iteration