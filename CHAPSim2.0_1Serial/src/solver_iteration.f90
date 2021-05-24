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
!      setting up 1/re, 1/re/prt, gravity, etc
!===============================================================================
       call Calculate_parameters_in_eqs(flow, thermo, iter)
!===============================================================================
!      setting up flow solver
!===============================================================================
      if ( (iter >= nIterFlow0) .and. (iter <=nIterFlow1)) then
        is_flow = .true.
        flow%t = flow%t + dt
        call Check_cfl_diffusion (domain%h2r(:), flow%rre)
        call Check_cfl_convection(flow%qx, flow%qy, flow%qz, domain)
      end if
!===============================================================================
!     setting up thermo solver
!===============================================================================
      if ( (iter >= nIterThermo0) .and. (iter <=nIterThermo1)) then
        is_thermo = .true.
        thermo%t = thermo%t + dt
      end if
!===============================================================================
!     main solver
!===============================================================================
      do isub = 1, nsubitr
        if(is_thermo) call Solve_energy_eq (flow, thermo, domain, isub)
        if(is_flow)   call Solve_momentum_eq (flow, domain, isub)
      end do

    end do



  end subroutine Solve_eqs_iteration


end module solver_iteration