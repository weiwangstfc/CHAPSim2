
module bc_convective_outlet_mod
  use parameters_constant_mod
  use decomp_2d
  use udf_type_mod
  implicit none 

  private :: get_convective_outlet_ux
  private :: calculate_fbcx_convective_outlet
  !public  :: update_flow_from_dyn_fbcx

  private  :: enforce_domain_mass_balance_dyn_fbc
  !private :: enforce_domain_energy_balance_dyn_fbc

  !public  :: update_dyn_fbcx_from_flow
  public  :: update_fbcx_convective_outlet_flow
  public  :: update_fbcx_convective_outlet_thermo

  contains

!==========================================================================================================
  subroutine get_convective_outlet_ux(fl, dm, uxdx)
    use wtformat_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow), intent(in) :: fl
    real(WP), intent(out) :: uxdx
    logical :: flg_bc_conv

    real(WP) :: uxmax, uxmin, uxmax_work, uxmin_work, uintf
    integer :: nn, k, j
    
    if(.not. dm%is_conv_outlet) return

    uxmax = MINN
    uxmin = MAXP
    nn = dm%dpcc%xsz(1) -  1
    do k = 1, dm%dpcc%xsz(3)
      do j = 1, dm%dpcc%xsz(2)
        uintf = fl%qx(nn, j, k) !( fl%qx(nn, j, k) + dm%fbcx_qx(2, j, k) ) * HALF ! at i=nc
        if(fl%qx(nn, j, k) > uxmax) uxmax = uintf
        if(fl%qx(nn, j, k) < uxmin) uxmin = uintf
      end do
    end do

    !write(*,*) 'outlet bc', fl%qx(nn, :, 2)

    call MPI_ALLREDUCE(uxmax, uxmax_work, 1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
    call MPI_ALLREDUCE(uxmin, uxmin_work, 1, MPI_REAL_WP, MPI_MIN, MPI_COMM_WORLD, ierror)

    uxdx = HALF * (uxmax_work + uxmin_work)
    uxdx = uxdx * dm%h1r(1)
    if(nrank == 0) write(*, '(10X, A, 3ES13.5, A, 1I5.1)') 'convective outlet uxmax, min, ave = ', &
      uxmax_work, uxmin_work, HALF * (uxmax_work + uxmin_work), ' at iter(real) =', fl%iteration

    return
  end subroutine
!==========================================================================================================
  subroutine calculate_fbcx_convective_outlet(fbcx_var, uxdx, fbc_rhs0, var, dtmp, dm, isub)

    type(DECOMP_INFO), intent(in) :: dtmp 
    type(t_domain), intent(in) :: dm
    real(WP), dimension(4,           dtmp%xsz(2), dtmp%xsz(3)), intent(inout) :: fbcx_var
    real(WP), dimension(             dtmp%xsz(2), dtmp%xsz(3)), intent(inout) :: fbc_rhs0
    real(WP), dimension(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3)), intent(inout) :: var
    real(WP), intent(in) :: uxdx
    integer, intent(in) :: isub

    integer :: j, k, nn
    logical :: is_x
    real(WP) :: rhs_explicit_current, rhs_explicit_last, rhs_total
    
    if(.not. dm%is_conv_outlet) return

    ! all based on x pencil
    ! dphi/dt + ux * dphi/dx = 0
    ! data storage:
    ! ux = cell centre
    ! 
    if(dtmp%xsz(1) == dm%dpcc%xsz(1)) then
      ! qx,  -----|-----||
      !          qx     bc2/bc4 
      is_x = .true.
      nn = dm%dpcc%xsz(1) - 1
    else
      ! any vars else, like v, w, phi, T, etc 
      ! qy,  --x--|--x--||--x--|
      !       qy    qy  bc2 bc4
      is_x = .false.
      nn = dm%dccc%xsz(1)
    end if

    do k = 1, dtmp%xsz(3)
      do j = 1, dtmp%xsz(2)
      ! add explicit terms : convection rhs
        rhs_explicit_current = fbcx_var(4, j, k) - var(nn, j, k) ! at cell centre for ux, and bc point for otherse
        rhs_explicit_current = - rhs_explicit_current * uxdx
        rhs_explicit_last    = fbc_rhs0(j, k)
        rhs_total = ( dm%tGamma(isub) * rhs_explicit_current + &
                      dm%tZeta (isub) * rhs_explicit_last ) * dm%dt
        fbc_rhs0(j, k) = rhs_explicit_current
      ! calculate updated b.c. values
        fbcx_var(4, j, k) = fbcx_var(4, j, k) + rhs_total
      end do
    end do

    if(is_x) then
      ! ux, fbc = var(last point)
      do k = 1, dtmp%xsz(3)
        do j = 1, dtmp%xsz(2)
          fbcx_var(2, j, k) = fbcx_var(4, j, k)
          var(dm%dpcc%xsz(1), j, k) = fbcx_var(2, j, k)
        end do
      end do
    else
      do k = 1, dtmp%xsz(3)
        do j = 1, dtmp%xsz(2)
          fbcx_var(2, j, k) = (fbcx_var(4, j, k) + var(nn, j, k)) * HALF
        end do
      end do
    end if

  end subroutine 

!==========================================================================================================
  ! subroutine update_dyn_fbcx_from_flow(dm, ux, uy, uz, fbcx1, fbcx2, fbcx3)
  !   use print_msg_mod
  !   implicit none 
  !   type(t_domain), intent(in) :: dm
  !   real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent (in) :: ux
  !   real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent (in) :: uy
  !   real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent (in) :: uz
  !   real(WP), dimension(4,              dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent (inout) :: fbcx1
  !   real(WP), dimension(4,              dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent (inout) :: fbcx2
  !   real(WP), dimension(4,              dm%dccp%xsz(2), dm%dccp%xsz(3)), intent (inout) :: fbcx3

  !   if( .not. dm%is_conv_outlet) return

  !   ! x - pencil 
  !   if(dm%ibcx_nominal(2, 1) == IBC_CONVECTIVE) then
  !     !fbcx1(2, :, :) = ux(dm%dpcc%xsz(1), :, :)
  !     fbcx1(4, :, :) = fbcx1(2, :, :)
  !   end if
  !   if(dm%ibcx_nominal(2, 2) == IBC_CONVECTIVE) then
  !     fbcx2(4, :, :) = TWO * fbcx2(2, :, :) - uy(dm%dcpc%xsz(1), :, :)
  !   end if
  !   if(dm%ibcx_nominal(2, 3) == IBC_CONVECTIVE) then
  !     fbcx3(4, :, :) = TWO * fbcx3(2, :, :) - uz(dm%dccp%xsz(1), :, :)
  !   end if

  !   return
  ! end subroutine


!==========================================================================================================
  ! subroutine update_flow_from_dyn_fbcx(dm, ux, uy, uz, fbcx1, fbcx2, fbcx3)
  !   use udf_type_mod
  !   use parameters_constant_mod
  !   use print_msg_mod
  !   implicit none 
  !   type(t_domain), intent(in) :: dm
  !   real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent (inout) :: ux
  !   real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent (inout) :: uy
  !   real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent (inout) :: uz
  !   real(WP), dimension(4,              dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent (in)    :: fbcx1
  !   real(WP), dimension(4,              dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent (in)    :: fbcx2
  !   real(WP), dimension(4,              dm%dccp%xsz(2), dm%dccp%xsz(3)), intent (in)    :: fbcx3

  !   if( .not. dm%is_conv_outlet) return

  !   ! x - pencil 
  !   if(dm%ibcx_nominal(2, 1) == IBC_CONVECTIVE) then
  !     ux(dm%dpcc%xsz(1), :, :) = fbcx1(2, :, :)
  !   end if
  !   if(dm%ibcx_nominal(2, 2) == IBC_CONVECTIVE) then
  !     uy(dm%dcpc%xsz(1), :, :) = TWO * fbcx2(2, :, :) - fbcx2(4, :, :)
  !   end if
  !   if(dm%ibcx_nominal(2, 3) == IBC_CONVECTIVE) then
  !     uz(dm%dccp%xsz(1), :, :) =  TWO * fbcx3(2, :, :) - fbcx3(4, :, :)
  !   end if
    

  !   return
  ! end subroutine
  !==========================================================================================================
  subroutine enforce_domain_mass_balance_dyn_fbc(fl, dm)
    use wtformat_mod
    use bc_dirichlet_mod
    use find_max_min_ave_mod
    implicit none
    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(inout) :: dm
    type(DECOMP_INFO) :: dtmp
    real(WP), dimension(4, dm%dpcc%xsz(2), dm%dpcc%xsz(3)) :: fbcx
    real(WP), dimension(dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3)) :: fbcy
    real(WP), dimension(dm%dccp%zsz(1), dm%dccp%zsz(2), 4) :: fbcz
    real(WP) :: mass_rate_net
    real(WP) :: mass_rate_scaling
    real(WP) :: mass_rate_iin(3), mass_rate_iin_work(3)
    real(WP) :: mass_rate_out(3), mass_rate_out_work(3)
    real(WP) :: mass_rate_core, mass_rate_core_work
    real(WP) :: mass_rate_iin_net, mass_rate_out_net
    real(WP) :: dummy(7), dummy_work(7)
    real(WP) :: dx, dy, dz
    integer :: nn, i, j, k, jj
    
    if (.not. dm%is_conv_outlet) return

    mass_rate_iin = ZERO
    mass_rate_out = ZERO
    mass_rate_core= ZERO

    mass_rate_net = ZERO
    mass_rate_scaling = ONE

    if(dm%is_thermo) then
      fbcx = dm%fbcx_gx
      fbcy = dm%fbcy_gy
      fbcz = dm%fbcz_gz
      call Get_volumetric_average_3d_for_var_xcx(dm, dm%dccc, fl%drhodt, mass_rate_core, LF3D_VOL_SUM, "mass_rate_core")
    else
      fbcx = dm%fbcx_qx
      fbcy = dm%fbcy_qy
      fbcz = dm%fbcz_qz
      mass_rate_core = ZERO
    end if

!----------------------------------------------------------------------------------------------------------
! x - inlet/outlet
!----------------------------------------------------------------------------------------------------------
    nn = 1
    if( dm%ibcx_nominal(2, nn) == IBC_CONVECTIVE) then
      dtmp = dm%dpcc
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1 !local2global_yid(j, dtmp)
        dy = dm%h(2)/dm%yMappingcc(jj, 1) !dm%yp(jj+1) - dm%yp(jj)
        do k = 1, dtmp%xsz(3)
          dz = dm%h(3) * dm%rc(jj)
          mass_rate_iin(nn) = mass_rate_iin(nn) + fbcx(1, j, k) * dy * dz
          mass_rate_out(nn) = mass_rate_out(nn) + fbcx(2, j, k) * dy * dz
        end do
      end do 
    end if
!----------------------------------------------------------------------------------------------------------
! y - inlet/outlet - still x-pencil
!----------------------------------------------------------------------------------------------------------
    nn = 2
    if( dm%ibcy_nominal(2, nn) == IBC_CONVECTIVE) then
      dtmp = dm%dcpc
      do i = 1, dtmp%xsz(1)
        dx = dm%h(1)
        do k = 1, dtmp%xsz(3)
          dz = dm%h(3)
          if(dtmp%xst(nn) ==         1) mass_rate_iin(nn) = mass_rate_iin(nn) + fbcy(i, 1, k) * dx * dz * dm%rc(1)
          if(dtmp%xen(nn) == dm%np(nn)) mass_rate_out(nn) = mass_rate_out(nn) + fbcy(i, 2, k) * dx * dz * dm%rc(dm%np(2))
        end do
      end do 
    end if
!----------------------------------------------------------------------------------------------------------
! z - inlet/outlet - still x-pencil
!----------------------------------------------------------------------------------------------------------
    nn = 3
    if( dm%ibcz_nominal(2, nn) == IBC_CONVECTIVE) then
      dtmp = dm%dccp
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1 !local2global_yid(j, dtmp)
        dy = dm%h(2)/dm%yMappingcc(jj, 1)
        do i = 1, dtmp%xsz(1)
          dx = dm%h(1)
          if(dtmp%xst(nn) ==        1)  mass_rate_iin(nn) = mass_rate_iin(nn) + fbcz(i, j, 1) * dx * dy
          if(dtmp%xen(nn) == dm%np(nn)) mass_rate_out(nn) = mass_rate_out(nn) + fbcz(i, j, 2) * dx * dy
        end do
      end do 
    end if
!----------------------------------------------------------------------------------------------------------
! add from all ranks
!----------------------------------------------------------------------------------------------------------
    dummy(1:3) = mass_rate_iin(1:3)  ! unit: kg/s
    dummy(4:6) = mass_rate_out(1:3)  ! unit: kg/s 
    dummy(7)   = mass_rate_core ! unit: kg/s

    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce( dummy,  dummy_work, 7, MPI_REAL_WP, MPI_SUM, MPI_COMM_WORLD, ierror)

    mass_rate_iin_work(1:3)  = dummy(1:3)
    mass_rate_out_work(1:3)  = dummy(4:6)
    mass_rate_core_work = dummy(7) 
!----------------------------------------------------------------------------------------------------------
! scaling factor for a mass conservation
!----------------------------------------------------------------------------------------------------------
    mass_rate_iin_net = mass_rate_iin_work(1) + mass_rate_iin_work(2) + mass_rate_iin_work(3)
    mass_rate_out_net = mass_rate_out_work(1) + mass_rate_out_work(2) + mass_rate_out_work(3)
    do nn = 1, 3
      mass_rate_net = mass_rate_core_work + mass_rate_iin_net - mass_rate_out_net! check 1st term plus or minus?                
      mass_rate_scaling = (mass_rate_core_work + mass_rate_iin_net) / mass_rate_out_net
    end do

!#ifdef DEBUG_STEPS 
    if(nrank == 0) then 
      !write (*, *) ">------convective outlet-------------->"
      !write (*, *) "mass_in(3) = ", mass_rate_iin_work(:)
      !write (*, *) "massout(3) = ", mass_rate_out_work(:)
      write (*, '(10X, A, 4ES13.5)') " mass_rate_iin_net, out_net, core = ", mass_rate_iin_net, mass_rate_out_net, mass_rate_core_work
      write (*, '(10X, A, 1ES13.5, 1F16.13)') " mass rate net change and scaling = ", mass_rate_net, mass_rate_scaling
      !write (*, *) "<------convective outlet--------------<"
    end if
!#endif
!----------------------------------------------------------------------------------------------------------
! scale the dynamic bc
!----------------------------------------------------------------------------------------------------------
    if( dm%ibcx_nominal(2, 1) == IBC_CONVECTIVE) then
      fbcx(2, :, :) = fbcx(2, :, :) * mass_rate_scaling
      fbcx(4, :, :) = fbcx(2, :, :)
    end if
    if( dm%ibcy_nominal(2, 2) == IBC_CONVECTIVE) then
      fbcy(:, 2, :) = fbcy(:, 2, :) * mass_rate_scaling
      fbcy(:, 4, :) = fbcy(:, 2, :)
    end if
    if( dm%ibcz_nominal(2, 3) == IBC_CONVECTIVE) then
      fbcz(:, :, 2) = fbcz(:, :, 2) * mass_rate_scaling
      fbcz(:, :, 2) = fbcz(:, :, 4)
    end if 
!----------------------------------------------------------------------------------------------------------
! back to real fbc
!----------------------------------------------------------------------------------------------------------
    if(dm%is_thermo) then
      if( dm%ibcx_nominal(2, 1) == IBC_CONVECTIVE) dm%fbcx_gx(2, :, :) = fbcx(2, :, :) 
      if( dm%ibcy_nominal(2, 2) == IBC_CONVECTIVE) dm%fbcy_gy(2, :, :) = fbcy(2, :, :) 
      if( dm%ibcz_nominal(2, 3) == IBC_CONVECTIVE) dm%fbcz_gz(2, :, :) = fbcz(2, :, :) 
      if( dm%ibcx_nominal(2, 1) == IBC_CONVECTIVE) dm%fbcx_qx(2, :, :) = dm%fbcx_gx(2, :, :) / dm%fbcx_ftp(2, :, :)%d
      if( dm%ibcy_nominal(2, 2) == IBC_CONVECTIVE) dm%fbcy_qy(2, :, :) = dm%fbcy_gy(2, :, :) / dm%fbcy_ftp(2, :, :)%d
      if( dm%ibcz_nominal(2, 3) == IBC_CONVECTIVE) dm%fbcz_qz(2, :, :) = dm%fbcz_gz(2, :, :) / dm%fbcz_ftp(2, :, :)%d
      !call update_flow_from_dyn_fbcx(dm, fl%gx, fl%gy, fl%gz, dm%fbcx_gx, dm%fbcx_gy, dm%fbcx_gz)
    else
      if( dm%ibcx_nominal(2, 1) == IBC_CONVECTIVE) dm%fbcx_qx(2, :, :) = fbcx(2, :, :)
      if( dm%ibcy_nominal(2, 2) == IBC_CONVECTIVE) dm%fbcy_qy(2, :, :) = fbcy(2, :, :)
      if( dm%ibcz_nominal(2, 3) == IBC_CONVECTIVE) dm%fbcz_qz(2, :, :) = fbcz(2, :, :)
    end if
    !call update_flow_from_dyn_fbcx(dm, fl%qx, fl%qy, fl%qz, dm%fbcx_qx, dm%fbcx_qy, dm%fbcx_qz)
    
    return
  end subroutine enforce_domain_mass_balance_dyn_fbc

!==========================================================================================================
  subroutine update_fbcx_convective_outlet_flow(fl, dm, isub)
    use bc_dirichlet_mod
    implicit none
    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(inout) :: dm
    integer,        intent(in)    :: isub
    
    real(WP) :: uxdx
    integer :: i

    if(.not. dm%is_conv_outlet) return
    if(nrank == 0) call Print_debug_mid_msg("Calculate convective outlet for flow ...")
    ! work on fbcx, not fl directly
    call get_convective_outlet_ux(fl, dm, uxdx)
    if ( .not. dm%is_thermo) then
      
      if(dm%ibcx_nominal(2, 1) == IBC_CONVECTIVE) then
        call calculate_fbcx_convective_outlet(dm%fbcx_qx(:, :, :), uxdx, fl%fbcx_qx_rhs0(:, :), fl%qx, dm%dpcc, dm, isub)
      end if
      if(dm%ibcx_nominal(2, 2) == IBC_CONVECTIVE) then
        call calculate_fbcx_convective_outlet(dm%fbcx_qy(:, :, :), uxdx, fl%fbcx_qy_rhs0(:, :), fl%qy, dm%dcpc, dm, isub)
      end if 
      if(dm%ibcx_nominal(2, 3) == IBC_CONVECTIVE) then
        call calculate_fbcx_convective_outlet(dm%fbcx_qz(:, :, :), uxdx, fl%fbcx_qz_rhs0(:, :), fl%qz, dm%dccp, dm, isub)
      end if
    else
      ! check , whether it is better to use qx = gx / density ?
      if(dm%ibcx_nominal(2, 1) == IBC_CONVECTIVE) then
        call calculate_fbcx_convective_outlet(dm%fbcx_qx(:, :, :), uxdx, fl%fbcx_qx_rhs0(:, :), fl%qx, dm%dpcc, dm, isub)
        call calculate_fbcx_convective_outlet(dm%fbcx_gx(:, :, :), uxdx, fl%fbcx_gx_rhs0(:, :), fl%gx, dm%dpcc, dm, isub)
      end if
      if(dm%ibcx_nominal(2, 2) == IBC_CONVECTIVE) then
        call calculate_fbcx_convective_outlet(dm%fbcx_qy(:, :, :), uxdx, fl%fbcx_qy_rhs0(:, :), fl%qy, dm%dcpc, dm, isub)
        call calculate_fbcx_convective_outlet(dm%fbcx_gy(:, :, :), uxdx, fl%fbcx_gy_rhs0(:, :), fl%gy, dm%dcpc, dm, isub)
      end if
      if(dm%ibcx_nominal(2, 3) == IBC_CONVECTIVE) then
        call calculate_fbcx_convective_outlet(dm%fbcx_qz(:, :, :), uxdx, fl%fbcx_qz_rhs0(:, :), fl%qz, dm%dccp, dm, isub)
        call calculate_fbcx_convective_outlet(dm%fbcx_gz(:, :, :), uxdx, fl%fbcx_gz_rhs0(:, :), fl%gz, dm%dccp, dm, isub)
      end if

    end if

    call enforce_domain_mass_balance_dyn_fbc(fl, dm)

    if ( .not. dm%is_thermo) then
      call enforce_velo_from_fbc(dm, fl%qx, fl%qy, fl%qz, dm%fbcx_qx, dm%fbcy_qy, dm%fbcz_qz)
    else
      ! check , whether it is better to use qx = gx / density ?
      call enforce_velo_from_fbc(dm, fl%qx, fl%qy, fl%qz, dm%fbcx_qx, dm%fbcy_qy, dm%fbcz_qz)
      call enforce_velo_from_fbc(dm, fl%gx, fl%gy, fl%gz, dm%fbcx_gx, dm%fbcy_gy, dm%fbcz_gz)
    end if

    return
  end subroutine

!==========================================================================================================
  subroutine update_fbcx_convective_outlet_thermo(fl, tm, dm, isub)
    use thermo_info_mod
    implicit none
    type(t_flow),   intent(inout) :: fl
    type(t_thermo), intent(inout) :: tm
    type(t_domain), intent(inout) :: dm
    integer,        intent(in)    :: isub
    
    real(WP) :: uxdx
    integer :: j, k
    real(WP), dimension(4, dm%dccc%xsz(2), dm%dccc%xsz(3)) :: a4cc_xpencil

    if ( .not. dm%is_thermo) return
    if ( .not. dm%is_conv_outlet) return

    if(nrank == 0) call Print_debug_mid_msg("Calculate convective outlet for thermo ...")
    call get_convective_outlet_ux(fl, dm, uxdx)

    if(dm%ibcx_nominal(2, 5) == IBC_CONVECTIVE) then
      a4cc_xpencil = dm%fbcx_ftp(:, :, :)%rhoh
      call calculate_fbcx_convective_outlet(a4cc_xpencil, uxdx, tm%fbcx_rhoh_rhs0(:, :), tm%rhoh, dm%dccc, dm, isub)
      dm%fbcx_ftp(:, :, :)%rhoh = a4cc_xpencil
      do j = 1, size(dm%fbcx_ftp, 2)
        do k = 1, size(dm%fbcx_ftp, 3)
          call ftp_refresh_thermal_properties_from_DH(dm%fbcx_ftp(2, j, k))
          call ftp_refresh_thermal_properties_from_DH(dm%fbcx_ftp(4, j, k))
        end do
      end do
    end if

    !call enforce_domain_energy_balance_dyn_fbc(fl, dm) ! to check necessary? 

    return
  end subroutine

end module
