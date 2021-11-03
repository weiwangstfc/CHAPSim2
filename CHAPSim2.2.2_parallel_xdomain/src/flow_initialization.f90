!-------------------------------------------------------------------------------
!                      CHAPSim version 2.0.0
!                      --------------------------
! This file is part of CHAPSim, a general-purpose CFD tool.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 3 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 0type t_domain=============================================
!> \file flow_initialisation.f90
!>
!> \brief Define and initialise flow and thermal variables.
!>
!===============================================================================
module flow_thermo_initialiasation
  use var_dft_mod
  implicit none

  
  private :: Generate_poiseuille_flow_profile
  private :: Initialize_poiseuille_flow
  private :: Initialize_vortexgreen_2dflow
  private :: Initialize_vortexgreen_3dflow
  private :: Initialize_thermal_variables
  private :: Generate_random_field

  private  :: Allocate_thermoflow_variables
  private  :: Initialize_flow_variables

  public  :: Calculate_xz_mean
  public  :: Check_maximum_velocity
  
  public  :: Calculate_RePrGr
  public  :: Validate_TGV2D_error
  public  :: Initialize_flow_thermal_fields

contains
!===============================================================================
!===============================================================================
!> \brief Allocate flow and thermal variables.     
!>
!> This subroutine is called once at beginning of solver.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
  subroutine Allocate_thermoflow_variables (d, f, t)
    use domain_decomposition_mod
    use input_general_mod,       only : ithermo
    use parameters_constant_mod, only : ZERO, ONE
    use mpi_mod
    implicit none

    type(t_domain), intent(in)    :: d
    type(t_flow),   intent(inout) :: f
    type(t_thermo), intent(inout) :: t

    if(nrank == 0) call Print_debug_start_msg("Allocating flow and thermal variables ...")

!_______________________________________________________________________________
! x pencil. 
! varaible index is LOCAL. means 1:xsize(1)
!_______________________________________________________________________________
    call alloc_x(f%qx,      d%dpcc) ; f%qx = ZERO
    call alloc_x(f%qy,      d%dcpc) ; f%qy = ZERO
    call alloc_x(f%qz,      d%dccp) ; f%qz = ZERO
    
    call alloc_x(f%gx,      d%dpcc) ; f%gx = ZERO
    call alloc_x(f%gy,      d%dcpc) ; f%gy = ZERO
    call alloc_x(f%gz,      d%dccp) ; f%gz = ZERO

    call alloc_x(f%pres,    d%dccc) ; f%pres = ZERO
    call alloc_x(f%pcor,    d%dccc) ; f%pcor = ZERO

    call alloc_x(f%dDens,   d%dccc) ; f%dDens = ONE
    call alloc_x(f%mVisc,   d%dccc) ; f%mVisc = ONE

    call alloc_x(f%dDensm1, d%dccc) ; f%dDensm1 = ONE
    call alloc_x(f%dDensm2, d%dccc) ; f%dDensm2 = ONE

    call alloc_x(f%mx_rhs,  d%dpcc) ; f%mx_rhs = ZERO
    call alloc_x(f%my_rhs,  d%dcpc) ; f%my_rhs = ZERO
    call alloc_x(f%mz_rhs,  d%dccp) ; f%mz_rhs = ZERO

    call alloc_x(f%mx_rhs0, d%dpcc) ; f%mx_rhs0 = ZERO
    call alloc_x(f%my_rhs0, d%dcpc) ; f%my_rhs0 = ZERO
    call alloc_x(f%mz_rhs0, d%dccp) ; f%mz_rhs0 = ZERO

    if(d%ithermo == 1) then
      call alloc_x(t%dh,    d%dccc) ; f%dh    = ZERO
      call alloc_x(t%hEnth, d%dccc) ; f%hEnth = ZERO
      call alloc_x(t%kCond, d%dccc) ; f%kCond = ONE
      call alloc_x(t%tTemp, d%dccc) ; f%tTemp = ONE
    end if
!_______________________________________________________________________________
! x pencil. 
!_______________________________________________________________________________

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine Allocate_thermoflow_variables
!===============================================================================
!===============================================================================
!> \brief Initialise thermal variables if ithermo = 1.     
!>
!> This subroutine is called once.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  f             flow type
!> \param[inout]  t             thermo type
!_______________________________________________________________________________
  subroutine Initialize_thermal_variables (f, t, d)
    use parameters_constant_mod
    use input_thermo_mod
    implicit none
    type(t_flow),   intent(inout) :: f
    type(t_thermo), intent(inout) :: t
    type(t_domain), intent(in)    :: d
    
    type(thermoProperty_t) :: tp_ini
    integer :: i
!-------------------------------------------------------------------------------
!   initialize thermo fields
!-------------------------------------------------------------------------------
    tp_ini = t%tpIni
    tp_ini%t = t%tiRef / t%t0Ref

    call tp_ini%Refresh_thermal_properties_from_T_undim

    f%dDens(:, :, :) = tp_ini%d
    f%mVisc(:, :, :) = tp_ini%m

    t%dh   (:, :, :) = tp_ini%dh
    t%hEnth(:, :, :) = tp_ini%h
    t%kCond(:, :, :) = tp_ini%k
    t%tTemp(:, :, :) = tp_ini%t

    return
  end subroutine Initialize_thermal_variables
!===============================================================================
!===============================================================================
!> \brief The main code for initializing flow variables
!>
!> not changing storage position, exclude b.c. values, for example, developing
!> flow.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  none          NA
!_______________________________________________________________________________
  subroutine Calculate_xz_mean(u, uxz_work, str, d)
    use mpi_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in) :: d
    real(WP),       intent(in) :: u(:, :, :)
    real(WP),    intent(inout) :: uxz_work(:)
    character(2),   intent(in) :: str

    real(wp) :: uxz( size(uxz_work) )
    integer(4) :: jj, kk, ii, ny, i, j, k
    integer(4) :: ist, ien, jst, jen, kst, ken
    integer(4) :: xst(3), xen(3), xsz(3)
!-------------------------------------------------------------------------------
!   Default X-pencil
!-------------------------------------------------------------------------------
    if(str == 'ux') then
      xst(:) = d%ux_xst(:)
      xen(:) = d%ux_xen(:)
      xsz(:) = d%ux_xsz(:)
      ist = 1
      ien = xsz(1)
      jst = 1
      jen = xsz(2)
      kst = 1
      ken = xsz(3)
      if(xst(1)==1 .and. d%bc(1, 1) == IBC_UDIRICHLET) then
        ist = ist + 1
      end if
      if(xen(1)==d%np(1) .and. d%bc(2, 1) == IBC_UDIRICHLET)then
        ien = ien - 1
      end if
      ny = d%nc(2)
    else if (str == 'uy') then
      xst(:) = d%uy_xst(:)
      xen(:) = d%uy_xen(:)
      xsz(:) = d%uy_xsz(:)
      ist = 1
      ien = xsz(1)
      jst = 1
      jen = xsz(2)
      kst = 1
      ken = xsz(3)
      if(xst(2)==1 .and. d%bc(1, 2) == IBC_UDIRICHLET) then
        jst = jst + 1
      end if
      if(xen(2)==d%np(2) .and. d%bc(2, 2) == IBC_UDIRICHLET)then
        jen = jen - 1
      end if
      ny = d%np(2)
    else if (str == 'uz') then
      xst(:) = d%uz_xst(:)
      xen(:) = d%uz_xen(:)
      xsz(:) = d%uz_xsz(:)
      ist = 1
      ien = xsz(1)
      jst = 1
      jen = xsz(2)
      kst = 1
      ken = xsz(3)
      if(xst(3)==1 .and. d%bc(1, 3) == IBC_UDIRICHLET) then
        kst = kst + 1
      end if
      if(xen(3)==d%np(3) .and. d%bc(2, 3) == IBC_UDIRICHLET)then
        ken = ken - 1
      end if
      ny = d%nc(2)
    else if (str == 'ps') then
      xst(:) = d%ps_xst(:)
      xen(:) = d%ps_xen(:)
      xsz(:) = d%ps_xsz(:)
      ist = 1
      ien = xsz(1)
      jst = 1
      jen = xsz(2)
      kst = 1
      ken = xsz(3)
      ny = d%nc(2)
    else
      call Print_error_msg("Error, input type is wrong.")
    end if
  
    uxz(:) = ZERO
    uxz_work(:) = ZERO
    kk = 0
    ii = 0
    do j = jst, jen
      jj = j - 1 + xst(2)
      do k = kst, ken
        kk = kk + 1
        do i = ist, ien
          ii = ii + 1
          uxz(jj) = uxz(jj) + u(i, j, k)
        end do
      end do
    end do
    uxz(:) = uxz(:) / real(kk * ii, wp)

    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(uxz, uxz_work, ny, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    uxz_work(:) = uxz_work(:) / real(ncol, wp)

    return
  end subroutine
!===============================================================================
!===============================================================================
  subroutine Calculate_xzmean_perturbation(u, uxz, d, str, uprofile)
    use mpi_mod
    use udf_type_mod
    implicit none
    type(t_domain),     intent(in) :: d
    real(WP),        intent(inout) :: u(:, :, :)
    real(WP),           intent(in) :: uxz(:)
    character(2),       intent(in) :: str
    real(WP), optional, intent(in) :: uprofile(:)

    integer(4) :: jj, kk, ii, ny, i, j, k
    integer(4) :: ist, ien, jst, jen, kst, ken
    integer(4) :: xst(3), xen(3), xsz(3)
!-------------------------------------------------------------------------------
!   Default X-pencil
!   excludes b.c. index for Dirichlet B.C.
!-------------------------------------------------------------------------------
    if(str == 'ux') then
      xst(:) = d%ux_xst(:)
      xen(:) = d%ux_xen(:)
      xsz(:) = d%ux_xsz(:)
      ist = 1
      ien = xsz(1)
      jst = 1
      jen = xsz(2)
      kst = 1
      ken = xsz(3)
      if(xst(1)==1 .and. d%bc(1, 1) == IBC_UDIRICHLET) then
        ist = ist + 1
      end if
      if(xen(1)==d%np(1) .and. d%bc(2, 1) == IBC_UDIRICHLET)then
        ien = ien - 1
      end if
      ny = d%nc(2)
    else if (str == 'uy') then
      xst(:) = d%uy_xst(:)
      xen(:) = d%uy_xen(:)
      xsz(:) = d%uy_xsz(:)
      ist = 1
      ien = xsz(1)
      jst = 1
      jen = xsz(2)
      kst = 1
      ken = xsz(3)
      if(xst(2)==1 .and. d%bc(1, 2) == IBC_UDIRICHLET) then
        jst = jst + 1
      end if
      if(xen(2)==d%np(2) .and. d%bc(2, 2) == IBC_UDIRICHLET)then
        jen = jen - 1
      end if
      ny = d%np(2)
    else if (str == 'uz') then
      xst(:) = d%uz_xst(:)
      xen(:) = d%uz_xen(:)
      xsz(:) = d%uz_xsz(:)
      ist = 1
      ien = xsz(1)
      jst = 1
      jen = xsz(2)
      kst = 1
      ken = xsz(3)
      if(xst(3)==1 .and. d%bc(1, 3) == IBC_UDIRICHLET) then
        kst = kst + 1
      end if
      if(xen(3)==d%np(3) .and. d%bc(2, 3) == IBC_UDIRICHLET)then
        ken = ken - 1
      end if
      ny = d%nc(2)
    else if (str == 'ps') then
      xst(:) = d%ps_xst(:)
      xen(:) = d%ps_xen(:)
      xsz(:) = d%ps_xsz(:)
      ist = 1
      ien = xsz(1)
      jst = 1
      jen = xsz(2)
      kst = 1
      ken = xsz(3)
      ny = d%nc(2)
    else
      call Print_error_msg("Error, input type is wrong.")
    end if
!-------------------------------------------------------------------------------
!   X-pencil
!   calculate perturbations
!-------------------------------------------------------------------------------
    do j = jst, jen
      jj = j - 1 + xst(2)
      do k = kst, ken
        do i = ist, ien
          if( present(uprofile) ) then
            u(:, j, :) = u(:, j, :) - uxz(jj) + ufyc(jj)
          else
            u(:, j, :) = u(:, j, :) - uxz(jj)
          end if
        end do
      end do
    end do

    return
  end subroutine
!===============================================================================
!===============================================================================
!> \brief Generate a flow profile for Poiseuille flow in channel or pipe.     
!>
!> This subroutine is called locally once.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!> \param[out]    ufyc          u(yc), velocity profile along wall-normal direction
!_______________________________________________________________________________
  subroutine Generate_poiseuille_flow_profile(ufyc, d)
    use parameters_constant_mod, only : ZERO, ONE, ONEPFIVE, TWO, MAXP, TRUNCERR
    use input_general_mod
    use udf_type_mod
    use math_mod
    implicit none

    type(t_domain), intent(in)  :: d
    real(WP),       intent(out) :: ufyc(:)
    
    real(WP)   :: a, b, c, yy, ymax, ymin
    integer(4) :: j
    

    ufyc (:) = ZERO

    ymax = d%yp( d%np_geo(2) )
    ymin = d%yp( 1 )
    if (d%icase == ICASE_CHANNEL) then
      a = (ymax - ymin) / TWO
      b = ZERO
      c = ONEPFIVE
    else if (d%icase == ICASE_PIPE) then
      a = (ymax - ymin)
      b = ZERO
      c = TWO
    else if (d%icase == ICASE_ANNUAL) then
      a = (ymax - ymin) / TWO
      b = (ymax + ymin) / TWO
      c = TWO
    else 
      a = MAXP
      b = ZERO
      c = ONE
    end if

    do j = 1, d%nc(2)
      yy = d%yc(j)
      ufyc(j) = ( ONE - ( (yy - b)**2 ) / a / a ) * c
    end do

    return
  end subroutine Generate_poiseuille_flow_profile
!===============================================================================
!===============================================================================
  subroutine Generate_random_field(ux, uy, uz, lnoise, d)
    use random_number_generation_mod
    use parameters_constant_mod, only : ZERO, ONE
    use input_general_mod
    use mpi_mod
    implicit none
    type(t_domain), intent(in) :: d
    real(WP),    intent(inout) :: ux(:, :, :)
    real(WP),    intent(inout) :: uy(:, :, :)
    real(WP),    intent(inout) :: uz(:, :, :)
    real(WP),    intent(in)    :: lnoise
    integer(4) :: seed, i, j, k, ii, jj, kk
!-------------------------------------------------------------------------------
!   Ensure it is in x-pencil
!-------------------------------------------------------------------------------
    if(d%ux_xsz(1) /= d%np(1)) call Print_error_msg("Error, not X-pencil")
!-------------------------------------------------------------------------------
!   Initialisation
!-------------------------------------------------------------------------------
    seed = 0
    ux(:, :, :) = ZERO
    uy(:, :, :) = ZERO
    uz(:, :, :) = ZERO
!-------------------------------------------------------------------------------
!   x-pencil : global index -->local index, ux
!-------------------------------------------------------------------------------
    do ii = 1, d%np(1)
      i = ii
      do jj = 1, d%nc(2)
        if(jj >= d%ux_xst(2) .and. jj <= d%ux_xen(2) ) then
          j = jj - d%ux_xst(2) + 1
          do kk = 1, d%nc(3)
            if(kk >= d%ux_xst(3) .and. kk <= d%ux_xen(3) ) then
              k = kk - d%ux_xst(3) + 1
              seed = ii + jj + kk
              call Initialize_random_number ( seed )
              call Generate_r_random( -ONE, ONE, rd)
              ux(i, j, k) = lnoise * rd
            end if
          end do
        end if
      end do
    end do

!-------------------------------------------------------------------------------
!   x-pencil : global index -->local index, uy
!-------------------------------------------------------------------------------
    do ii = 1, d%nc(1)
      i = ii
      do jj = 1, d%np(2)
        if(jj >= d%uy_xst(2) .and. jj <= d%uy_xen(2) ) then
          j = jj - d%uy_xst(2) + 1
          do kk = 1, d%nc(3)
            if(kk >= d%uy_xst(3) .and. kk <= d%uy_xen(3) ) then
              k = kk - d%uy_xst(3) + 1
              seed = seed + ii + jj + kk
              call Initialize_random_number ( seed )
              call Generate_r_random( -ONE, ONE, rd)
              uy(i, j, k) = lnoise * rd
            end if
          end do
        end if
      end do
    end do

!-------------------------------------------------------------------------------
!   x-pencil : global index -->local index, uz
!-------------------------------------------------------------------------------
    do ii = 1, d%nc(1)
      i = ii
      do jj = 1, d%nc(2)
        if(jj >= d%uz_xst(2) .and. jj <= d%uz_xen(2) ) then
          j = jj - d%uz_xst(2) + 1
          do kk = 1, d%np(3)
            if(kk >= d%uz_xst(3) .and. kk <= d%uz_xen(3) ) then
              k = kk - d%uz_xst(3) + 1
              seed = seed + ii + jj + kk
              call Initialize_random_number ( seed )
              call Generate_r_random( -ONE, ONE, rd)
              uz(i, j, k) = lnoise * rd
            end if
          end do
        end if
      end do
    end do

    return
  end subroutine
!===============================================================================
!===============================================================================
!> \brief Initialize Poiseuille flow in channel or pipe.     
!>
!> This subroutine is called locally once.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!> \param[out]    f             flow
!_______________________________________________________________________________
  subroutine Initialize_poiseuille_flow(ux, uy, uz, p, lnoise, d)
    use input_general_mod
    use udf_type_mod
    use boundary_conditions_mod
    implicit none
    type(t_domain), intent(in   ) :: d
    real(WP),       intent(inout) :: ux(:, :, :) , &
                                     uy(:, :, :) , &
                                     uz(:, :, :) , &
                                     p (:, :, :)    
    real(WP), intent(in) :: lnoise        
    real(WP) :: ufyc(d%nc(2))
    integer :: seed
    real(WP) :: rd
    integer :: pf_unit
    real(WP) :: uxxza(d%nc(2))
    real(WP) :: uyxza(d%np(2))
    real(WP) :: uzxza(d%nc(2))
    real(WP) :: ux_ypencil(d%ux_ysz(1), d%ux_ysz(2), d%ux_ysz(3))
!-------------------------------------------------------------------------------
!   Ensure it is in x-pencil
!-------------------------------------------------------------------------------
    if(d%ux_xsz(1) /= d%np(1)) call Print_error_msg("Error, not X-pencil")
!-------------------------------------------------------------------------------
!   x-pencil : to get Poiseuille profile
!-------------------------------------------------------------------------------
    ufyc(:) = ZERO
    call Generate_poiseuille_flow_profile ( ufyc, d )
!-------------------------------------------------------------------------------
!   x-pencil : to get random fields [-1,1] for ux, uy, uz
!-------------------------------------------------------------------------------
    p (:, :, :) = ZERO
    ux(:, :, :) = ZERO
    uy(:, :, :) = ZERO
    uz(:, :, :) = ZERO
    seed = 0
    call Generate_random_field(ux, uy, uz, lnoise, d)
!-------------------------------------------------------------------------------
!   x-pencil : build up boundary
!-------------------------------------------------------------------------------
    call Apply_BC_velocity(ux, uy, uz, d)
!-------------------------------------------------------------------------------
!   x-pencil : Get the u, v, w, averged in x and z directions
!-------------------------------------------------------------------------------
    call Calculate_xz_mean(ux, uxxza, 'ux', d)
    call Calculate_xz_mean(uy, uyxza, 'uy', d)
    call Calculate_xz_mean(uz, uzxza, 'uz', d)
!-------------------------------------------------------------------------------
!   x-pencil : Ensure u, v, w, averaged in x and z direction is zero.
!-------------------------------------------------------------------------------
    call Calculate_xzmean_perturbation(ux, uxxza, 'ux', d, ufyc)
    call Calculate_xzmean_perturbation(uy, uyxza, 'uy', d)
    call Calculate_xzmean_perturbation(uz, uzxza, 'uz', d)

!-------------------------------------------------------------------------------
!   x-pencil : Ensure u, v, w, averaged in x and z direction is zero.
!-------------------------------------------------------------------------------
    call Get_volumetric_average_3d(ux, 'ux', d, ubulk)
    ux(:, :, :) = ux(:, :, :) / ubulk
    call Apply_BC_velocity(ux, uy, uz, d)
    call Get_volumetric_average_3d(ux, 'ux', d, ubulk)
!-------------------------------------------------------------------------------
!   X-pencil ==> Y-pencil
!-------------------------------------------------------------------------------
      call transpose_x_to_y(ux, ux_ypencil, d%dpcc)
!-------------------------------------------------------------------------------
!   Y-pencil : write out velocity profile
!-------------------------------------------------------------------------------
    if(nrank == 0) then
      open ( newunit = pf_unit,     &
              file    = 'output_check_poiseuille_profile.dat', &
              status  = 'replace',         &
              action  = 'write')
      write(pf_unit, '(A)') "# :yc, ux_laminar, ux, uy, uz"
      do j = 1, d%nc(2)
        write(pf_unit, '(5ES13.5)') d%yc(j), ufyc(j), ux_ypencil(d%ux_yen(1), j, d%ux_yen(3)) )
      end do
      close(pf_unit)
    end if
    
    return
  end subroutine  Initialize_poiseuille_flow
!===============================================================================
!===============================================================================
!> \brief Initialize Vortex Green flow
!>
!> This subroutine is called locally once.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!> \param[out]    f             flow
!_______________________________________________________________________________
  subroutine  Initialize_vortexgreen_2dflow(ux, uy, uz, p, d)
    use parameters_constant_mod!, only : HALF, ZERO, SIXTEEN, TWO
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in   ) :: d
    real(WP),       intent(inout) :: ux(:, :, :) , &
                                     uy(:, :, :) , &
                                     uz(:, :, :) , &
                                     p (:, :, :)
    real(WP) :: xc, yc
    real(WP) :: xp, yp
    integer(4) :: i, j, ii, jj
!-------------------------------------------------------------------------------
!   Ensure it is in x-pencil
!-------------------------------------------------------------------------------
    if(d%ux_xsz(1) /= d%np(1)) call Print_error_msg("Error, not X-pencil")
!-------------------------------------------------------------------------------
!   ux in x-pencil
!------------------------------------------------------------------------------- 
    do j = 1, d%ux_xsz(2)
      jj = d%ux_xst(2) + j - 1
      yc = d%yc(jj)
      do i = 1, d%ux_xsz(1)
        ii = d%ux_xst(1) + i - 1
        xp = d%h(1) * real(ii - 1, WP)
        ux(i, j, :) =  sin_wp ( xp ) * cos_wp ( yc )
      end do
    end do
!-------------------------------------------------------------------------------
!   uy in x-pencil
!------------------------------------------------------------------------------- 
    do j = 1, d%uy_xsz(2)
      jj = d%uy_xst(2) + j - 1
      yp = d%yp(jj)
      do i = 1, d%uy_xsz(1)
        ii = d%uy_xst(1) + i - 1
        xc = d%h(1) * (real(ii - 1, WP) + HALF)
        uy(i, j, :) = -cos_wp ( xc ) * sin_wp ( yp )
      end do
    end do
!-------------------------------------------------------------------------------
!   uz in x-pencil
!------------------------------------------------------------------------------- 
    uz(:, :, :) =  ZERO
!-------------------------------------------------------------------------------
!   p in x-pencil
!------------------------------------------------------------------------------- 
    do j = 1, d%ps_xsz(2)
      jj = d%ps_xst(2) + j - 1
      yc = d%yc(jj)
      do i = 1, d%ps_xsz(1)
        ii = d%ps_xst(1) + i - 1
        xc = d%h(1) * (real(ii - 1, WP) + HALF)
        p(i, j, :)= ( cos_wp(TWO * xc) + sin(TWO * yc) ) / FOUR
      end do
    end do
    
    return
  end subroutine Initialize_vortexgreen_2dflow
!===============================================================================
!===============================================================================
  subroutine  Validate_TGV2D_error(ux, uy, p, rre, tt, d)
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none

    type(t_domain),    intent(in) :: d
    real(WP),       intent(inout) :: ux(:, :, :) , &
                                     uy(:, :, :) , &
                                     p (:, :, :)
    real(WP),          intent(in) :: rre
    real(WP),          intent(in) :: tt
    integer :: k, i, j
    real(wp) :: uerr, ue, uc, verr, perr
    real(wp) :: xc, yc, xp, yp
    real(wp) :: uerrmax, verrmax, perrmax

    integer :: output_unit
    character( len = 128) :: filename
    logical :: file_exists = .FALSE.
!-------------------------------------------------------------------------------
!   Ensure it is in x-pencil
!-------------------------------------------------------------------------------
    if(d%ux_xsz(1) /= d%np(1)) call Print_error_msg("Error, not X-pencil")
!-------------------------------------------------------------------------------
!   X-pencil : Find Max. error of ux
!-------------------------------------------------------------------------------
    uerr = ZERO
    uerrmax = ZERO
    do k = 1, d%ux_xsz(3)
      do j = 1, d%ux_xsz(2)
        jj = d%ux_xst(2) + j - 1
        yc = d%yc(jj)
        do i = 1, d%ux_xsz(1)
          ii = d%ux_xst(1) + i - 1
          xp = d%h(1) * real(ii - 1, WP)
          uc = ux(i, j, k)
          ue = sin_wp ( xp ) * cos_wp ( yc ) * exp(- TWO * rre * tt)
          uerr = uerr + (uc - ue)**2
          if(dabs(uc - ue) > uerrmax) uerrmax = dabs(uc - ue)
        end do
      end do
    end do
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(uerr,    uerr_work,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    call mpi_allreduce(uerrmax, uerrmax_work, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
    uerr_work = uerr_work / real(d%np(1), wp) / real(d%nc(2), wp) / real(d%nc(3), wp)
    uerr_work = sqrt_wp(uerr_work)
!-------------------------------------------------------------------------------
!   X-pencil : Find Max. error of uy
!-------------------------------------------------------------------------------
    verr = ZERO
    verrmax = ZERO
    do k = 1, d%uy_xsz(3)
      do j = 1, d%uy_xsz(2)
        jj = d%uy_xst(2) + j - 1
        yp = d%yp(jj)
        do i = 1, d%uy_xsz(1)
          ii = d%uy_xst(1) + i - 1
          xc = d%h(1) * (real(ii - 1, WP) + HALF)
          uc = uy(i, j, k)
          ue = - cos_wp ( xc ) * sin_wp ( yp ) * exp(- TWO * f%rre * f%time)
          verr = verr + (uc - ue)**2
          if(dabs(uc - ue) > verrmax) verrmax = dabs(uc - ue)
        end do
      end do
    end do
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(verr,    verr_work,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    call mpi_allreduce(verrmax, verrmax_work, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
    verr_work = verr_work / real(d%nc(1), wp) / real(d%np(2), wp) / real(d%nc(3), wp)
    verr_work = sqrt_wp(verr_work)
!-------------------------------------------------------------------------------
!   X-pencil : Find Max. error of p
!-------------------------------------------------------------------------------
    perr = ZERO
    perrmax = ZERO
    do k = 1, d%ps_xsz(3)
      do j = 1, d%ps_xsz(2)
        jj = d%ps_xst(2) + j - 1
        yc = d%yc(jj)
        do i = 1, d%ps_xsz(1)
          ii = d%ps_xst(1) + i - 1
          xc = d%h(1) * (real(ii - 1, WP) + HALF)
          uc = p(i, j, k)
          ue = ( cos_wp ( TWO * xc ) + sin_wp ( TWO * yc ) ) / FOUR * (exp(- TWO * f%rre * f%time))**2
          perr = perr + (uc - ue)**2
          if(dabs(uc - ue) > perrmax) perrmax = dabs(uc - ue)
        end do
      end do
    end do
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(perr,    perr_work,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    call mpi_allreduce(perrmax, perrmax_work, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
    perr_work = perr_work / real(d%nc(1), wp) / real(d%nc(2), wp) / real(d%nc(3), wp)
    perr_work = sqrt_wp(perr_work)
!-------------------------------------------------------------------------------
!   X-pencil : write data in rank=0
!-------------------------------------------------------------------------------
    if(nrank == 0) then
      filename = 'Validation_TGV2d.dat'
      INQUIRE(FILE = trim(filename), exist = file_exists)
      if(.not.file_exists) then
        open(newunit = output_unit, file = trim(filename), action = "write", status = "new")
        write(output_unit, '(A)') 'Time, SD(u), SD(v), SD(p)'
      else
        open(newunit = output_unit, file = trim(filename), action = "write", status = "old", position="append")
      end if
      write(output_unit, '(1F10.4, 6ES15.7)') tt, uerr_work, verr_work, perr_work, &
            uerrmax_work, verrmax_work, perrmax_work
      close(output_unit)
    end if

    return
  end subroutine
!===============================================================================
!===============================================================================
!> \brief Initialize Vortex Green flow
!>
!> This subroutine is called locally once.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!> \param[out]    f             flow
!_______________________________________________________________________________
  subroutine  Initialize_vortexgreen_3dflow(ux, uy, uz, p, d)
    use parameters_constant_mod, only : HALF, ZERO, SIXTEEN, TWO
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in   ) :: d
    real(WP),       intent(inout) :: ux(:, :, :) , &
                                     uy(:, :, :) , &
                                     uz(:, :, :) , &
                                     p (:, :, :)
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    integer(4) :: i, j, k
!-------------------------------------------------------------------------------
!   Ensure it is in x-pencil
!-------------------------------------------------------------------------------
    if(d%ux_xsz(1) /= d%np(1)) call Print_error_msg("Error, not X-pencil")
!-------------------------------------------------------------------------------
!   ux in x-pencil
!------------------------------------------------------------------------------- 
    do k = 1, d%ux_xsz(3)
      kk = d%ux_xst(3) + k - 1
      zc = d%h(3) * (real(kk - 1, WP) + HALF)
      do j = 1, d%ux_xsz(2)
        jj = d%ux_xst(2) + j - 1
        yc = d%yc(jj)
        do i = 1, d%ux_xsz(1)
          ii = d%ux_xst(1) + i - 1
          xp = d%h(1) * real(ii - 1, WP)
          ux(i, j, k) =  sin_wp ( xp ) * cos_wp ( yc ) * cos_wp ( zc )
        end do
      end do
    end do
!-------------------------------------------------------------------------------
!   uy in x-pencil
!------------------------------------------------------------------------------- 
    do k = 1, d%uy_xsz(3)
      kk = d%uy_xst(3) + k - 1
      zc = d%h(3) * (real(kk - 1, WP) + HALF)
      do j = 1, d%uy_xsz(2)
        jj = d%uy_xst(2) + j - 1
        yp = d%yp(jj)
        do i = 1, d%uy_xsz(1)
          ii = d%uy_xst(1) + i - 1
          xc = d%h(1) * (real(ii - 1, WP) + HALF)
          uy(i, j, k) = -cos_wp ( xc ) * sin_wp ( yp ) * cos_wp ( zc )
        end do
      end do
    end do
!-------------------------------------------------------------------------------
!   uz in x-pencil
!------------------------------------------------------------------------------- 
    uz(:, :, :) =  ZERO
!-------------------------------------------------------------------------------
!   p in x-pencil
!------------------------------------------------------------------------------- 
    do k = 1, d%ps_xsz(3)
      kk = d%ps_xst(3) + k - 1
      zc = d%h(3) * (real(kk - 1, WP) + HALF)
      do j = 1, d%ps_xsz(2)
        jj = d%ps_xst(2) + j - 1
        yc = d%yc(jj)
        do i = 1, d%ps_xsz(1)
          ii = d%ps_xst(1) + i - 1
          xc = d%h(1) * (real(ii - 1, WP) + HALF)
          p(i, j, k)= ( cos_wp( TWO * xc       ) + &
                        cos_wp( TWO * yc       ) ) * &
                      ( cos_wp( TWO * zc + TWO ) ) / SIXTEEN
        end do
      end do
    end do
    
    return
  end subroutine Initialize_vortexgreen_3dflow
!===============================================================================
!===============================================================================
!> \brief Initialize Sine signal for test only
!>
!> This subroutine is called locally once.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!> \param[out]    f             flow
!_______________________________________________________________________________
  subroutine  Initialize_sinetest_flow(ux, uy, uz, p, d)
    use udf_type_mod, only : t_domain, t_flow
    use math_mod, only : sin_wp
    use parameters_constant_mod, only : HALF, ZERO, SIXTEEN, TWO
    
    implicit none
    type(t_domain), intent(in )   :: d
    real(WP),       intent(inout) :: ux(:, :, :), &
                                     uy(:, :, :), &
                                     uz(:, :, :), &
                                     p (:, :, :)

    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    integer(4) :: i, j, k
!-------------------------------------------------------------------------------
!   Ensure it is in x-pencil
!-------------------------------------------------------------------------------
    if(d%ux_xsz(1) /= d%np(1)) call Print_error_msg("Error, not X-pencil")
!-------------------------------------------------------------------------------
!   ux in x-pencil
!------------------------------------------------------------------------------- 
    do k = 1, d%ux_xsz(3)
      kk = d%ux_xst(3) + k - 1
      zc = d%h(3) * (real(kk - 1, WP) + HALF)
      do j = 1, d%ux_xsz(2)
        jj = d%ux_xst(2) + j - 1
        yc = d%yc(jj)
        do i = 1, d%ux_xsz(1)
          ii = d%ux_xst(1) + i - 1
          xp = d%h(1) * real(ii - 1, WP)
          ux(i, j, k) =  sin_wp ( xp ) + sin_wp(yc) + sin_wp(zc)
        end do 
      end do
    end do
!-------------------------------------------------------------------------------
!   uy in x-pencil
!------------------------------------------------------------------------------- 
    do k = 1, d%uy_xsz(3)
      kk = d%uy_xst(3) + k - 1
      zc = d%h(3) * (real(kk - 1, WP) + HALF)
      do i = 1, d%uy_xsz(1)
        ii = d%uy_xst(1) + i - 1
        xc = d%h(1) * (real(ii - 1, WP) + HALF)
        do j = 1, d%uy_xsz(2)
          jj = d%uy_xst(2) + j - 1
          yp = d%yp(jj)
          uy(i, j, k) = sin_wp ( xc ) + sin_wp(yp) + sin_wp(zc)
        end do
      end do
    end do
!-------------------------------------------------------------------------------
!   uz in x-pencil
!------------------------------------------------------------------------------- 
    do j = 1, d%uz_xsz(2)
      jj = d%uz_xst(2) + j - 1
      yc = d%yc(jj)
      do i = 1, d%uz_xsz(1)
        ii = d%uz_xst(1) + i - 1
        xc = d%h(1) * (real(ii - 1, WP) + HALF)
        do k = 1, d%uz_xsz(3)
          kk = d%uz_xst(3) + k - 1
          zp = d%h(3) * real(kk - 1, WP)
          uz(i, j, k) = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zp)
        end do
      end do
    end do
!-------------------------------------------------------------------------------
!   p in x-pencil
!------------------------------------------------------------------------------- 
    do j = 1, d%ps_xsz(2)
      jj = d%ps_xst(2) + j - 1
      yc = d%yc(jj)
      do i = 1, d%ps_xsz(1)
        ii = d%ps_xst(1) + i - 1
        xc = d%h(1) * (real(ii - 1, WP) + HALF)
        do k = 1, d%ps_xsz(3)
          kk = d%ps_xst(3) + k - 1
          zc = d%h(3) * (real(kk - 1, WP) + HALF)
          p(i, j, k) = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
        end do
      end do
    end do
    
    return
  end subroutine Initialize_sinetest_flow
!===============================================================================
!===============================================================================
  subroutine Check_maximum_velocity(ux, uy, uz)
    use precision_mod
    use math_mod
    use mpi_mod
    implicit none

    real(WP), intent(in) :: ux(:, :, :), uy(:, :, :), uz(:, :, :)

    real(WP)   :: u(3), u_work(3)

    u(1) = MAXVAL( abs_wp( ux(:, :, :) ) )
    u(2) = MAXVAL( abs_wp( uy(:, :, :) ) )
    u(3) = MAXVAL( abs_wp( uz(:, :, :) ) )

    call mpi_barrier(MPI_COMM_WORLD, ierror)
    call mpi_allreduce(u, u_work, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)

    if(nrank == 0) then
      Call Print_debug_mid_msg("  The maximum velocities are:")
      write (OUTPUT_UNIT, '(5X, A, 1ES13.5)') 'Umax : ', u_work(1)
      write (OUTPUT_UNIT, '(5X, A, 1ES13.5)') 'Umax : ', u_work(2)
      write (OUTPUT_UNIT, '(5X, A, 1ES13.5)') 'Umax : ', u_work(3)
    end if

    return
  end subroutine

!===============================================================================
!===============================================================================
!> \brief The main code for initializing flow variables
!>
!> This subroutine is called once in \ref Initialize_chapsim.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  none          NA
!_______________________________________________________________________________
  subroutine Initialize_flow_variables( d, f, t )
    use input_general_mod
    use parameters_constant_mod
    use boundary_conditions_mod
    use continuity_eq_mod
    use test_algrithms_mod
    use solver_tools_mod
    use mpi_mod
    implicit none
    type(t_domain), intent(in   ) :: d
    type(t_flow),   intent(inout) :: f
    type(t_thermo), intent(inout) :: t

    interface 
       subroutine Display_vtk_slice(d, str, varnm, vartp, var0, iter)
        use udf_type_mod
        type(t_domain), intent( in ) :: d
        integer(4) :: vartp
        character( len = *), intent( in ) :: str
        character( len = *), intent( in ) :: varnm
        real(WP), intent( in ) :: var0(:, :, :)
        integer(4), intent( in ) :: iter
       end subroutine Display_vtk_slice
    end interface

    if(nrank == 0) call Print_debug_start_msg("Initializing flow and thermal fields ...")

!-------------------------------------------------------------------------------
! to initialize thermal variables 
!-------------------------------------------------------------------------------
    if(nrank == 0) call Print_debug_mid_msg("Initializing thermal field ...")
    if (d%ithermo == 1) then
      call Initialize_thermal_variables (f, t)
    else
      f%dDens(:, :, :) = ONE
      f%mVisc(:, :, :) = ONE
    end if
!-------------------------------------------------------------------------------
! to initialize flow velocity and pressure
!-------------------------------------------------------------------------------
    if(nrank == 0) call Print_debug_mid_msg("Initializing flow field ...")
    if ( (d%icase == ICASE_CHANNEL) .or. &
         (d%icase == ICASE_PIPE) .or. &
         (d%icase == ICASE_ANNUAL) ) then
      call Initialize_poiseuille_flow    (f%qx, f%qy, f%qz, f%pres, f%initNoise, d)
    else if (icase == ICASE_TGV2D) then
      call Initialize_vortexgreen_2dflow (f%qx, f%qy, f%qz, f%pres, d)
    else if (icase == ICASE_TGV3D) then
      call Initialize_vortexgreen_3dflow (f%qx, f%qy, f%qz, f%pres, d)
    else if (icase == ICASE_SINETEST) then
      call Initialize_sinetest_flow      (f%qx, f%qy, f%qz, f%pres, d)
    else 
      if(nrank == 0) call Print_error_msg("No such case defined" )
    end if
!-------------------------------------------------------------------------------
! to initialize pressure correction term
!-------------------------------------------------------------------------------
    f%pcor(:, :, :) = ZERO
!-------------------------------------------------------------------------------
! to check maximum velocity
!-------------------------------------------------------------------------------
    call Check_maximum_velocity(f%qx, f%qy, f%qz)
!-------------------------------------------------------------------------------
! to update mass flux terms 
!-------------------------------------------------------------------------------
    if (d%ithermo == 1) then
      call Calculate_massflux_from_velocity (f, d)
    else
      f%gx(:, :, :) = f%qx(:, :, :)
      f%gy(:, :, :) = f%qy(:, :, :)
      f%gz(:, :, :) = f%qz(:, :, :)
    end if
!-------------------------------------------------------------------------------
! to set up old arrays 
!-------------------------------------------------------------------------------
    f%dDensm1(:, :, :) = f%dDens(:, :, :)
    f%dDensm2(:, :, :) = f%dDens(:, :, :)
!-------------------------------------------------------------------------------
! to write and display the initial fields
!-------------------------------------------------------------------------------
    !call Display_vtk_slice(d, 'xy', 'u', 1, qx)
    !call Display_vtk_slice(d, 'xy', 'v', 2, qy)
    !call Display_vtk_slice(d, 'xy', 'p', 0, pres)

    !call Display_vtk_slice(d, 'yz', 'v', 2, qy)
    !call Display_vtk_slice(d, 'yz', 'w', 3, qz)
    !call Display_vtk_slice(d, 'yz', 'p', 0, pres)

    !call Display_vtk_slice(d, 'zx', 'u', 1, qx)
    !call Display_vtk_slice(d, 'zx', 'w', 3, qz)
    !call Display_vtk_slice(d, 'zx', 'p', 0, pres)

    call Display_vtk_slice(d, 'xy', 'u', 1, f%qx, 0)
    call Display_vtk_slice(d, 'xy', 'v', 2, f%qy, 0)
    call Display_vtk_slice(d, 'xy', 'w', 0, f%pres, 0)

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine
!===============================================================================
!===============================================================================
!> \brief The main code for initializing flow variables
!>
!> This subroutine is called once in \ref Initialize_chapsim.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  
!_______________________________________________________________________________
  subroutine Calculate_RePrGr(f, t, d, iter)
    use input_general_mod
    use input_thermo_mod, only : tpRef0
    use parameters_constant_mod, only : GRAVITY, ONE
    use udf_type_mod, only : t_flow, t_thermo
    implicit none
    type(t_domain), intent(in   ) :: d
    type(t_flow),   intent(inout) :: f
    type(t_thermo), intent(inout) :: t
    integer(4),     intent(in   ) :: iter  
  
    real(WP) :: u0
  
    if(iter < f%nIterIniRen) then
      f%rre = ONE / f%renIni
    else
      f%rre = ONE / f%ren
    end if
  
    if(d%ithermo == 1) then
  
      t%rPrRen = f%rre * tpRef0%k / tpRef0%m / tpRef0%cp
  
      u0 = ONE / f%rre * tpRef0%m / tpRef0%d / t%lenRef
      if (t%igravity == 0) then
        ! no gravity
        f%fgravity = ZERO
      else if (t%igravity == 1 .or. t%igravity == 2 .or. t%igravity == 3 ) then 
        ! flow/gravity same dirction
        f%fgravity =  t%lenRef / u0 / u0 * GRAVITY
      else if (t%igravity == -1 .or. t%igravity == -2 .or. t%igravity == -3 ) then 
        ! flow/gravity opposite dirction
        f%fgravity = -t%lenRef / u0 / u0 * GRAVITY
      else
        ! no gravity
        f%fgravity = ZERO
      end if
  
    end if
  
    return
  end subroutine Calculate_RePrGr

  !===============================================================================
!===============================================================================
!> \brief Initialisation and preprocessing of the flow field
!>
!> This subroutine is called at beginning of the main program
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
  subroutine Initialize_flow_thermal_fields ()
    use input_general_mod
    implicit none

    logical :: itest = .false.
    integer :: i

    do i = 1, nxdomain
      call Allocate_thermoflow_variables (domain(i), flow(i), thermo(i))
      call Calculate_RePrGr(flow(i), thermo(i), 0)
      if (irestart == INITIAL_RANDOM) then
        call Initialize_flow_variables ( domain(i), flow(i), thermo(i) )
        flow(i)%time = ZERO
        thermo(i)%time = ZERO 
      else if (irestart == INITIAL_RESTART) then

      else if (irestart == INITIAL_INTERPL) then

      else
        call Print_error_msg("Error in flow initialisation flag.")
      end if
    end do

  !-------------------------------------------------------------------------------
  ! to test algorithms based on given values.
  !-------------------------------------------------------------------------------
    if(itest) call Test_schemes()
    
    return
  end subroutine Initialize_flow_thermal_fields

end module flow_thermo_initialiasation