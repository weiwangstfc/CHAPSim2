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
module flow_variables_mod
  use save_vars_mod
  implicit none

  private :: Calculate_xbulk_velocity
  private :: Check_maximum_velocity
  private :: Generate_poiseuille_flow_profile
  private :: Initialize_poiseuille_flow
  private :: Initialize_vortexgreen_flow
  private :: Initialize_thermal_variables
  private :: Unify_initial_velocity

  public  :: Allocate_variables
  public  :: Initialize_flow_variables

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
  subroutine Allocate_variables
    use input_general_mod, only : ithermo
    use parameters_constant_mod, only : ZERO, ONE
    use save_vars_mod
    implicit none

    allocate ( flow%qx( 1 : domain%np(1), 1 : domain%nc(2), 1 : domain%nc(3) ) )
    allocate ( flow%qy( 1 : domain%nc(1), 1 : domain%np(2), 1 : domain%nc(3) ) )
    allocate ( flow%qz( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%np(3) ) )
    flow%qx = ZERO
    flow%qy = ZERO
    flow%qz = ZERO

    allocate ( flow%gx ( 1 : domain%np(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    allocate ( flow%gy ( 1 : domain%nc(1), 1 : domain%np(2), 1 : domain%nc(3) )  )
    allocate ( flow%gz ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%np(3) )  )
    flow%gx = ZERO
    flow%gy = ZERO
    flow%gz = ZERO

    allocate ( flow%pres ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    allocate ( flow%pcor ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    flow%pres = ZERO
    flow%pcor = ZERO

    allocate ( flow%dDens ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    allocate ( flow%mVisc ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    flow%dDens = ONE
    flow%mVisc = ONE

    allocate ( flow%dDensm1 ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    allocate ( flow%dDensm2 ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    flow%dDensm1 = ONE
    flow%dDensm2 = ONE

    if(ithermo == 1) then
      allocate ( thermo%dh    ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
      allocate ( thermo%hEnth ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
      allocate ( thermo%kCond ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
      allocate ( thermo%tTemp ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
      thermo%dh    = ZERO
      thermo%hEnth = ZERO
      thermo%kCond = ONE
      thermo%tTemp = ONE
    end if

    return
  end subroutine Allocate_variables
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
  subroutine Initialize_thermal_variables (f, t)
    use input_general_mod, only : tiRef, t0Ref
    use input_thermo_mod, only : tpIni
    implicit none
    type(t_flow),   intent(inout) :: f
    type(t_thermo), intent(inout) :: t
  
    tpIni%t = tiRef / t0Ref
    call tpIni%Refresh_thermal_properties_from_T()

    f%dDens(:, :, :)  = tpIni%d
    f%mVisc(:, :, :)  = tpIni%m

    t%dh    = tpIni%dh
    t%hEnth = tpIni%h
    t%kCond = tpIni%k
    t%tTemp = tpIni%t

    return
  end subroutine Initialize_thermal_variables
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
!> \param[out]    u_laminar     velocity profile along wall-normal direction
!_______________________________________________________________________________
  subroutine Generate_poiseuille_flow_profile(u_laminar, d)
    use parameters_constant_mod, only : ZERO, ONE, ONEPFIVE, TWO, MAXP, TRUNCERR
    use input_general_mod
    use udf_type_mod
    use math_mod
    implicit none

    type(t_domain), intent(in)  :: d
    real(WP),       intent(out) :: u_laminar(:)
    
    real(WP) :: a, b, c, yy, ymax, ymin, umean
    integer(4) :: j

    u_laminar (:) = ZERO

    ymax = d%yp( d%np_geo(2) )
    ymin = d%yp( 1 )
    if (d%case == ICASE_CHANNEL) then
      a = (ymax - ymin) / TWO
      b = ZERO
      c = ONEPFIVE
    else if (d%case == ICASE_PIPE) then
      a = (ymax - ymin)
      b = ZERO
      c = TWO
    else if (d%case == ICASE_ANNUAL) then
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
      u_laminar(j) = ( ONE - ( (yy - b)**2 ) / a / a ) * c
    end do

    ! scale the bulk velocity to be one
    umean = ZERO
    do j = 1, d%nc(2)
      umean = umean + u_laminar(j) * (d%yp(j + 1) - d%yp(j) )
    end do
    umean = umean / (ymax - ymin)

    u_laminar(:) = u_laminar(:) / umean

    ! check the bulk velocity is one
    umean = ZERO
    do j = 1, d%nc(2)
      umean = umean + u_laminar(j) * (d%yp(j + 1) - d%yp(j) )
    end do
    umean = umean / (ymax - ymin)
    if ( abs_wp(umean - ONE) > TRUNCERR) then
      write(*, *) umean
      call Print_error_msg("Error in poiseuille_flow_profile in Subroutine" &
            // "Generate_poiseuille_flow_profile")
    end if

    return
  end subroutine Generate_poiseuille_flow_profile
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
  subroutine Initialize_poiseuille_flow(f, d)
    use random_number_generation_mod
    use parameters_constant_mod, only : ZERO, ONE
    use input_general_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in   ) :: d
    type(t_flow  ), intent(inout) :: f
    
    real(WP), allocatable, dimension(:) :: u_laminar
    integer :: i, j, k
    integer :: seed
    real(WP) :: rd(3)

    ! to get the profile
    allocate ( u_laminar ( d%nc(2) ) ); u_laminar(:) = ZERO
    call Generate_poiseuille_flow_profile ( u_laminar, d )

    f%pres(:, :, :) =  ZERO
    seed = 0 ! real random
    do k = 1, d%nc(3)
      do j = 1, d%nc(2)
        do i = 1, d%nc(1)
          seed = seed + k + j + i ! repeatable random
          call Initialize_random_number ( seed )
          call Generate_rvec_random( -ONE, ONE, 3, rd)
          f%qx(i, j, k) = initNoise * rd(1) + u_laminar (j)
          f%qy(i, j, k) = initNoise * rd(2)
          f%qz(i, j, k) = initNoise * rd(3)
        end do
      end do
    end do

    f%qy(:, 1,       :) = d%ubc(1, 2)
    f%qy(:, d%np(2), :) = d%ubc(2, 2)
    
    deallocate (u_laminar)
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
  subroutine  Initialize_vortexgreen_flow(f, d)
    use parameters_constant_mod, only : HALF, ZERO, SIXTEEN, TWO
    use udf_type_mod
    use math_mod
    implicit none

    type(t_domain), intent(in   ) :: d
    type(t_flow  ), intent(inout) :: f

    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    integer(4) :: i, j, k

    do k = 1, d%nc(3)
      zp = d%h(3) * real(k - 1, WP)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, d%nc(2)
        yp = d%yp(j)
        yc = d%yc(j)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          f%qx(i, j, k) =  sin_wp ( xp ) * cos_wp ( yc ) * cos_wp ( zc )
        end do
      end do
    end do

    do k = 1, d%nc(3)
      zp = d%h(3) * real(k - 1, WP)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, d%np(2)
        yp = d%yp(j)
        yc = d%yc(j)
        do i = 1, d%nc(1)
          xp = d%h(1) * real(i - 1, WP)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          f%qy(i, j, k) = -cos_wp ( xc ) * sin_wp ( yp ) * cos_wp ( zc )
        end do
      end do
    end do

    do k = 1, d%np(3)
      zp = d%h(3) * real(k - 1, WP)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, d%nc(2)
        yp = d%yp(j)
        yc = d%yc(j)
        do i = 1, d%nc(1)
          xp = d%h(1) * real(i - 1, WP)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          f%qz(i, j, k) =  ZERO
        end do
      end do
    end do

    do k = 1, d%nc(3)
      zp = d%h(3) * real(k - 1, WP)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, d%nc(2)
        yp = d%yp(j)
        yc = d%yc(j)
        do i = 1, d%nc(1)
          xp = d%h(1) * real(i - 1, WP)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          f%pres(i, j, k)= ( cos_wp( TWO * xc       ) + &
                             cos_wp( TWO * yc       ) ) * &
                           ( cos_wp( TWO * zc + TWO ) ) / SIXTEEN
        end do
      end do
    end do
    
    return
  end subroutine Initialize_vortexgreen_flow
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
  subroutine  Initialize_sinetest_flow(f, d)
    use udf_type_mod, only: t_domain, t_flow
    use math_mod, only: sin_wp
    use parameters_constant_mod, only : HALF, ZERO, SIXTEEN, TWO
    
    implicit none

    type(t_domain), intent(in )   :: d
    type(t_flow  ), intent(inout) :: f

    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    integer(4) :: i, j, k

    do k = 1, d%nc(3)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, d%nc(2)
        yc = d%yc(j)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          f%qx(i, j, k) =  sin_wp ( xp ) + sin_wp(yc) + sin_wp(zc)
        end do 
      end do
    end do

    do k = 1, d%nc(3)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do i = 1, d%nc(1)
        xc = d%h(1) * (real(i - 1, WP) + HALF)
        do j = 1, d%np(2)
          yp = d%yp(j)
          f%qy(i, j, k) = sin_wp ( xc ) + sin_wp(yp) + sin_wp(zc)
        end do
      end do
    end do

    
    do j = 1, d%nc(2)
      yc = d%yc(j)
      do i = 1, d%nc(1)
        xc = d%h(1) * (real(i - 1, WP) + HALF)
        do k = 1, d%np(3)
          zp = d%h(3) * real(k - 1, WP)
          f%qz(i, j, k) = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zp)
        end do
      end do
    end do

    do j = 1, d%nc(2)
      yc = d%yc(j)
      do i = 1, d%nc(1)
        xc = d%h(1) * (real(i - 1, WP) + HALF)
        do k = 1, d%nc(3)
          zc = d%h(3) * (real(k - 1, WP) + HALF)
          f%pres(i, j, k) = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
        end do
      end do
    end do
    
    return
  end subroutine Initialize_sinetest_flow

  subroutine Unify_initial_velocity(f, d)
    use precision_mod
    implicit none

    type(t_domain), intent(in )   :: d
    type(t_flow  ), intent(inout) :: f

    real(WP)   :: uxa, uya, uza
    integer(4) :: i, j, k

    do j = 1, d%nc(2)
      uxa = sum( f%qx(:, j, :) )
      uza = sum( f%qz(:, j, :) )
      uxa = uxa / real(d%np(1) * d%nc(3), WP)
      uza = uza / real(d%nc(1) * d%np(3), WP)
      f%qx(:, j, :) = f%qx(:, j, :) - uxa
      f%qz(:, j, :) = f%qz(:, j, :) - uza
    end do

    do j = 1, d%np(2)
      uya = sum( f%qy(:, j, :) )
      uya = uya / real(d%nc(1) * d%nc(3), WP)
      f%qy(:, j, :) = f%qy(:, j, :) - uya
    end do

    return
  end subroutine

  subroutine Check_maximum_velocity(f, d, str)
    use precision_mod
    use math_mod
    implicit none

    type(t_domain), intent(in ) :: d
    type(t_flow  ), intent(in ) :: f
    character(*),   intent(in ) :: str 

    real(WP)   :: u(3)
    integer(4) :: i, j, k


    u(1) = MAXVAL( abs_wp( f%qx(:, :, :) ) )
    u(2) = MAXVAL( abs_wp( f%qy(:, :, :) ) )
    u(3) = MAXVAL( abs_wp( f%qz(:, :, :) ) )

    write(*, *) "The maximum velocity (u, v, w) "//str//" is"
    write(*, '(3ES15.7)') u(:)

    return
  end subroutine
!===============================================================================
!===============================================================================
!> \brief The main code for initializing flow variables
!>
!> This subroutine is only for pre-processing/post-processing 2nd order only.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  none          NA
!_______________________________________________________________________________
  subroutine Calculate_xbulk_velocity(f, d, umean)
    use parameters_constant_mod, only : ZERO, HALF
    use operations, only: Get_midp_interpolation
    implicit none

    type(t_domain), intent(in ) :: d
    type(t_flow  ), intent(in ) :: f
    real(WP),       intent(out) :: umean

    real(WP), allocatable :: fi(:), fo(:)
    real(WP) :: lmean
    integer(4) :: i, j, k

    umean = ZERO
    lmean = ZERO
    allocate ( fi( d%nc(2) ) ); fi = ZERO
    allocate ( fo( d%np(2) ) ); fo = ZERO
    do k = 1, d%nc(3)
      do i = 1, d%nc(1)
        fi(:) = f%qx(i, :, k)
        call Get_midp_interpolation( 'y', 'C2P', d, fi(:), fo(:) )
        do j = 1, d%nc(2)
          lmean = lmean + ( d%yp(j+1) - d%yc(j) ) + &
                          ( d%yc(j  ) - d%yp(j) )
          umean = umean + &
                  (fo(j + 1) + fi(j) ) * ( d%yp(j+1) - d%yc(j) ) * HALF + &
                  (fo(j    ) + fi(j) ) * ( d%yc(j  ) - d%yp(j) ) * HALF
        end do
      end do
    end do
    deallocate (fi)
    deallocate (fo)
    umean = umean / real(d%nc(1) * d%nc(3), WP) / lmean 

    write(*, *) "The bulk velocity is ", umean

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
  subroutine Initialize_flow_variables( )
    use save_vars_mod
    use input_general_mod
    use parameters_constant_mod
    use solver_tools_mod
    use boundary_conditions_mod
    use continuity_eq_mod
    use test_algrithms_mod
    implicit none

    real(WP) :: umean
    logical :: itest = .false.

    interface 
       subroutine Display_vtk_slice(d, str, varnm, vartp, var0)
        use udf_type_mod
        type(t_domain), intent( in ) :: d
        integer(4) :: vartp
        character( len = *), intent( in ) :: str
        character( len = *), intent( in ) :: varnm
        real(WP), intent( in ) :: var0(:, :, :)
       end subroutine Display_vtk_slice
    end interface

!-------------------------------------------------------------------------------
! global flow info 
!-------------------------------------------------------------------------------    
    flow%rre = ONE / REN
!-------------------------------------------------------------------------------
! to initialize thermal variables 
!-------------------------------------------------------------------------------
    if (ithermo == 1) then
      call Initialize_thermal_variables (flow, thermo)
    else
      flow%dDens(:, :, :) = ONE
      flow%mVisc(:, :, :) = ONE
    end if
!-------------------------------------------------------------------------------
! to initialize flow velocity and pressure
!-------------------------------------------------------------------------------
    if ( (icase == ICASE_CHANNEL) .or. &
         (icase == ICASE_PIPE) .or. &
         (icase == ICASE_ANNUAL) ) then
      call Initialize_poiseuille_flow  (flow, domain)
    else if (icase == ICASE_TGV) then
      call Initialize_vortexgreen_flow (flow, domain)
    else if (icase == ICASE_SINETEST) then
      call Initialize_sinetest_flow    (flow, domain)
    else 
      call Print_error_msg("No such case defined" )
    end if
!-------------------------------------------------------------------------------
! to initialize pressure correction term
!-------------------------------------------------------------------------------
    flow%pcor(:, :, :) = ZERO
!-------------------------------------------------------------------------------
! to test algorithms based on given values.
!-------------------------------------------------------------------------------
    if(itest) call Test_schemes()
!-------------------------------------------------------------------------------
! to adjust the initial velocity to match a zero mean of x, z direction 
!-------------------------------------------------------------------------------
    call Check_maximum_velocity(flow, domain, "after initialisation")
    call Unify_initial_velocity(flow, domain)
    call Check_maximum_velocity(flow, domain, "after unifying")
!-------------------------------------------------------------------------------
! to apply the b.c. 
!-------------------------------------------------------------------------------
    call Apply_BC_velocity (flow, domain)
!-------------------------------------------------------------------------------
! to unify the bulk velocity to be one for poiseuille flow. 
!-------------------------------------------------------------------------------
    if ( (icase == ICASE_CHANNEL) .or. &
         (icase == ICASE_PIPE) .or. &
         (icase == ICASE_ANNUAL) ) then
      call Calculate_xbulk_velocity(flow, domain, umean)
      flow%qx(:, :, :) = flow%qx(:, :, :) / umean
      call Calculate_xbulk_velocity(flow, domain, umean)
    end if
!-------------------------------------------------------------------------------
! to update mass flux terms 
!-------------------------------------------------------------------------------
    if (ithermo == 1) then
      call Calculate_massflux_from_velocity (flow, domain)
    else
      flow%gx(:, :, :) = flow%qx(:, :, :)
      flow%gy(:, :, :) = flow%qy(:, :, :)
      flow%gz(:, :, :) = flow%qz(:, :, :)
    end if
!-------------------------------------------------------------------------------
! to set up old arrays 
!-------------------------------------------------------------------------------
    flow%dDensm1(:, :, :) = flow%dDens(:, :, :)
    flow%dDensm2(:, :, :) = flow%dDens(:, :, :)
!-------------------------------------------------------------------------------
! to check divergence
!-------------------------------------------------------------------------------
    call Check_divergence(flow, domain, 0)
!-------------------------------------------------------------------------------
! to write and display the initial fields
!-------------------------------------------------------------------------------
    !call Display_vtk_slice(domain, 'xy', 'u', 1, qx)
    !call Display_vtk_slice(domain, 'xy', 'v', 2, qy)
    !call Display_vtk_slice(domain, 'xy', 'p', 0, pres)

    !call Display_vtk_slice(domain, 'yz', 'v', 2, qy)
    !call Display_vtk_slice(domain, 'yz', 'w', 3, qz)
    !call Display_vtk_slice(domain, 'yz', 'p', 0, pres)

    !call Display_vtk_slice(domain, 'zx', 'u', 1, qx)
    !call Display_vtk_slice(domain, 'zx', 'w', 3, qz)
    !call Display_vtk_slice(domain, 'zx', 'p', 0, pres)



    return
  end subroutine

end module flow_variables_mod