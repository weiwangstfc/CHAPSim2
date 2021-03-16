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
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------
!===============================================================================
!> \file flow_initialisation.f90
!>
!> \brief Define and initialise flow and thermal variables.
!>
!===============================================================================
module flow_variables_mod
  use precision_mod
  implicit none

  real(WP), save, allocatable, dimension(:, :, :) :: qx, qy, qz
  real(WP), save, allocatable, dimension(:, :, :) :: gx, gy, gz
  real(WP), save, allocatable, dimension(:, :, :) :: pres
  real(WP), save, allocatable, dimension(:, :, :) :: pcor

  real(WP), save, allocatable, dimension(:, :, :) :: dDens
  real(WP), save, allocatable, dimension(:, :, :) :: mVisc

  real(WP), save, allocatable, dimension(:, :, :) :: dh
  real(WP), save, allocatable, dimension(:, :, :) :: hEnth
  real(WP), save, allocatable, dimension(:, :, :) :: kCond
  real(WP), save, allocatable, dimension(:, :, :) :: tTemp
  
  private :: Generate_poiseuille_flow_profile
  private :: Initialize_poiseuille_flow
  private :: Initialize_vortexgreen_flow
  private :: Initialize_thermal_variables

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
    use geometry_mod
    use parameters_constant_mod, only : ZERO, ONE
    implicit none

    allocate ( qx ( 1 : domain%np(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    allocate ( qy ( 1 : domain%nc(1), 1 : domain%np(2), 1 : domain%nc(3) )  )
    allocate ( qz ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%np(3) )  )
    qx = ZERO
    qy = ZERO
    qz = ZERO

    allocate ( gx ( 1 : domain%np(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    allocate ( gy ( 1 : domain%nc(1), 1 : domain%np(2), 1 : domain%nc(3) )  )
    allocate ( gz ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%np(3) )  )
    gx = ZERO
    gy = ZERO
    gz = ZERO

    allocate ( pres ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    allocate ( pcor ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    pres = ZERO
    pcor = ZERO

    allocate ( dDens ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    allocate ( mVisc ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
    dDens = ONE
    mVisc = ONE


    if(ithermo == 1) then
      allocate ( dh    ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
      allocate ( hEnth ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
      allocate ( kCond ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
      allocate ( tTemp ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  )
      dh    = ZERO
      hEnth = ZERO
      kCond = ONE
      tTemp = ONE
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
!> \param[inout]  d_dummy       density
!> \param[inout]  m_dummy       dynamic viscousity
!> \param[inout]  dh_dummy      density * enthalpy
!> \param[inout]  h_dummy       enthalpy
!> \param[inout]  k_dummy       thermal conductivity
!> \param[inout]  t_dummy       temperature
!_______________________________________________________________________________
  subroutine Initialize_thermal_variables (d_dummy, m_dummy, dh_dummy, h_dummy, k_dummy, t_dummy)
    use input_general_mod, only : tiRef, t0Ref
    use input_thermo_mod, only : tpIni
    implicit none

    real(WP), intent(inout) :: d_dummy(:, :, :)
    real(WP), intent(inout) :: m_dummy(:, :, :)
    real(WP), intent(inout) :: dh_dummy(:, :, :)
    real(WP), intent(inout) :: h_dummy(:, :, :)
    real(WP), intent(inout) :: k_dummy(:, :, :)
    real(WP), intent(inout) :: t_dummy(:, :, :)
  
    tpIni%t = tiRef / t0Ref
    call tpIni%Refresh_thermal_properties_from_T()

    d_dummy(:, :, :)  = tpIni%d
    m_dummy(:, :, :)  = tpIni%m

    dh_dummy(:, :, :) = tpIni%dh
    h_dummy(:, :, :)  = tpIni%h
    k_dummy(:, :, :)  = tpIni%k
    t_dummy(:, :, :)  = tpIni%t
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
!> \param[out]    ux_dummy      velocity in the streamwise direction
!> \param[out]    uy_dummy      velocity in the wall-normal direction
!> \param[out]    uz_dummy      velocity in the spanwise direction
!> \param[out]    pres_dummy    pressure
!_______________________________________________________________________________
  subroutine Initialize_poiseuille_flow(ux_dummy, uy_dummy, uz_dummy, pres_dummy, d)
    use random_number_generation_mod
    use parameters_constant_mod, only : ZERO, ONE
    use input_general_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in) :: d
    real(WP), intent(out) :: ux_dummy(:, :, :)
    real(WP), intent(out) :: uy_dummy(:, :, :)
    real(WP), intent(out) :: uz_dummy(:, :, :)
    real(WP), intent(out) :: pres_dummy(:, :, :)
    
    real(WP), allocatable, dimension(:) :: u_laminar
    integer :: i, j, k
    integer :: seed
    real(WP) :: rd(3)

    ! to get the profile
    allocate ( u_laminar ( d%nc(2) ) ); u_laminar(:) = ZERO
    call Generate_poiseuille_flow_profile ( u_laminar, d )

    pres_dummy(:, :, :) =  ZERO
    seed = 0 ! real random
    do k = 1, d%nc(3)
      do j = 1, d%nc(2)
        do i = 1, d%nc(1)
          seed = seed + k + j + i ! repeatable random
          call Initialize_random_number ( seed )
          call Generate_rvec_random( -ONE, ONE, 3, rd)
          ux_dummy(i, j, k) = initNoise * rd(1) + u_laminar (j)
          uy_dummy(i, j, k) = initNoise * rd(2)
          uz_dummy(i, j, k) = initNoise * rd(3)
        end do
      end do
    end do

    uy_dummy(:, 1,       :) = d%ubc(1, 2)
    uy_dummy(:, d%np(2), :) = d%ubc(2, 2)
    
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
!> \param[out]    ux_dummy      velocity in the streamwise direction
!> \param[out]    uy_dummy      velocity in the wall-normal direction
!> \param[out]    uz_dummy      velocity in the spanwise direction
!> \param[out]    pres_dummy    pressure
!_______________________________________________________________________________
  subroutine  Initialize_vortexgreen_flow(ux_dummy, uy_dummy, uz_dummy, pres_dummy, d)
    use parameters_constant_mod, only : HALF, ZERO, SIXTEEN, TWO
    use udf_type_mod
    use math_mod
    implicit none

    type(t_domain), intent(in) :: d
    real(WP), intent(out) :: ux_dummy(:, :, :)
    real(WP), intent(out) :: uy_dummy(:, :, :)
    real(WP), intent(out) :: uz_dummy(:, :, :)
    real(WP), intent(out) :: pres_dummy(:, :, :)

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
          ux_dummy(i, j, k) =  sin_wp ( xp ) * cos_wp ( yc ) * cos_wp ( zc )
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
          uy_dummy(i, j, k) = -cos_wp ( xc ) * sin_wp ( yp ) * cos_wp ( zc )
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
          uz_dummy(i, j, k) =  ZERO
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
          pres_dummy(i, j, k)= ( cos_wp( TWO * xc       ) + cos_wp( TWO * yc       ) ) * &
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
!> \param[out]    ux_dummy      velocity in the streamwise direction
!> \param[out]    uy_dummy      velocity in the wall-normal direction
!> \param[out]    uz_dummy      velocity in the spanwise direction
!> \param[out]    pres_dummy    pressure
!_______________________________________________________________________________
  subroutine  Initialize_sinetest_flow(ux_dummy,   &
                                       uy_dummy,   &
                                       uz_dummy,   &
                                       pres_dummy, &
                                       d)
    use parameters_constant_mod, only : HALF, ZERO, SIXTEEN, TWO
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in) :: d
    real(WP), intent(out) :: ux_dummy(:, :, :)
    real(WP), intent(out) :: uy_dummy(:, :, :)
    real(WP), intent(out) :: uz_dummy(:, :, :)
    real(WP), intent(out) :: pres_dummy(:, :, :)
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    integer(4) :: i, j, k

    do k = 1, d%nc(3)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, d%nc(2)
        yc = d%yc(j)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          ux_dummy(i, j, k) =  sin_wp ( xp ) + sin_wp(yc) + sin_wp(zc)
        end do 
      end do
    end do

    do k = 1, d%nc(3)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do i = 1, d%nc(1)
        xc = d%h(1) * (real(i - 1, WP) + HALF)
        do j = 1, d%np(2)
          yp = d%yp(j)
          uy_dummy(i, j, k) = sin_wp ( xc ) + sin_wp(yp) + sin_wp(zc)
        end do
      end do
    end do

    
    do j = 1, d%nc(2)
      yc = d%yc(j)
      do i = 1, d%nc(1)
        xc = d%h(1) * (real(i - 1, WP) + HALF)
        do k = 1, d%np(3)
          zp = d%h(3) * real(k - 1, WP)
          uz_dummy(i, j, k) = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zp)
        end do
      end do
    end do

    do j = 1, d%nc(2)
      yc = d%yc(j)
      do i = 1, d%nc(1)
        xc = d%h(1) * (real(i - 1, WP) + HALF)
        do k = 1, d%nc(3)
          zc = d%h(3) * (real(k - 1, WP) + HALF)
          pres_dummy(i, j, k) = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
        end do
      end do
    end do
    
    return
  end subroutine Initialize_sinetest_flow
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
    use geometry_mod
    use input_general_mod
    use parameters_constant_mod
    use solver_tools_mod
    implicit none

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
! to initialize thermal variables 
!-------------------------------------------------------------------------------
    if (ithermo == 1) then
      call Initialize_thermal_variables (dDens, mVisc, dh, hEnth, kCond, tTemp)
    else
      dDens(:, :, :) = ONE
      mVisc(:, :, :) = ONE
    end if
!-------------------------------------------------------------------------------
! to initialize flow velocity and pressure
!-------------------------------------------------------------------------------
    if ( (icase == ICASE_CHANNEL) .or. &
         (icase == ICASE_PIPE) .or. &
         (icase == ICASE_ANNUAL) ) then

      call Initialize_poiseuille_flow (qx, qy, qz, pres, domain)

    else if (icase == ICASE_TGV) then
      
      call Initialize_vortexgreen_flow (qx, qy, qz, pres, domain)
    else if (icase == ICASE_SINETEST) then
      call Initialize_sinetest_flow (qx, qy, qz, pres, domain)
    else 
      call Print_error_msg("No such case defined" )
    end if
!-------------------------------------------------------------------------------
! to initialize pressure correction term
!-------------------------------------------------------------------------------
    pcor(:, :, :) = ZERO
!-------------------------------------------------------------------------------
! to update mass flux terms 
!-------------------------------------------------------------------------------
    if (ithermo == 1) then
      call Calculate_massflux_from_velocity (qx, qy, qz, gx, gy, gz, dDens, d)
    else
      gx(:, :, :) = qx(:, :, :)
      gy(:, :, :) = qy(:, :, :)
      gz(:, :, :) = qz(:, :, :)
    end if
!-------------------------------------------------------------------------------
! to set up old arrays 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! to write and display the initial fields
!-------------------------------------------------------------------------------
    !call Display_vtk_slice(domain, 'xy', 'u', 1, qx)
    !call Display_vtk_slice(domain, 'xy', 'v', 2, qy)
    call Display_vtk_slice(domain, 'xy', 'p', 0, pres)

    !call Display_vtk_slice(domain, 'yz', 'v', 2, qy)
    !call Display_vtk_slice(domain, 'yz', 'w', 3, qz)
    call Display_vtk_slice(domain, 'yz', 'p', 0, pres)

    !call Display_vtk_slice(domain, 'zx', 'u', 1, qx)
    !call Display_vtk_slice(domain, 'zx', 'w', 3, qz)
    call Display_vtk_slice(domain, 'zx', 'p', 0, pres)



    return
  end subroutine

end module flow_variables_mod