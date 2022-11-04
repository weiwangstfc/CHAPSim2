!----------------------------------------------------------------------------------------------------------
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

!----------------------------------------------------------------------------------------------------------
!==========================================================================================================
!> \file geometry.f90
!>
!> \brief Building up the geometry and mesh information.
!>
!==========================================================================================================
module geometry_mod
  use vars_df_mod, only : domain
  use precision_mod
  implicit none

  integer, save :: ndm = 0
  
  real(WP) :: alpha, beta, gamma, delta
  !private
  private :: Buildup_grid_mapping_1D
  public  :: Buildup_geometry_mesh_info
  
contains
!==========================================================================================================
!==========================================================================================================
!> \brief Building up the mesh mapping relation between physical domain and mesh
!>  to a computational domain and mesh.   
!>
!> This subroutine is used locally for 1D only.
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     str          string to indicate mapping of cell centre or nodes
!> \param[in]     n            number of mapping points
!> \param[out]    y            the physical coordinate array
!> \param[out]    mp           the mapping relations for 1st and 2nd deriviatives
!_______________________________________________________________________________
  subroutine Buildup_grid_mapping_1D (str, n, dm, y, mp)
    use math_mod
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    character(len = *), intent(in)   :: str
    integer,            intent(in)   :: n
    type(t_domain),     intent(in)   :: dm
    real(WP), intent( out )        :: y(n)
    real(WP), intent( out )        :: mp(n, 3)
    integer :: j
    real(WP) :: eta_shift
    real(WP) :: eta_delta
    
    real(WP) :: cc, dd, ee, st1, st2, mm
    real(WP), dimension(n) :: eta

    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    eta_shift = ZERO
    eta_delta = ONE
    if ( trim( str ) == 'nd' ) then
      eta_shift = ZERO
      eta_delta = ONE / real( n - 1, WP )
    else if ( trim( str ) == 'cl' ) then
      eta_shift = ONE / ( real(n, WP) ) * HALF
      eta_delta = ONE / real( n, WP )
    else 
      call Print_error_msg('Grid stretching location not defined.')
    end if
    
    !----------------------------------------------------------------------------------------------------------
    ! to build up the computational domain \eta \in [0, 1] uniform mesh
    !----------------------------------------------------------------------------------------------------------
    eta(1) = ZERO + eta_shift

    do j = 2, n
      eta(j) = eta(1) + real(j - 1, WP) * eta_delta
    end do
    !----------------------------------------------------------------------------------------------------------
    ! to build up the physical domain y stretching grids based on Eq(53) of Leizet2009JCP
    ! and to build up the derivates based on Eq(53) and (47) in Leizet2009JCP
    !----------------------------------------------------------------------------------------------------------
    gamma = ONE
    delta = ZERO
    if (dm%istret == ISTRET_NO) then
      do j = 1, n
        y(j) = eta(j)
        y(j) = y(j) * (dm%lyt - dm%lyb) + dm%lyb
        mp(j, 1) = ONE
        mp(j, 2) = ONE
        mp(j, 3) = ONE
      end do
      return
    else if (dm%istret == ISTRET_CENTRE) then
      gamma = ONE
      delta = ZERO
    else if (dm%istret == ISTRET_2SIDES) then
      gamma = ONE
      delta = HALF
    else if (dm%istret == ISTRET_BOTTOM) then
      gamma = HALF
      delta = HALF
    else if (dm%istret == ISTRET_TOP) then
      gamma = HALF
      delta = ZERO
    else
      call Print_error_msg('Grid stretching flag is not valid.')
    end if

    beta = dm%rstret
    alpha =  ( -ONE + sqrt_wp( ONE + FOUR * PI * PI * beta * beta ) ) / beta * HALF

    cc = sqrt_wp( alpha * beta + ONE ) / sqrt_wp( beta )
    dd = cc / sqrt_wp( alpha )
    ee = cc * sqrt_wp( alpha )

    st1 = (ONE   - TWO * delta) / gamma * HALF
    st2 = (THREE - TWO * delta) / gamma * HALF

    do j = 1, n
      mm = PI * (gamma * eta(j) + delta)
      !----------------------------------------------------------------------------------------------------------
      ! y \in [0, 1]
      !----------------------------------------------------------------------------------------------------------
      y(j) = atan_wp ( dd * tan_wp( mm ) ) - &
            atan_wp ( dd * tan_wp( PI * delta) ) + &
            PI * ( heaviside_step( eta(j) - st1 ) + heaviside_step( eta(j) - st2 ) )
      y(j) = ONE / (gamma * ee) * y(j)
      !----------------------------------------------------------------------------------------------------------
      ! y \in [lyb, lyt]
      !----------------------------------------------------------------------------------------------------------
      y(j) = y(j) * (dm%lyt - dm%lyb) + dm%lyb
      !----------------------------------------------------------------------------------------------------------
      ! 1/h'
      !----------------------------------------------------------------------------------------------------------
      mp(j, 1) = (alpha / PI + sin_wp(mm) * sin_wp(mm) / PI / beta)  / (dm%lyt - dm%lyb)
      !----------------------------------------------------------------------------------------------------------
      ! (1/h')^2
      !----------------------------------------------------------------------------------------------------------
      mp(j, 2) = mp(j, 1) * mp(j, 1)
      !----------------------------------------------------------------------------------------------------------
      ! -h"/(h'^3) = 1/h' * [ d(1/h') / d\eta]
      !----------------------------------------------------------------------------------------------------------
      mp(j, 3) = gamma / (dm%lyt - dm%lyb) / beta * sin_wp(TWO * mm) * mp(j, 1)

    end do

    return
  end subroutine Buildup_grid_mapping_1D
!==========================================================================================================
  subroutine Buildup_geometry_mesh_info (dm)
    use mpi_mod
    use math_mod
    use parameters_constant_mod
    use udf_type_mod
    use typeconvert_mod
    use mpi_mod
    use wtformat_mod
    implicit none

    type(t_domain), intent(inout) :: dm

    integer    :: i, j, k
    integer    :: outputunit
    logical    :: file_exists = .FALSE.
    character( len = 128) :: filename
    
    if(nrank == 0) call Print_debug_start_msg("Initializing domain geometric ...")

    !----------------------------------------------------------------------------------------------------------
    ! set up node number in geometry domain
    !----------------------------------------------------------------------------------------------------------
    do i = 1, NDIM
      dm%np_geo(i) = dm%nc(i) + 1 
    end do
    !----------------------------------------------------------------------------------------------------------
    ! set up node number in computational domain
    !----------------------------------------------------------------------------------------------------------
    do i = 1, NDIM
      if ( dm%is_periodic(i) ) then
        dm%np(i) = dm%nc(i)
      else 
        dm%np(i) = dm%np_geo(i)
      end if
    end do
    !----------------------------------------------------------------------------------------------------------
    ! set dx, dz for uniform grids
    !----------------------------------------------------------------------------------------------------------
    dm%h(1) = dm%lxx / real(dm%nc(1), WP)
    dm%h(3) = dm%lzz / real(dm%nc(3), WP)
    !----------------------------------------------------------------------------------------------------------
    ! allocate  variables for mapping physical domain to computational domain
    !----------------------------------------------------------------------------------------------------------
    allocate ( dm%yp( dm%np_geo(2) ) ) ! yp(1:np_geo)
    allocate ( dm%yc( dm%nc    (2) ) ) ! yc(1:nc)
    dm%yp(:) = ZERO
    dm%yc(:) = ZERO

    if(dm%is_stretching(2)) then
      allocate ( dm%yMappingpt( dm%np_geo(2), 3 ) )
      allocate ( dm%yMappingcc( dm%nc    (2), 3 ) )
      dm%yMappingpt(:, :) = ONE
      dm%yMappingcc(:, :) = ONE
      dm%h(2) = ONE / real(dm%nc(2), WP)
      call Buildup_grid_mapping_1D ('nd', dm%np_geo(2), dm, dm%yp(:), dm%yMappingpt(:, :))
      call Buildup_grid_mapping_1D ('cl', dm%nc(2),     dm, dm%yc(:), dm%yMappingcc(:, :))
    else
      dm%h(2) = (dm%lyt - dm%lyb) / real(dm%nc(2), WP)
      do i = 1, dm%np_geo(2)
        dm%yp(i) = real(i - 1, WP) * dm%h(2) + dm%lyb
      end do
      do i = 1, dm%nc(2)
        dm%yc(i) = real(i - 1, WP) * dm%h(2) + dm%h(2) * HALF + dm%lyb
      end do
    end if
    !----------------------------------------------------------------------------------------------------------
    ! set 1/dx, 1/(dx)^2
    !----------------------------------------------------------------------------------------------------------
    do i = 1, NDIM
      dm%h2r(i) = ONE / (dm%h(i) * dm%h(i))
      dm%h1r(i) = ONE / dm%h(i)
    end do
    !----------------------------------------------------------------------------------------------------------
    ! print out data for debugging
    !----------------------------------------------------------------------------------------------------------
#ifdef DEBUG_STEPS
    if(nrank == 0) then
      open(221, file = 'mesh_yp.dat')
      do i = 1, dm%np_geo(2)
        write (221, *) i, dm%yp(i)
      end do
    end if
#endif

    ! print out for postprocessing
    ndm = ndm + 1
    if(nrank == 0) then
      filename = 'display_mesh_domain'//trim(int2str(ndm))//'.vtk'
      INQUIRE(FILE = trim(filename), exist = file_exists)
      if(.not.file_exists) then
        open(newunit = outputunit, file = trim(filename), action = "write", status = "new")
        write(outputunit, '(A)') '# vtk DataFile Version 4.0'
        write(outputunit, '(A)') 'vtk mesh single precision'
        write(outputunit, '(A)') 'ASCII'
        write(outputunit, '(A)') 'DATASET RECTILINEAR_GRID'
        write(outputunit, '(A, 3I10.1)') 'DIMENSIONS', dm%np_geo(1:3)
        write(outputunit, '(A, 1I10.1, A)') 'X_COORDINATES', dm%np_geo(1), 'double'
        write(outputunit, *) (dm%h(1) * real(i - 1, WP), i = 1, dm%np_geo(1))
        write(outputunit, '(A, 1I10.1, A)') 'Y_COORDINATES', dm%np_geo(2), 'double'
        write(outputunit, *) (dm%yp(j), j = 1, dm%np_geo(2))
        write(outputunit, '(A, 1I10.1, A)') 'Z_COORDINATES', dm%np_geo(3), 'double'
        write(outputunit, *) (dm%h(3) * real(k - 1, WP), k = 1, dm%np_geo(3))
        close(outputunit)
      end if
    end if

    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine  Buildup_geometry_mesh_info
!==========================================================================================================
  subroutine Buildup_geometry_mesh_info_all_domains
    use mpi_mod
    use vars_df_mod
    implicit none

    integer :: i

    do i = 1, nxdomain
      call Buildup_geometry_mesh_info(domain(i)) 
    end do

  end subroutine 
end module geometry_mod

