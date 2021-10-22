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
!> \file geometry.f90
!>
!> \brief Building up the geometry and mesh information.
!>
!===============================================================================
module geometry_mod
  use type_vars_mod, only : domain
  implicit none

  !private
  private :: Buildup_grid_mapping_1D
  private :: Buildup_neibour_index
  public  :: Initialize_geometry_variables
  
contains
!===============================================================================
!===============================================================================
!> \brief Building up the mesh mapping relation between physical domain and mesh
!>  to a computational domain and mesh.   
!>
!> This subroutine is used locally for 1D only.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     str          string to indicate mapping of cell centre or nodes
!> \param[in]     n            number of mapping points
!> \param[out]    y            the physical coordinate array
!> \param[out]    mp           the mapping relations for 1st and 2nd deriviatives
!_______________________________________________________________________________
  subroutine Buildup_grid_mapping_1D (str, n, y, mp)
!===============================================================================
! Module files
!===============================================================================
    use math_mod
    use input_general_mod
    use parameters_constant_mod
    implicit none
!===============================================================================
! Arguments
!===============================================================================
    character(len = *), intent(in) :: str
    integer(4), intent( in )       :: n
    real(WP), intent( out )        :: y(n)
    real(WP), intent( out )        :: mp(n, 3)
!===============================================================================
! Local Arguments
!===============================================================================    
    integer(4) :: j
    real(WP) :: eta_shift
    real(WP) :: eta_delta
    real(WP) :: alpha, beta, gamma, delta, cc, dd, ee, st1, st2, mm
    real(WP), dimension(n) :: eta

    eta_shift = ZERO
    eta_delta = ONE
    if ( trim( str ) == 'nd' ) then
      eta_shift = ZERO
      eta_delta = ONE / real( n - 1, WP )
    else if ( trim( str ) == 'cl' ) then
      eta_shift = ONE / ( real(n, WP) ) * HALF
      eta_delta = ONE / real( n, WP )
    else 
      call Print_error_msg('Grid stretching location not defined in Subroutine: '// &
      "Buildup_grid_mapping_1D")
    end if

    ! to build up the computational domain \eta \in [0, 1] uniform mesh
    eta(1) = ZERO + eta_shift

    do j = 2, n
      eta(j) = eta(1) + real(j - 1, WP) * eta_delta
    end do

    ! to build up the physical domain y stretching grids based on Eq(53) of Leizet2009JCP
    ! and to build up the derivates based on Eq(53) and (47) in Leizet2009JCP
    gamma = ONE
    delta = ZERO
    if (istret == ISTRET_NO) then
      y(:) = eta(:)
      y(:) = y(:) * (lyt - lyb) + lyb
      mp(:, 1) = ONE
      mp(:, 2) = ONE
      mp(:, 3) = ONE
      return
    else if (istret == ISTRET_CENTRE) then
      gamma = ONE
      delta = ZERO
    else if (istret == ISTRET_2SIDES) then
      gamma = ONE
      delta = HALF
    else if (istret == ISTRET_BOTTOM) then
      gamma = HALF
      delta = HALF
    else if (istret == ISTRET_TOP) then
      gamma = HALF
      delta = ZERO
    else
      call Print_error_msg('Grid stretching flag is not valid in Subroutine: '// &
      "Buildup_grid_mapping_1D")
    end if

    beta = rstret
    alpha =  ( -ONE + sqrt_wp( ONE + FOUR * PI * PI * beta * beta ) ) / beta * HALF

    cc = sqrt_wp( alpha * beta + ONE ) / sqrt_wp( beta )
    dd = cc / sqrt_wp( alpha )
    ee = cc * sqrt_wp( alpha )

    st1 = (ONE   - TWO * delta) / gamma * HALF
    st2 = (THREE - TWO * delta) / gamma * HALF

    do j = 1, n
      mm = PI * (gamma * eta(j) + delta)

      ! y \in [0, 1]
      y(j) = atan_wp ( dd * tan_wp( mm ) ) - &
            atan_wp ( dd * tan_wp( PI * delta) ) + &
            PI * ( heaviside_step( eta(j) - st1 ) + heaviside_step( eta(j) - st2 ) )
      y(j) = ONE / (gamma * ee) * y(j)
      ! y \in [lyb, lyt]
      y(j) = y(j) * (lyt - lyb) + lyb

      ! 1/h'
      mp(j, 1) = (alpha / PI + sin_wp(mm) * sin_wp(mm) / PI / beta)  / (lyt - lyb)

      ! (1/h')^2
      mp(j, 2) = mp(j, 1) * mp(j, 1)

      ! -h"/(h'^3) = 1/h' * [ d(1/h') / d\eta]
      mp(j, 3) = gamma / (lyt - lyb) / beta * sin_wp(TWO * mm) * mp(j, 1)

    end do

    return
  end subroutine Buildup_grid_mapping_1D
!===============================================================================
!===============================================================================
!> \brief Building up the neibouring index of a given index array.   
!>
!> This subroutine is used locally for the bulk part of the grids. The two points
!> near the boundary are not considered except periodic b.c.
!> The neibouring index reduces the repeated calculation of index increase
!> /decrease for a 5-point stencil. 
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n            number of the given index range
!> \param[in]     is_peri      whether the given index array is of periodic b.c.
!> \param[out]    nbr          the neibouring index in order of -2, -1, +1, +2
!_______________________________________________________________________________
  subroutine Buildup_neibour_index(n, is_peri, nbr)
    implicit none
!===============================================================================
! Arguments
!===============================================================================
    integer(4), intent(in)  :: n
    logical,    intent(in)  :: is_peri
    integer(4), intent(out) :: nbr(4, n)
!===============================================================================
! Local arguments
!===============================================================================
    integer(4) :: i
!===============================================================================
! Code
!===============================================================================
    nbr(:, :) = huge(i)

    do i = 3, n
      nbr(1, i) = i - 2
    end do

    do i = 2, n
      nbr(2, i) = i - 1
    end do

    do i = 1, n-1
      nbr(3, i) = i + 1
    end do

    do i = 1, n-2
      nbr(4, i) = i + 2
    end do

    if(is_peri) then
      nbr(1, 1)   = n-1  ! -2, i = 1
      nbr(1, 2)   = n    ! -2, i = 2
      nbr(2, 1)   = n    ! -1, i = 1
      
      nbr(3, n)   = 1    ! +1, i = n
      nbr(4, n)   = 2    ! +2, i = n
      nbr(4, n-1) = 1    ! +2, i = n-1
    end if

    return
  end subroutine
!===============================================================================
!===============================================================================
!> \brief The main code for initializing the geometry, mesh and index.
!>
!> This subroutine is used once in \ref Initialize_chapsim. It builds up the udf
!> domain. Currently only one domain is defined, and it could extend to multiple
!> domain simulations.
!> [mpi] all ranks
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  none          NA
!_______________________________________________________________________________
  subroutine Initialize_geometry_variables ()
!===============================================================================
! Module files
!===============================================================================
    use mpi_mod
    use input_general_mod
    use math_mod
    use parameters_constant_mod, only : ONE, HALF, ZERO, MAXP, MINP, TRUNCERR
    implicit none
!===============================================================================
! Local arguments
!===============================================================================
    integer(4) :: i
    logical    :: dbg = .false.
!===============================================================================
! Code
!===============================================================================
    call Print_debug_start_msg("Initializing domain geometric information ...")
    ! Build up domain info
    domain%case   = icase

    domain%bc(:, 1) = ifbcx(:)
    domain%bc(:, 2) = ifbcy(:)
    domain%bc(:, 3) = ifbcz(:)

    domain%ubc(:, 1) = uxinf(:)
    domain%ubc(:, 2) = uyinf(:)
    domain%ubc(:, 3) = uzinf(:)

    domain%is_periodic(:) = is_periodic(:)

    domain%nc(1) = ncx
    domain%nc(2) = ncy
    domain%nc(3) = ncz

    domain%np_geo(1) = ncx + 1 
    domain%np_geo(2) = ncy + 1
    domain%np_geo(3) = ncz + 1

    do i = 1, 3
      if ( domain%is_periodic(i) ) then
        domain%np(i) = domain%nc(i)
      else 
        domain%np(i) = domain%np_geo(i)
      end if
    end do

    domain%is_stretching(:) = .false.
    if (istret /= ISTRET_NO) domain%is_stretching(2) = .true.
    
    if(domain%is_stretching(2)) then
      domain%h(2) = ONE / real(domain%nc(2), WP)
    else 
      domain%h(2) = (lyt - lyb) / real(domain%nc(2), WP) ! mean dy
    end if
    domain%h(1) = lxx / real(domain%nc(1), WP)
    domain%h(3) = lzz / real(domain%nc(3), WP)
    domain%h2r(:) = ONE / domain%h(:) / domain%h(:)
    domain%h1r(:) = ONE / domain%h(:)

    !build up index sequence for bulk part (no b.c. except periodic)
    allocate ( domain%iNeighb( 4, domain%np(1) ) ); domain%iNeighb =  0
    allocate ( domain%jNeighb( 4, domain%np(2) ) ); domain%jNeighb =  0
    allocate ( domain%kNeighb( 4, domain%np(3) ) ); domain%kNeighb =  0

    call Buildup_neibour_index (domain%np(1), domain%is_periodic(1), domain%iNeighb(:, :) )
    call Buildup_neibour_index (domain%np(2), domain%is_periodic(2), domain%jNeighb(:, :) )
    call Buildup_neibour_index (domain%np(3), domain%is_periodic(3), domain%kNeighb(:, :) )

    ! allocate  variables for mapping physical domain to computational domain
    allocate ( domain%yp( domain%np_geo(2) ) ); domain%yp(:) = ZERO
    allocate ( domain%yc( domain%nc(2) ) ); domain%yc(:) = ZERO

    allocate ( domain%yMappingpt( domain%np_geo(2), 3 ) ); domain%yMappingpt(:, :) = ONE
    allocate ( domain%yMappingcc( domain%nc(2),     3 ) ); domain%yMappingcc(:, :) = ONE

    call Buildup_grid_mapping_1D ('nd', domain%np_geo(2), domain%yp(:), domain%yMappingPt(:, :))
    call Buildup_grid_mapping_1D ('cl', domain%nc(2),     domain%yc(:), domain%yMappingcc(:, :))

    ! print out for debugging
    if(dbg) then
      do i = 1, domain%np_geo(2)
        write(*, '(I5, 1F8.4)') i, domain%yp(i)
      end do
    end if

    call Print_debug_end_msg
    return
  end subroutine  Initialize_geometry_variables
end module geometry_mod

