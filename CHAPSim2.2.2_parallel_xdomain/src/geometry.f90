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
!> \brief Building up the geometry and mesh inwrtfmt1ion.
!>
!===============================================================================
module geometry_mod
  use vars_df_mod, only : domain
  implicit none

  !private
  private :: Buildup_grid_mapping_1D
  private :: Buildup_npneibour_index
  private :: Buildup_ncneibour_index
  public  :: Buildup_geometry_mesh_info
  
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
  subroutine Buildup_grid_mapping_1D (str, n, y, dm, mp)
    use math_mod
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    character(len = *), intent(in) :: str
    integer, intent( in )          :: n
    type(t_domain), intent(in)     :: dm
    real(WP), intent( out )        :: y(n)
    real(WP), intent( out )        :: mp(n, 3)
    integer :: j
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
    if (dm%istret == ISTRET_NO) then
      y(:) = eta(:)
      y(:) = y(:) * (dm%lyt - dm%lyb) + dm%lyb
      mp(:, 1) = ONE
      mp(:, 2) = ONE
      mp(:, 3) = ONE
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
      call Print_error_msg('Grid stretching flag is not valid in Subroutine: '// &
      "Buildup_grid_mapping_1D")
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

      ! y \in [0, 1]
      y(j) = atan_wp ( dd * tan_wp( mm ) ) - &
            atan_wp ( dd * tan_wp( PI * delta) ) + &
            PI * ( heaviside_step( eta(j) - st1 ) + heaviside_step( eta(j) - st2 ) )
      y(j) = ONE / (gamma * ee) * y(j)
      ! y \in [lyb, lyt]
      y(j) = y(j) * (dm%lyt - dm%lyb) + dm%lyb

      ! 1/h'
      mp(j, 1) = (alpha / PI + sin_wp(mm) * sin_wp(mm) / PI / beta)  / (dm%lyt - dm%lyb)

      ! (1/h')^2
      mp(j, 2) = mp(j, 1) * mp(j, 1)

      ! -h"/(h'^3) = 1/h' * [ d(1/h') / d\eta]
      mp(j, 3) = gamma / (dm%lyt - dm%lyb) / beta * sin_wp(TWO * mm) * mp(j, 1)

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
  subroutine Buildup_npneibour_index(n, ibc, nbr, str)
    use parameters_constant_mod
    implicit none
    integer, intent(in)  :: n ! np
    integer, intent(in)  :: ibc(2)
    integer, intent(inout) :: nbr(4, 4)
    character(3), intent(in) :: str
    integer :: i, j
    integer :: npmax
!-------------------------------------------------------------------------------
! nbr(i, j)
!   i = 1, 2, 3, 4 : left two index, right two index
!   j = 1, 2, 3, 4 : i - 2, i - 1, i + 1, i + 2
!-------------------------------------------------------------------------------
!   for non-periodic:
!   ---(-1')---(0')---(|1')---(2')---(3')---(i')---(n'-1)---(n'|)--(n'+1)---(n'+2)---
!   for periodic:
!   ---(-1')---(0')---(|1')---(2')---(3')---(i')---(n'-1)---(n')--(n'+1|)---(n'+2)---
!-------------------------------------------------------------------------------
    nbr(:, :) = huge(i)
    if(str == 'p2p') then
      npmax = n
    else if(str == 'p2c')then
      if(ibc(1) == IBC_PERIODIC .or. ibc(2) == IBC_PERIODIC) then
        npmax = n
      else
        npmax = n + 1
      end if
    else
    end if
!-------------------------------------------------------------------------------
! left
!-------------------------------------------------------------------------------
    do i = 1, 2
      nbr(i, 1) = i - 2
      nbr(i, 2) = i - 1
      nbr(i, 3) = i + 1
      nbr(i, 4) = i + 2
    end do

    do i = 1, 2
      do j = 1, 4
        if(ibc(1) == IBC_PERIODIC) then
          if( nbr(i, j) == -1 ) nbr(i, j) = npmax - 1
          if( nbr(i, j) == 0  ) nbr(i, j) = npmax
        else if (ibc(1) == IBC_SYMMETRIC .or. ibc(1) == IBC_ASYMMETRIC) then
          if( nbr(i, j) == -1 ) nbr(i, j) = 3
          if( nbr(i, j) == 0  ) nbr(i, j) = 2
        else
        end if
      end do
    end do
!-------------------------------------------------------------------------------
! right
!-------------------------------------------------------------------------------
    do i = 3, 4
      if (i == 3) j = npmax - 1
      if (i == 4) j = npmax 
      nbr(i, 1) = j - 2
      nbr(i, 2) = j - 1
      nbr(i, 3) = j + 1
      nbr(i, 4) = j + 2
    end do
    
    do i = 3, 4
      do j = 1, 4
        if(ibc(2) == IBC_PERIODIC) then
          if( nbr(i, j) == npmax + 1 ) nbr(i, j) = 1
          if( nbr(i, j) == npmax + 2 ) nbr(i, j) = 2
        else if (ibc(2) == IBC_SYMMETRIC .or. ibc(2) == IBC_ASYMMETRIC) then
          if( nbr(i, j) == npmax + 1 ) nbr(i, j) = npmax - 1
          if( nbr(i, j) == npmax + 2 ) nbr(i, j) = npmax - 2
        else
        end if
      end do
    end do

    return
  end subroutine
!_______________________________________________________________________________
  subroutine Buildup_ncneibour_index(n, ibc, nbr, str)
    use parameters_constant_mod
    implicit none
    integer, intent(in)  :: n ! nc
    integer, intent(in)  :: ibc(2)
    integer, intent(inout) :: nbr(4, 4)
    character(3), intent(in) :: str
    integer :: i, j
    integer :: ncmax
!-------------------------------------------------------------------------------
! nbr(i, j)
!   i = 1, 2, 3, 4 : left two index, right two index
!   j = 1, 2, 3, 4 : i - 2, i - 1, i + 1, i + 2
!-------------------------------------------------------------------------------
!   for non-periodic:
!   ---(-1)---(0)--|--(1)---(2)---(3)---(i)---(n-1)---(n)--|--(n+1)---(n+2)---
!   for periodic:
!   ---(-1)---(0)--|--(1)---(2)---(3)---(i)---(n-1)---(n)--|--(n+1)---(n+2)---
!-------------------------------------------------------------------------------
    nbr(:, :) = huge(i)
    if(str == 'c2c') then
      ncmax = n
    else if(str == 'c2p')then
      if(ibc(1) == IBC_PERIODIC .or. ibc(2) == IBC_PERIODIC) then
        ncmax = n
      else
        ncmax = n - 1
      end if
    else
    end if
!-------------------------------------------------------------------------------
! left
!-------------------------------------------------------------------------------
    do i = 1, 2
      nbr(i, 1) = i - 2
      nbr(i, 2) = i - 1
      nbr(i, 3) = i + 1
      nbr(i, 4) = i + 2
    end do

    do i = 1, 2
      do j = 1, 4
        if(ibc(1) == IBC_PERIODIC) then
          if( nbr(i, j) == -1 ) nbr(i, j) = ncmax - 1
          if( nbr(i, j) == 0  ) nbr(i, j) = ncmax
        else if (ibc(1) == IBC_SYMMETRIC .or. ibc(1) == IBC_ASYMMETRIC) then
          if( nbr(i, j) == -1 ) nbr(i, j) = 2
          if( nbr(i, j) == 0  ) nbr(i, j) = 1
        else
        end if
      end do
    end do
!-------------------------------------------------------------------------------
! right
!-------------------------------------------------------------------------------
    do i = 3, 4
      if (i == 3) j = ncmax - 1
      if (i == 4) j = ncmax 
      nbr(i, 1) = j - 2
      nbr(i, 2) = j - 1
      nbr(i, 3) = j + 1
      nbr(i, 4) = j + 2
    end do
    
    do i = 3, 4
      do j = 1, 4
        if(ibc(2) == IBC_PERIODIC) then
          if( nbr(i, j) == ncmax + 1 ) nbr(i, j) = 1
          if( nbr(i, j) == ncmax + 2 ) nbr(i, j) = 2
        else if (ibc(2) == IBC_SYMMETRIC .or. ibc(2) == IBC_ASYMMETRIC) then
          if( nbr(i, j) == ncmax + 1 ) nbr(i, j) = ncmax
          if( nbr(i, j) == ncmax + 2 ) nbr(i, j) = ncmax - 1
        else
        end if
      end do
    end do

    return
  end subroutine
!===============================================================================
  subroutine Buildup_geometry_mesh_info (dm)
    use mpi_mod
    use math_mod
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm
    integer :: i, j
    logical    :: dbg = .true.

    if(nrank == 0) call Print_debug_start_msg("Initializing domain geometric ...")
    ! Build up domain info

    dm%is_periodic(:) = .false.
    if(dm%ibcx(1, 1) == IBC_PERIODIC) dm%is_periodic(1) = .true.
    if(dm%ibcy(1, 1) == IBC_PERIODIC) dm%is_periodic(2) = .true.
    if(dm%ibcz(1, 1) == IBC_PERIODIC) dm%is_periodic(3) = .true.

    dm%np_geo(1) = dm%nc(1) + 1 
    dm%np_geo(2) = dm%nc(2) + 1
    dm%np_geo(3) = dm%nc(3) + 1

    do i = 1, 3
      if ( dm%is_periodic(i) ) then
        dm%np(i) = dm%nc(i)
      else 
        dm%np(i) = dm%np_geo(i)
      end if
    end do

    dm%is_stretching(:) = .false.
    if (dm%istret /= ISTRET_NO) dm%is_stretching(2) = .true.
    
    if(dm%is_stretching(2)) then
      dm%h(2) = ONE / real(dm%nc(2), WP)
    else 
      dm%h(2) = (dm%lyt - dm%lyb) / real(dm%nc(2), WP) ! mean dy
    end if
    dm%h(1) = dm%lxx / real(dm%nc(1), WP)
    dm%h(3) = dm%lzz / real(dm%nc(3), WP)
    dm%h2r(:) = ONE / dm%h(:) / dm%h(:)
    dm%h1r(:) = ONE / dm%h(:)

    !build up index sequence for boundary part

    call Buildup_npneibour_index (dm%np(1), dm%ibcx(1:2,1), dm%ipnbr_p2p(:, :), 'p2p')
    call Buildup_npneibour_index (dm%np(2), dm%ibcy(1:2,1), dm%jpnbr_p2p(:, :), 'p2p' )
    call Buildup_npneibour_index (dm%np(3), dm%ibcz(1:2,1), dm%kpnbr_p2p(:, :), 'p2p' )

    call Buildup_npneibour_index (dm%nc(1), dm%ibcx(1:2,1), dm%ipnbr_p2c(:, :), 'p2c' )
    call Buildup_npneibour_index (dm%nc(2), dm%ibcy(1:2,1), dm%jpnbr_p2c(:, :), 'p2c' )
    call Buildup_npneibour_index (dm%nc(3), dm%ibcz(1:2,1), dm%kpnbr_p2c(:, :), 'p2c' )

    call Buildup_ncneibour_index (dm%nc(1), dm%ibcx(1:2,1), dm%icnbr_c2c(:, :), 'c2c' )
    call Buildup_ncneibour_index (dm%nc(2), dm%ibcy(1:2,1), dm%jcnbr_c2c(:, :), 'c2c' )
    call Buildup_ncneibour_index (dm%nc(3), dm%ibcz(1:2,1), dm%kcnbr_c2c(:, :), 'c2c' )

    call Buildup_ncneibour_index (dm%nc(1), dm%ibcx(1:2,1), dm%icnbr_c2p(:, :), 'c2p' )
    call Buildup_ncneibour_index (dm%nc(2), dm%ibcy(1:2,1), dm%jcnbr_c2p(:, :), 'c2p' )
    call Buildup_ncneibour_index (dm%nc(3), dm%ibcz(1:2,1), dm%kcnbr_c2p(:, :), 'c2p' )


    ! allocate  variables for mapping physical domain to computational domain
    allocate ( dm%yp( dm%np_geo(2) ) ); dm%yp(:) = ZERO
    allocate ( dm%yc( dm%nc(2) ) ); dm%yc(:) = ZERO

    allocate ( dm%yMappingpt( dm%np_geo(2), 3 ) ); dm%yMappingpt(:, :) = ONE
    allocate ( dm%yMappingcc( dm%nc(2),     3 ) ); dm%yMappingcc(:, :) = ONE

    call Buildup_grid_mapping_1D ('nd', dm%np_geo(2), dm%yp(:), dm, dm%yMappingPt(:, :))
    call Buildup_grid_mapping_1D ('cl', dm%nc(2),     dm%yc(:), dm, dm%yMappingcc(:, :))

    ! print out for debugging
    if(dbg) then
      do i = 1, dm%np_geo(2)
        write (OUTPUT_UNIT, '(I5, 1F8.4)') i, dm%yp(i)
      end do
      
      write (OUTPUT_UNIT, '(A)') 'For Point, p2p'
      do i = 1, 4
        j = i
        if(i == 3) j = dm%np(1) - 1
        if(i == 4) j = dm%np(1)
        write (OUTPUT_UNIT, '(A, I4.1, A)') 'For ip =', j, ' its neighbours at given bc'
        write (OUTPUT_UNIT, '(5I7.1)') dm%ipnbr_p2p(i, 1), dm%ipnbr_p2p(i, 2), j, dm%ipnbr_p2p(i, 3), dm%ipnbr_p2p(i, 4)
      end do

      write (OUTPUT_UNIT, '(A)') 'For Point, p2c'
      do i = 1, 4
        j = i
        if(i == 3) j = dm%np(1) - 1
        if(i == 4) j = dm%np(1)
        write (OUTPUT_UNIT, '(A, I4.1, A)') 'For ip =', j, ' its neighbours at given bc'
        write (OUTPUT_UNIT, '(5I7.1)') dm%ipnbr_p2c(i, 1), dm%ipnbr_p2c(i, 2), j, dm%ipnbr_p2c(i, 3), dm%ipnbr_p2c(i, 4)
      end do

      write (OUTPUT_UNIT, '(A)') 'For CC, c2c'
      do i = 1, 4
        j = i
        if(i == 3) j = dm%nc(1) - 1
        if(i == 4) j = dm%nc(1)
        write (OUTPUT_UNIT, '(A, I4.1, A)') 'For ic =', j, ' its neighbours at given bc'
        write (OUTPUT_UNIT, '(5I7.1)') dm%icnbr_c2c(i, 1), dm%icnbr_c2c(i, 2), j, dm%icnbr_c2c(i, 3), dm%icnbr_c2c(i, 4)
      end do

      write (OUTPUT_UNIT, '(A)') 'For CC, c2p'
      do i = 1, 4
        j = i
        if(i == 3) j = dm%nc(1) - 1
        if(i == 4) j = dm%nc(1)
        write (OUTPUT_UNIT, '(A, I4.1, A)') 'For ic =', j, ' its neighbours at given bc'
        write (OUTPUT_UNIT, '(5I7.1)') dm%icnbr_c2p(i, 1), dm%icnbr_c2p(i, 2), j, dm%icnbr_c2p(i, 3), dm%icnbr_c2p(i, 4)
      end do
    end if
    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine  Buildup_geometry_mesh_info
end module geometry_mod

