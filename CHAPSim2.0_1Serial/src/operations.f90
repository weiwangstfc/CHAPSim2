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
!> \file operations.f90
!>
!> \brief A general operation of derivative and interpolation in 1D.
!>
!===============================================================================
module operations
  use precision_mod, only : WP
  implicit none
  private
!----------------------------------------------------------------
! coefficients for TDMA of 1st deriviative  
! to store coefficients for TDMA
! eg, d1fC2C(5, 3, 4)
!     First column: 1:2 for one side b.c.
!                   4:5 for the other side b.c.
!                   3   for interior
!     Second column: 1 for coefficients of f^(1)_{i-1}
!                    2 for coefficients of f^(1)_{i}
!                    3 for coefficients of f^(1)_{i+1}
!     Third column:  for b.c. flags
!     Fourth Column (interpolation only): 1 for orthognal like u in y
!                                         2 for parallel like v in y
!----------------------------------------------------------------
!-------------------------------------------------------------------------------
! for 1st derivative
!-------------------------------------------------------------------------------
  ! collocated C2C
  real(WP) :: d1fC2C(5, 3, 4)
  real(WP) :: d1rC2C(5, 3, 4)
  
  ! collocated P2P
  real(WP) :: d1fP2P(5, 3, 4)
  real(WP) :: d1rP2P(5, 3, 4)

  ! staggered C2P
  real(WP) :: d1fC2P(5, 3, 4)
  real(WP) :: d1rC2P(5, 3, 4)

  ! staggered P2C
  real(WP) :: d1fP2C(5, 3, 4)
  real(WP) :: d1rP2C(5, 3, 4)
!-------------------------------------------------------------------------------
! for 2nd derivative
!-------------------------------------------------------------------------------
  ! collocated C2C
  real(WP), public :: d2fC2C(5, 3, 4)
  real(WP), public :: d2rC2C(5, 4, 4)
  
  ! collocated P2P
  real(WP) :: d2fP2P(5, 3, 4)
  real(WP) :: d2rP2P(5, 4, 4)
!-------------------------------------------------------------------------------
! for iterpolation
!-------------------------------------------------------------------------------
  ! interpolation P2C
  real(WP) :: m1fP2C(5, 3, 4)
  real(WP) :: m1rP2C(5, 3, 4)

  ! interpolation C2P
  real(WP) :: m1fC2P(5, 3, 4)
  real(WP) :: m1rC2P(5, 3, 4)

!-------------------------------------------------------------------------------
! pre-processed TDMA LHS Matrix for 1st deriviative
!-------------------------------------------------------------------------------
  ! x
  real(WP), allocatable :: ad1x_P2P(:)
  real(WP), allocatable :: bd1x_P2P(:)
  real(WP), allocatable :: cd1x_P2P(:)
  real(WP), allocatable :: dd1x_P2P(:)

  real(WP), allocatable :: ad1x_C2C(:)
  real(WP), allocatable :: bd1x_C2C(:)
  real(WP), allocatable :: cd1x_C2C(:)
  real(WP), allocatable :: dd1x_C2C(:)

  real(WP), allocatable :: ad1x_P2C(:)
  real(WP), allocatable :: bd1x_P2C(:)
  real(WP), allocatable :: cd1x_P2C(:)
  real(WP), allocatable :: dd1x_P2C(:)

  real(WP), allocatable :: ad1x_C2P(:)
  real(WP), allocatable :: bd1x_C2P(:)
  real(WP), allocatable :: cd1x_C2P(:)
  real(WP), allocatable :: dd1x_C2P(:)
  ! y
  real(WP), allocatable :: ad1y_P2P(:)
  real(WP), allocatable :: bd1y_P2P(:)
  real(WP), allocatable :: cd1y_P2P(:)
  real(WP), allocatable :: dd1y_P2P(:)

  real(WP), allocatable :: ad1y_C2C(:)
  real(WP), allocatable :: bd1y_C2C(:)
  real(WP), allocatable :: cd1y_C2C(:)
  real(WP), allocatable :: dd1y_C2C(:)

  real(WP), allocatable :: ad1y_P2C(:)
  real(WP), allocatable :: bd1y_P2C(:)
  real(WP), allocatable :: cd1y_P2C(:)
  real(WP), allocatable :: dd1y_P2C(:)

  real(WP), allocatable :: ad1y_C2P(:)
  real(WP), allocatable :: bd1y_C2P(:)
  real(WP), allocatable :: cd1y_C2P(:)
  real(WP), allocatable :: dd1y_C2P(:)

  ! z
  real(WP), allocatable :: ad1z_P2P(:)
  real(WP), allocatable :: bd1z_P2P(:)
  real(WP), allocatable :: cd1z_P2P(:)
  real(WP), allocatable :: dd1z_P2P(:)

  real(WP), allocatable :: ad1z_C2C(:)
  real(WP), allocatable :: bd1z_C2C(:)
  real(WP), allocatable :: cd1z_C2C(:)
  real(WP), allocatable :: dd1z_C2C(:)

  real(WP), allocatable :: ad1z_P2C(:)
  real(WP), allocatable :: bd1z_P2C(:)
  real(WP), allocatable :: cd1z_P2C(:)
  real(WP), allocatable :: dd1z_P2C(:)

  real(WP), allocatable :: ad1z_C2P(:)
  real(WP), allocatable :: bd1z_C2P(:)
  real(WP), allocatable :: cd1z_C2P(:)
  real(WP), allocatable :: dd1z_C2P(:)


!-------------------------------------------------------------------------------
! pre-processed TDMA LHS Matrix for 2nd deriviative
!-------------------------------------------------------------------------------
! x
  real(WP), allocatable :: ad2x_P2P(:)
  real(WP), allocatable :: bd2x_P2P(:)
  real(WP), allocatable :: cd2x_P2P(:)
  real(WP), allocatable :: dd2x_P2P(:)

  real(WP), allocatable :: ad2x_C2C(:)
  real(WP), allocatable :: bd2x_C2C(:)
  real(WP), allocatable :: cd2x_C2C(:)
  real(WP), allocatable :: dd2x_C2C(:)
! y
  real(WP), allocatable :: ad2y_P2P(:)
  real(WP), allocatable :: bd2y_P2P(:)
  real(WP), allocatable :: cd2y_P2P(:)
  real(WP), allocatable :: dd2y_P2P(:)

  real(WP), allocatable :: ad2y_C2C(:)
  real(WP), allocatable :: bd2y_C2C(:)
  real(WP), allocatable :: cd2y_C2C(:)
  real(WP), allocatable :: dd2y_C2C(:)
! z
  real(WP), allocatable :: ad2z_P2P(:)
  real(WP), allocatable :: bd2z_P2P(:)
  real(WP), allocatable :: cd2z_P2P(:)
  real(WP), allocatable :: dd2z_P2P(:)

  real(WP), allocatable :: ad2z_C2C(:)
  real(WP), allocatable :: bd2z_C2C(:)
  real(WP), allocatable :: cd2z_C2C(:)
  real(WP), allocatable :: dd2z_C2C(:)

!-------------------------------------------------------------------------------
! pre-processed TDMA LHS Matrix for mid-point interpolation
!-------------------------------------------------------------------------------
  real(WP), allocatable :: am1x_P2C(:)
  real(WP), allocatable :: bm1x_P2C(:)
  real(WP), allocatable :: cm1x_P2C(:)
  real(WP), allocatable :: dm1x_P2C(:)

  real(WP), allocatable :: am1x_C2P(:)
  real(WP), allocatable :: bm1x_C2P(:)
  real(WP), allocatable :: cm1x_C2P(:)
  real(WP), allocatable :: dm1x_C2P(:)
  
  real(WP), allocatable :: am1y_P2C(:)
  real(WP), allocatable :: bm1y_P2C(:)
  real(WP), allocatable :: cm1y_P2C(:)
  real(WP), allocatable :: dm1y_P2C(:)

  real(WP), allocatable :: am1y_C2P(:)
  real(WP), allocatable :: bm1y_C2P(:)
  real(WP), allocatable :: cm1y_C2P(:)
  real(WP), allocatable :: dm1y_C2P(:)
  
  real(WP), allocatable :: am1z_P2C(:)
  real(WP), allocatable :: bm1z_P2C(:)
  real(WP), allocatable :: cm1z_P2C(:)
  real(WP), allocatable :: dm1z_P2C(:)

  real(WP), allocatable :: am1z_C2P(:)
  real(WP), allocatable :: bm1z_C2P(:)
  real(WP), allocatable :: cm1z_C2P(:)
  real(WP), allocatable :: dm1z_C2P(:)

!-------------------------------------------------------------------------------
! processures
!-------------------------------------------------------------------------------
  private :: Assign_TDMA_coeffs
  private :: Buildup_TDMA_LHS_array
  private :: Prepare_TDMA_interp_RHS_array
  private :: Prepare_TDMA_1deri_RHS_array
  private :: Prepare_TDMA_2deri_RHS_array


  public  :: Prepare_coeffs_for_operations
  public  :: Get_midp_interpolation_1D
  public  :: Get_1st_derivative_1D
  public  :: Get_2nd_derivative_1D

  public  :: Get_x_midp_C2P_3dArray
  public  :: Get_x_midp_P2C_3dArray

  public  :: Get_y_midp_C2P_3dArray
  public  :: Get_y_midp_P2C_3dArray

  public  :: Get_z_midp_C2P_3dArray
  public  :: Get_z_midp_P2C_3dArray

  public  :: Get_x_1st_derivative_C2C_3dArray
  public  :: Get_x_1st_derivative_P2P_3dArray
  public  :: Get_x_1st_derivative_C2P_3dArray
  public  :: Get_x_1st_derivative_P2C_3dArray

  public  :: Get_y_1st_derivative_C2C_3dArray
  public  :: Get_y_1st_derivative_P2P_3dArray
  public  :: Get_y_1st_derivative_C2P_3dArray
  public  :: Get_y_1st_derivative_P2C_3dArray

  public  :: Get_z_1st_derivative_C2C_3dArray
  public  :: Get_z_1st_derivative_P2P_3dArray
  public  :: Get_z_1st_derivative_C2P_3dArray
  public  :: Get_z_1st_derivative_P2C_3dArray

  public  :: Get_x_2nd_derivative_C2C_3dArray
  public  :: Get_x_2nd_derivative_P2P_3dArray

  public  :: Get_y_2nd_derivative_C2C_3dArray
  public  :: Get_y_2nd_derivative_P2P_3dArray

  public  :: Get_z_2nd_derivative_C2C_3dArray
  public  :: Get_z_2nd_derivative_P2P_3dArray

  public  :: Get_volumetric_average_3d

contains
!===============================================================================
!===============================================================================
!> \brief Assigned the cooefficients for the compact schemes     
!>
!> This subroutine is called once locally.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iaccu         the accuracy given by user
!_______________________________________________________________________________
  subroutine Assign_TDMA_coeffs(iaccu)
    use parameters_constant_mod
    use input_general_mod
    implicit none

    integer(4), intent(in) :: iaccu

    real(WP) :: alpha,  a,  b,  c,  d
    real(WP) :: alpha1, a1, b1, c1, d1
    real(WP) :: alpha2, a2, b2, c2, d2

    call Print_debug_start_msg &
         ("Assigning coefficient matrix for the compact FD ...")
!_____________________________________________________________________________!
!1st derivative on collocated grids
!_______________________________________________________________________________!
    ! C2C/P2P coefficients
    if (iaccu == IACCU_CD2) then
      alpha = ZERO
      a = ONE
      b = ZERO
      c = ZERO ! not used
    else if (iaccu == IACCU_CD4) then
      alpha = ZERO
      a = FOUR / THREE
      b = -ONE / THREE
      c = ZERO ! not used
    else if (iaccu == IACCU_CP4) then
      alpha = ONE / FOUR
      a = THREE / TWO
      b = ZERO
      c = ZERO ! not used
    else if (iaccu == IACCU_CP6) then
      alpha = ONE / THREE
      a = FOURTEEN / NINE
      b = ONE / NINE
      c = ZERO ! not used
    else ! default 2nd CD
      alpha = ZERO
      a = ONE
      b = ZERO
      c = ZERO ! not used
    end if

    !C2C for periodic b.c.
    d1fC2C(1:5, 1, IBC_PERIODIC) = alpha
    d1fC2C(1:5, 2, IBC_PERIODIC) = ONE
    d1fC2C(1:5, 3, IBC_PERIODIC) = alpha
    d1rC2C(1:5, 1, IBC_PERIODIC) = a / TWO ! a/2
    d1rC2C(1:5, 2, IBC_PERIODIC) = b / FOUR ! b/4
    d1rC2C(1:5, 3, IBC_PERIODIC) = c ! not used

    !P2P for periodic b.c.
    d1fP2P(:, :, IBC_PERIODIC) = d1fC2C(:, :, IBC_PERIODIC)
    d1rP2P(:, :, IBC_PERIODIC) = d1rC2C(:, :, IBC_PERIODIC)

    !C2C for symmetric b.c.
    d1fC2C(1, 1, IBC_SYMMETRIC) = ZERO ! not used
    d1fC2C(1, 2, IBC_SYMMETRIC) = ONE - alpha
    d1fC2C(1, 3, IBC_SYMMETRIC) = alpha

    d1fC2C(2:4, 1, IBC_SYMMETRIC) = alpha
    d1fC2C(2:4, 2, IBC_SYMMETRIC) = ONE
    d1fC2C(2:4, 3, IBC_SYMMETRIC) = alpha

    d1fC2C(5, 1, IBC_SYMMETRIC) = alpha
    d1fC2C(5, 2, IBC_SYMMETRIC) = ONE - alpha
    d1fC2C(5, 3, IBC_SYMMETRIC) = ZERO ! not used

    d1rC2C(1:5, 1, IBC_SYMMETRIC) = a / TWO ! a/2
    d1rC2C(1:5, 2, IBC_SYMMETRIC) = b / FOUR ! b/4
    d1rC2C(1:5, 3, IBC_SYMMETRIC) = c ! not used

    !P2P for symmetric b.c.
    d1fP2P(1, 1, IBC_SYMMETRIC) = ZERO ! not used
    d1fP2P(1, 2, IBC_SYMMETRIC) = ONE
    d1fP2P(1, 3, IBC_SYMMETRIC) = ZERO

    d1fP2P(2:4, :, IBC_SYMMETRIC) = d1fC2C(2:4, :, IBC_SYMMETRIC)

    d1fP2P(5, 1, IBC_SYMMETRIC) = ZERO
    d1fP2P(5, 2, IBC_SYMMETRIC) = ONE
    d1fP2P(5, 3, IBC_SYMMETRIC) = ZERO ! not used

    d1rP2P(:, :, IBC_SYMMETRIC) = d1rC2C(:, :, IBC_SYMMETRIC)

    !C2C for asymmetric b.c.
    d1fC2C(:, :, IBC_ASYMMETRIC) = d1fC2C(:, :, IBC_SYMMETRIC)
    d1rC2C(:, :, IBC_ASYMMETRIC) = d1rC2C(:, :, IBC_SYMMETRIC)
    !P2P for asymmetric b.c.
    d1fP2P(:, :, IBC_ASYMMETRIC) = d1fP2P(:, :, IBC_SYMMETRIC)
    d1rP2P(:, :, IBC_ASYMMETRIC) = d1rP2P(:, :, IBC_SYMMETRIC)

    !C2C/P2P for Dirichlet B.C.
    if (iaccu == IACCU_CD2) then

      alpha1 = ZERO
      a1 = -THREE / TWO
      b1 = TWO
      c1 = -ONE / TWO

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CD4) then ! degrade to 2nd CD

      alpha1 = ZERO
      a1 = -THREE / TWO
      b1 = TWO
      c1 = -ONE / TWO

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CP4) then ! degrade to 3rd CP

      alpha1 = TWO
      a1 = -FIVE / TWO
      b1 = TWO
      c1 = ONE / TWO

      alpha2 = ONE / FOUR
      a2 = THREE / TWO
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CP6) then ! degrade to 3rd CP

      alpha1 = TWO
      a1 = -FIVE / TWO
      b1 = TWO
      c1 = ONE / TWO

      alpha2 = ONE / FOUR
      a2 = THREE / TWO
      b2 = ZERO ! not used
      c2 = ZERO ! not used
      
    else ! default 2nd CD
      alpha1 = ZERO
      a1 = -THREE / TWO
      b1 = TWO
      c1 = -ONE / TWO

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used
      
    end if

    d1fC2C(1, 1, IBC_UDIRICHLET) = alpha1 ! not used
    d1fC2C(1, 2, IBC_UDIRICHLET) = ONE
    d1fC2C(1, 3, IBC_UDIRICHLET) = alpha1
    d1rC2C(1, 1, IBC_UDIRICHLET) = a1
    d1rC2C(1, 2, IBC_UDIRICHLET) = b1
    d1rC2C(1, 3, IBC_UDIRICHLET) = c1

    d1fC2C(2, 1, IBC_UDIRICHLET) = alpha2
    d1fC2C(2, 2, IBC_UDIRICHLET) = ONE
    d1fC2C(2, 3, IBC_UDIRICHLET) = alpha2
    d1rC2C(2, 1, IBC_UDIRICHLET) = a2 / TWO
    d1rC2C(2, 2, IBC_UDIRICHLET) = b2 / FOUR ! not used
    d1rC2C(2, 3, IBC_UDIRICHLET) = c2 ! not used

    d1fC2C(3, 1, IBC_UDIRICHLET) = alpha
    d1fC2C(3, 2, IBC_UDIRICHLET) = ONE
    d1fC2C(3, 3, IBC_UDIRICHLET) = alpha
    d1rC2C(3, 1, IBC_UDIRICHLET) = a / TWO ! a/2
    d1rC2C(3, 2, IBC_UDIRICHLET) = b / FOUR ! b/4
    d1rC2C(3, 3, IBC_UDIRICHLET) = c ! not used

    d1fC2C(4, 1, IBC_UDIRICHLET) = alpha2
    d1fC2C(4, 2, IBC_UDIRICHLET) = ONE
    d1fC2C(4, 3, IBC_UDIRICHLET) = alpha2
    d1rC2C(4, 1, IBC_UDIRICHLET) = a2 / TWO
    d1rC2C(4, 2, IBC_UDIRICHLET) = b2 / FOUR ! not used
    d1rC2C(4, 3, IBC_UDIRICHLET) = c2 ! not used

    d1fC2C(5, 1, IBC_UDIRICHLET) = alpha1
    d1fC2C(5, 2, IBC_UDIRICHLET) = ONE
    d1fC2C(5, 3, IBC_UDIRICHLET) = alpha1 ! not used
    d1rC2C(5, 1, IBC_UDIRICHLET) = -a1
    d1rC2C(5, 2, IBC_UDIRICHLET) = -b1
    d1rC2C(5, 3, IBC_UDIRICHLET) = -c1

    d1fP2P(:, :, IBC_UDIRICHLET) = d1fC2C(:, :, IBC_UDIRICHLET)
    d1rP2P(:, :, IBC_UDIRICHLET) = d1rC2C(:, :, IBC_UDIRICHLET)
!______________________________________________________________________________!
!1st derivative on staggered grids P2C and C2P
!______________________________________________________________________________!
    if (iaccu == IACCU_CD2) then
      alpha = ZERO
      a = ONE
      b = ZERO
      c = ZERO ! not used
    else if (iaccu == IACCU_CD4) then
      alpha = ZERO
      a = NINE / EIGHT
      b = -ONE / EIGHT
      c = ZERO ! not used
    else if (iaccu == IACCU_CP4) then
      alpha = ONE / TWENTYTWO
      a = TWELVE / ELEVEN
      b = ZERO
      c = ZERO ! not used
    else if (iaccu == IACCU_CP6) then
      alpha = NINE / SIXTYTWO
      a = SIXTYTHREE / SIXTYTWO
      b = SEVENTEEN / SIXTYTWO
      c = ZERO ! not used
    else  ! default 2nd CD
      alpha = ZERO
      a = ONE
      b = ZERO
      c = ZERO ! not used
      
    end if

    !C2P for periodic b.c.
    d1fC2P(1:5, 1, IBC_PERIODIC) = alpha
    d1fC2P(1:5, 2, IBC_PERIODIC) = ONE
    d1fC2P(1:5, 3, IBC_PERIODIC) = alpha
    d1rC2P(1:5, 1, IBC_PERIODIC) = a ! a
    d1rC2P(1:5, 2, IBC_PERIODIC) = b / THREE ! b/3
    d1rC2P(1:5, 3, IBC_PERIODIC) = c ! not used

    !P2C for periodic b.c.
    d1fP2C(:, :, IBC_PERIODIC) = d1fC2P(:, :, IBC_PERIODIC)
    d1rP2C(:, :, IBC_PERIODIC) = d1rC2P(:, :, IBC_PERIODIC)

    !C2P for symmetric 
    d1fC2P(:, :, IBC_SYMMETRIC) = d1fP2P(:, :, IBC_SYMMETRIC)
    d1rC2P(:, :, IBC_SYMMETRIC) = d1rC2P(:, :, IBC_PERIODIC)

    !P2C for symmetric 
    d1fP2C(:, :, IBC_SYMMETRIC) = d1fC2C(:, :, IBC_SYMMETRIC)
    d1rP2C(:, :, IBC_SYMMETRIC) = d1rC2P(:, :, IBC_SYMMETRIC)

    !C2P for asymmetric 
    d1fC2P(:, :, IBC_ASYMMETRIC) = d1fC2P(:, :, IBC_SYMMETRIC)
    d1rC2P(:, :, IBC_ASYMMETRIC) = d1rC2P(:, :, IBC_SYMMETRIC)

    !P2C for asymmetric 
    d1fP2C(:, :, IBC_ASYMMETRIC) = d1fP2C(:, :, IBC_SYMMETRIC)
    d1rP2C(:, :, IBC_ASYMMETRIC) = d1rP2C(:, :, IBC_SYMMETRIC)

    !P2C for Dirichlet B.C.
    if (iaccu == IACCU_CD2) then
      alpha1 = ZERO
      a1 = -ONE
      b1 = ONE
      c1 = ZERO

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CD4) then ! degrade to 2nd CD

      alpha1 = ZERO
      a1 = -ONE
      b1 = ONE
      c1 = ZERO

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CP4) then ! degrade to 3rd CP

      alpha1 = -ONE
      a1 = -ONE
      b1 = TWO
      c1 = -ONE

      alpha2 = ONE / TWENTYTWO
      a2 = TWELVE / ELEVEN
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CP6) then ! degrade to 3rd CP

      alpha1 = -ONE
      a1 = -ONE
      b1 = TWO
      c1 = -ONE

      alpha2 = ONE / TWENTYTWO
      a2 = TWELVE / ELEVEN
      b2 = ZERO ! not used
      c2 = ZERO ! not used
      
    else  ! default 2nd CD
      alpha1 = ZERO
      a1 = -ONE
      b1 = ONE
      c1 = ZERO

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used
    end if

    d1fP2C(1, 1, IBC_UDIRICHLET) = alpha1 ! not used
    d1fP2C(1, 2, IBC_UDIRICHLET) = ONE
    d1fP2C(1, 3, IBC_UDIRICHLET) = alpha1
    d1rP2C(1, 1, IBC_UDIRICHLET) = a1
    d1rP2C(1, 2, IBC_UDIRICHLET) = b1
    d1rP2C(1, 3, IBC_UDIRICHLET) = c1

    d1fP2C(2, 1, IBC_UDIRICHLET) = alpha2
    d1fP2C(2, 2, IBC_UDIRICHLET) = ONE
    d1fP2C(2, 3, IBC_UDIRICHLET) = alpha2
    d1rP2C(2, 1, IBC_UDIRICHLET) = a2
    d1rP2C(2, 2, IBC_UDIRICHLET) = b2 / THREE ! not used
    d1rP2C(2, 3, IBC_UDIRICHLET) = c2 ! not used

    d1fP2C(3, 1, IBC_UDIRICHLET) = alpha
    d1fP2C(3, 2, IBC_UDIRICHLET) = ONE
    d1fP2C(3, 3, IBC_UDIRICHLET) = alpha
    d1rP2C(3, 1, IBC_UDIRICHLET) = a
    d1rP2C(3, 2, IBC_UDIRICHLET) = b / THREE
    d1rP2C(3, 3, IBC_UDIRICHLET) = c ! not used

    d1fP2C(4, 1, IBC_UDIRICHLET) = alpha2
    d1fP2C(4, 2, IBC_UDIRICHLET) = ONE
    d1fP2C(4, 3, IBC_UDIRICHLET) = alpha2
    d1rP2C(4, 1, IBC_UDIRICHLET) = a2
    d1rP2C(4, 2, IBC_UDIRICHLET) = b2 / THREE ! not used
    d1rP2C(4, 3, IBC_UDIRICHLET) = c2 ! not used

    d1fP2C(5, 1, IBC_UDIRICHLET) = alpha1
    d1fP2C(5, 2, IBC_UDIRICHLET) = ONE
    d1fP2C(5, 3, IBC_UDIRICHLET) = alpha1 ! not used
    d1rP2C(5, 1, IBC_UDIRICHLET) = -a1
    d1rP2C(5, 2, IBC_UDIRICHLET) = -b1
    d1rP2C(5, 3, IBC_UDIRICHLET) = -c1

    !C2P for Dirichlet B.C.
    if (iaccu == IACCU_CD2) then

      alpha1 = ZERO
      a1 = -TWO
      b1 = THREE
      c1 = -ONE

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CD4) then ! degrade to 2nd CD

      alpha1 = ZERO
      a1 = -TWO
      b1 = THREE
      c1 = -ONE

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CP4) then ! degrade to 3rd CP

      alpha1 = TWENTYTHREE
      a1 = -TWENTYFIVE
      b1 = TWENTYSIX
      c1 = -ONE

      alpha2 = ONE / TWENTYTWO
      a2 = TWELVE / ELEVEN
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CP6) then ! degrade to 3rd CP

      alpha1 = TWENTYTHREE
      a1 = -TWENTYFIVE
      b1 = TWENTYSIX
      c1 = -ONE

      alpha2 = ONE / TWENTYTWO
      a2 = TWELVE / ELEVEN
      b2 = ZERO ! not used
      c2 = ZERO ! not used
      
    else  ! default 2nd CD
      alpha1 = ZERO
      a1 = -TWO
      b1 = THREE
      c1 = -ONE

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    end if

    d1fC2P(1, 1, IBC_UDIRICHLET) = alpha1 ! not used
    d1fC2P(1, 2, IBC_UDIRICHLET) = ONE
    d1fC2P(1, 3, IBC_UDIRICHLET) = alpha1
    d1rC2P(1, 1, IBC_UDIRICHLET) = a1
    d1rC2P(1, 2, IBC_UDIRICHLET) = b1
    d1rC2P(1, 3, IBC_UDIRICHLET) = c1

    d1fC2P(2, 1, IBC_UDIRICHLET) = alpha2
    d1fC2P(2, 2, IBC_UDIRICHLET) = ONE
    d1fC2P(2, 3, IBC_UDIRICHLET) = alpha2
    d1rC2P(2, 1, IBC_UDIRICHLET) = a2
    d1rC2P(2, 2, IBC_UDIRICHLET) = b2 / THREE ! not used
    d1rC2P(2, 3, IBC_UDIRICHLET) = c2 ! not used

    d1fC2P(3, 1, IBC_UDIRICHLET) = alpha
    d1fC2P(3, 2, IBC_UDIRICHLET) = ONE
    d1fC2P(3, 3, IBC_UDIRICHLET) = alpha
    d1rC2P(3, 1, IBC_UDIRICHLET) = a
    d1rC2P(3, 2, IBC_UDIRICHLET) = b / THREE
    d1rC2P(3, 3, IBC_UDIRICHLET) = c ! not used

    d1fC2P(4, 1, IBC_UDIRICHLET) = alpha2
    d1fC2P(4, 2, IBC_UDIRICHLET) = ONE
    d1fC2P(4, 3, IBC_UDIRICHLET) = alpha2
    d1rC2P(4, 1, IBC_UDIRICHLET) = a2
    d1rC2P(4, 2, IBC_UDIRICHLET) = b2 / THREE ! not used
    d1rC2P(4, 3, IBC_UDIRICHLET) = c2 ! not used

    d1fC2P(5, 1, IBC_UDIRICHLET) = alpha1
    d1fC2P(5, 2, IBC_UDIRICHLET) = ONE
    d1fC2P(5, 3, IBC_UDIRICHLET) = alpha1 ! not used
    d1rC2P(5, 1, IBC_UDIRICHLET) = -a1
    d1rC2P(5, 2, IBC_UDIRICHLET) = -b1
    d1rC2P(5, 3, IBC_UDIRICHLET) = -c1

!______________________________________________________________________________!
!interpolation. P2C and C2P
!______________________________________________________________________________!
    if (iaccu == IACCU_CD2) then
      alpha = ZERO
      a = ONE
      b = ZERO
      c = ZERO ! not used
    else if (iaccu == IACCU_CD4) then
      alpha = ZERO
      a = NINE / EIGHT
      b = -ONE / EIGHT
      c = ZERO ! not used
    else if (iaccu == IACCU_CP4) then
      alpha = ONE / SIX
      a = FOUR / THREE
      b = ZERO
      c = ZERO ! not used
    else if (iaccu == IACCU_CP6) then
      alpha = THREE / TEN
      a = THREE / TWO
      b = ONE / TEN
      c = ZERO ! not used
    else  ! default 2nd CD
      alpha = ZERO
      a = ONE
      b = ZERO
      c = ZERO ! not used
    end if

    !C2P: i'_max = np
    !     alpha * f_{i'-1} + f_{i'} + f_{i'+1} = b/2 * (f_{i+1} + f_{i-2}) + a/2 * (f_{i} + f_{i-1})
    !P2C: i_max = nc
    !     alpha * f_{i-1} + f_{i} + f_{i+1} = b/2 * (f_{i'+2} + f_{i'-1}) + a/2 * (f_{i'} + f_{i'+1})

    !C2P for periodic b.c.
    m1fC2P(1:5, 1, IBC_PERIODIC) = alpha
    m1fC2P(1:5, 2, IBC_PERIODIC) = ONE
    m1fC2P(1:5, 3, IBC_PERIODIC) = alpha
    m1rC2P(1:5, 1, IBC_PERIODIC) = a / TWO
    m1rC2P(1:5, 2, IBC_PERIODIC) = b / TWO
    m1rC2P(1:5, 3, IBC_PERIODIC) = c ! not used

    !P2C for periodic b.c.
    m1fP2C(:, :, IBC_PERIODIC) = m1fC2P(:, :, IBC_PERIODIC)
    m1rP2C(:, :, IBC_PERIODIC) = m1rC2P(:, :, IBC_PERIODIC)

    !C2P for symmetric, orthogonal, eg. u in y direction.
    m1fC2P(1, 1, IBC_SYMMETRIC) = ZERO ! not used
    m1fC2P(1, 2, IBC_SYMMETRIC) = ONE
    m1fC2P(1, 3, IBC_SYMMETRIC) = alpha + alpha

    m1fC2P(2:4, 1, IBC_SYMMETRIC) = alpha
    m1fC2P(2:4, 2, IBC_SYMMETRIC) = ONE
    m1fC2P(2:4, 3, IBC_SYMMETRIC) = alpha

    m1fC2P(5, 1, IBC_SYMMETRIC) = alpha + alpha
    m1fC2P(5, 2, IBC_SYMMETRIC) = ONE
    m1fC2P(5, 3, IBC_SYMMETRIC) = ZERO ! not used.

    m1rC2P(1:5, 1, IBC_SYMMETRIC) = a / TWO
    m1rC2P(1:5, 2, IBC_SYMMETRIC) = b / TWO
    m1rC2P(1:5, 3, IBC_SYMMETRIC) = ZERO ! not used. 

    !C2P for symmetric, parallel, eg. v in y direction.
    m1fC2P(1, 1, IBC_ASYMMETRIC) = ZERO ! not used
    m1fC2P(1, 2, IBC_ASYMMETRIC) = ONE
    m1fC2P(1, 3, IBC_ASYMMETRIC) = alpha - alpha

    m1fC2P(2:4, 1, IBC_ASYMMETRIC) = alpha
    m1fC2P(2:4, 2, IBC_ASYMMETRIC) = ONE
    m1fC2P(2:4, 3, IBC_ASYMMETRIC) = alpha

    m1fC2P(5, 1, IBC_ASYMMETRIC) = alpha - alpha
    m1fC2P(5, 2, IBC_ASYMMETRIC) = ONE
    m1fC2P(5, 3, IBC_ASYMMETRIC) = ZERO ! not used.

    m1rC2P(1:5, 1, IBC_ASYMMETRIC) = a / TWO
    m1rC2P(1:5, 2, IBC_ASYMMETRIC) = b / TWO
    m1rC2P(1:5, 3, IBC_ASYMMETRIC) = ZERO ! not used. 

    !P2C for symmetric, orthogonal, eg. u in y direction.
    m1fP2C(1, 1, IBC_SYMMETRIC) = ZERO ! not used
    m1fP2C(1, 2, IBC_SYMMETRIC) = ONE + alpha
    m1fP2C(1, 3, IBC_SYMMETRIC) = alpha

    m1fP2C(2:4, 1, IBC_SYMMETRIC) = alpha
    m1fP2C(2:4, 2, IBC_SYMMETRIC) = ONE
    m1fP2C(2:4, 3, IBC_SYMMETRIC) = alpha

    m1fP2C(5, 1, IBC_SYMMETRIC) = alpha
    m1fP2C(5, 2, IBC_SYMMETRIC) = ONE + alpha
    m1fP2C(5, 3, IBC_SYMMETRIC) = ZERO ! not used.

    m1rP2C(1:5, 1, IBC_SYMMETRIC) = a / TWO
    m1rP2C(1:5, 2, IBC_SYMMETRIC) = b / TWO
    m1rP2C(1:5, 3, IBC_SYMMETRIC) = ZERO ! not used. 

    !P2C for symmetric, parallel, eg. v in y direction.
    m1fP2C(1, 1, IBC_ASYMMETRIC) = ZERO ! not used
    m1fP2C(1, 2, IBC_ASYMMETRIC) = ONE - alpha
    m1fP2C(1, 3, IBC_ASYMMETRIC) = alpha

    m1fP2C(2:4, 1, IBC_ASYMMETRIC) = alpha
    m1fP2C(2:4, 2, IBC_ASYMMETRIC) = ONE
    m1fP2C(2:4, 3, IBC_ASYMMETRIC) = alpha

    m1fP2C(5, 1, IBC_ASYMMETRIC) = alpha
    m1fP2C(5, 2, IBC_ASYMMETRIC) = ONE - alpha
    m1fP2C(5, 3, IBC_ASYMMETRIC) = ZERO ! not used.

    m1rP2C(1:5, 1, IBC_ASYMMETRIC) = a / TWO
    m1rP2C(1:5, 2, IBC_ASYMMETRIC) = b / TWO
    m1rP2C(1:5, 3, IBC_ASYMMETRIC) = ZERO ! not used. 

    !P2C for Dirichlet B.C.
    if (iaccu == IACCU_CD2) then
      alpha1 = ZERO
      a1 = THREE / EIGHT
      b1 = THREE / FOUR
      c1 = -ONE / EIGHT

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CD4) then ! degrade to 2nd CD

      alpha1 = ZERO
      a1 = THREE / EIGHT
      b1 = THREE / FOUR
      c1 = -ONE / EIGHT

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CP4) then ! degrade to 3rd CP

      alpha1 = ONE
      a1 = ONE / FOUR
      b1 = THREE / TWO
      c1 = ONE / FOUR

      alpha2 = ONE / SIX
      a2 = FOUR / THREE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CP6) then ! degrade to 3rd CP

      alpha1 = ONE
      a1 = ONE / FOUR
      b1 = THREE / TWO
      c1 = ONE / FOUR

      alpha2 = ONE / SIX
      a2 = FOUR / THREE
      b2 = ZERO ! not used
      c2 = ZERO ! not used
      
    else  ! default 2nd CD
      alpha1 = ZERO
      a1 = THREE / EIGHT
      b1 = THREE / FOUR
      c1 = -ONE / EIGHT

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used
      
    end if
    !P2C
    m1fP2C(1, 1, IBC_UDIRICHLET) = ZERO ! not used
    m1fP2C(1, 2, IBC_UDIRICHLET) = ONE
    m1fP2C(1, 3, IBC_UDIRICHLET) = alpha1
    m1rP2C(1, 1, IBC_UDIRICHLET) = a1
    m1rP2C(1, 2, IBC_UDIRICHLET) = b1
    m1rP2C(1, 3, IBC_UDIRICHLET) = c1

    m1fP2C(2, 1, IBC_UDIRICHLET) = alpha2
    m1fP2C(2, 2, IBC_UDIRICHLET) = ONE
    m1fP2C(2, 3, IBC_UDIRICHLET) = alpha2
    m1rP2C(2, 1, IBC_UDIRICHLET) = a2 / TWO
    m1rP2C(2, 2, IBC_UDIRICHLET) = ZERO ! not used
    m1rP2C(2, 3, IBC_UDIRICHLET) = ZERO ! not used

    m1fP2C(3, 1, IBC_UDIRICHLET) = alpha
    m1fP2C(3, 2, IBC_UDIRICHLET) = ONE
    m1fP2C(3, 3, IBC_UDIRICHLET) = alpha
    m1rP2C(3, 1, IBC_UDIRICHLET) = a / TWO
    m1rP2C(3, 2, IBC_UDIRICHLET) = b / TWO
    m1rP2C(3, 3, IBC_UDIRICHLET) = ZERO ! not used

    m1fP2C(4, 1, IBC_UDIRICHLET) = alpha2
    m1fP2C(4, 2, IBC_UDIRICHLET) = ONE
    m1fP2C(4, 3, IBC_UDIRICHLET) = alpha2
    m1rP2C(4, 1, IBC_UDIRICHLET) = a2 / TWO
    m1rP2C(4, 2, IBC_UDIRICHLET) = ZERO ! not used
    m1rP2C(4, 3, IBC_UDIRICHLET) = ZERO ! not used

    m1fP2C(5, 1, IBC_UDIRICHLET) = alpha1
    m1fP2C(5, 2, IBC_UDIRICHLET) = ONE
    m1fP2C(5, 3, IBC_UDIRICHLET) = ZERO ! not used
    m1rP2C(5, 1, IBC_UDIRICHLET) = a1
    m1rP2C(5, 2, IBC_UDIRICHLET) = b1
    m1rP2C(5, 3, IBC_UDIRICHLET) = c1

    !C2P for Dirichlet B.C.
    if (iaccu == IACCU_CD2) then

      alpha1 = ZERO
      a1 = FIFTEEN / EIGHT
      b1 = - FIVE / FOUR
      c1 = THREE / EIGHT

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CD4) then ! degrade to 2nd CD

      alpha1 = ZERO
      a1 = FIFTEEN / EIGHT
      b1 = - FIVE / FOUR
      c1 = THREE / EIGHT

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CP4) then ! degrade to 3rd CP

      alpha1 = FIVE
      a1 = FIFTEEN / FOUR
      b1 = FIVE / TWO
      c1 = -ONE / FOUR

      alpha2 = ONE / SIX 
      a2 = FOUR / THREE
      b2 = ZERO ! not used
      c2 = ZERO ! not used

    else if (iaccu == IACCU_CP6) then ! degrade to 3rd CP

      alpha1 = FIVE
      a1 = FIFTEEN / FOUR
      b1 = FIVE / TWO
      c1 = -ONE / FOUR

      alpha2 = ONE / SIX
      a2 = FOUR / THREE
      b2 = ZERO ! not used
      c2 = ZERO ! not used
      
    else  ! default 2nd CD

      alpha1 = ZERO
      a1 = FIFTEEN / EIGHT
      b1 = - FIVE / FOUR
      c1 = THREE / EIGHT

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO ! not used
      c2 = ZERO ! not used
     
    end if

    m1fC2P(1, 1, IBC_UDIRICHLET) = ZERO ! not used
    m1fC2P(1, 2, IBC_UDIRICHLET) = ONE
    m1fC2P(1, 3, IBC_UDIRICHLET) = alpha1
    m1rC2P(1, 1, IBC_UDIRICHLET) = a1
    m1rC2P(1, 2, IBC_UDIRICHLET) = b1
    m1rC2P(1, 3, IBC_UDIRICHLET) = c1

    m1fC2P(2, 1, IBC_UDIRICHLET) = alpha2
    m1fC2P(2, 2, IBC_UDIRICHLET) = ONE
    m1fC2P(2, 3, IBC_UDIRICHLET) = alpha2
    m1rC2P(2, 1, IBC_UDIRICHLET) = a2 / TWO
    m1rC2P(2, 2, IBC_UDIRICHLET) = ZERO ! not used
    m1rC2P(2, 3, IBC_UDIRICHLET) = ZERO ! not used

    m1fC2P(3, 1, IBC_UDIRICHLET) = alpha
    m1fC2P(3, 2, IBC_UDIRICHLET) = ONE
    m1fC2P(3, 3, IBC_UDIRICHLET) = alpha
    m1rC2P(3, 1, IBC_UDIRICHLET) = a / TWO
    m1rC2P(3, 2, IBC_UDIRICHLET) = b / TWO
    m1rC2P(3, 3, IBC_UDIRICHLET) = ZERO ! not used

    m1fC2P(4, 1, IBC_UDIRICHLET) = alpha2
    m1fC2P(4, 2, IBC_UDIRICHLET) = ONE
    m1fC2P(4, 3, IBC_UDIRICHLET) = alpha2
    m1rC2P(4, 1, IBC_UDIRICHLET) = a2 / TWO
    m1rC2P(4, 2, IBC_UDIRICHLET) = ZERO ! not used
    m1rC2P(4, 3, IBC_UDIRICHLET) = ZERO ! not used

    m1fC2P(5, 1, IBC_UDIRICHLET) = alpha1
    m1fC2P(5, 2, IBC_UDIRICHLET) = ONE
    m1fC2P(5, 3, IBC_UDIRICHLET) = ZERO ! not used
    m1rC2P(5, 1, IBC_UDIRICHLET) = a1
    m1rC2P(5, 2, IBC_UDIRICHLET) = b1
    m1rC2P(5, 3, IBC_UDIRICHLET) = c1
    
!______________________________________________________________________________!
! 2nd diriviative P2P and C2C
!______________________________________________________________________________!
    if (iaccu == IACCU_CD2) then
      alpha = ZERO
      a = ONE
      b = ZERO
      c = ZERO ! not used
      d = ZERO ! not used
    else if (iaccu == IACCU_CD4) then
      alpha = ZERO
      a = FOUR / THREE
      b = -ONE / THREE
      c = ZERO ! not used
      d = ZERO ! not used
    else if (iaccu == IACCU_CP4) then
      alpha = ONE / TEN
      a = SIX / FIVE
      b = ZERO
      c = ZERO ! not used
      d = ZERO ! not used
    else if (iaccu == IACCU_CP6) then
      alpha = TWO / ELEVEN
      a = TWELVE / ELEVEN
      b = THREE / ELEVEN
      c = ZERO ! not used
      d = ZERO ! not used
    else  ! default 2nd CD
      alpha = ZERO
      a = ONE
      b = ZERO
      c = ZERO ! not used
      d = ZERO ! not used
    end if

    !C2C for periodic b.c.
    d2fC2C(1:5, 1, IBC_PERIODIC) = alpha
    d2fC2C(1:5, 2, IBC_PERIODIC) = ONE
    d2fC2C(1:5, 3, IBC_PERIODIC) = alpha

    d2rC2C(1:5, 1, IBC_PERIODIC) = a / ONE
    d2rC2C(1:5, 2, IBC_PERIODIC) = b / FOUR
    d2rC2C(1:5, 3, IBC_PERIODIC) = c ! not used
    d2rC2C(1:5, 4, IBC_PERIODIC) = d ! not used

    !P2P for periodic b.c.
    d2fP2P(:, :, IBC_PERIODIC) = d2fC2C(:, :, IBC_PERIODIC)
    d2rP2P(:, :, IBC_PERIODIC) = d2rC2C(:, :, IBC_PERIODIC)

    !P2C for Dirichlet B.C.
    if (iaccu == IACCU_CD2) then
      alpha1 = ZERO
      a1 = TWO
      b1 = -FIVE
      c1 = FOUR
      d1 = -ONE

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO
      c2 = ZERO ! not used
      d2 = ZERO ! not used

    else if (iaccu == IACCU_CD4) then ! degrade to 2nd CD

      alpha1 = ZERO
      a1 = TWO
      b1 = -FIVE
      c1 = FOUR
      d1 = -ONE

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO
      c2 = ZERO ! not used
      d2 = ZERO ! not used

    else if (iaccu == IACCU_CP4) then ! degrade to 3rd CP

      alpha1 = ELEVEN
      a1 = THIRTEEN
      b1 = -TWENTYSEVEN
      c1 = FIFTEEN
      d1 = -ONE

      alpha2 = ONE / TEN
      a2 = SIX / FIVE
      b2 = ZERO
      c2 = ZERO ! not used
      d2 = ZERO ! not used

    else if (iaccu == IACCU_CP6) then ! degrade to 3rd CP

      alpha1 = ELEVEN
      a1 = THIRTEEN
      b1 = -TWENTYSEVEN
      c1 = FIFTEEN
      d1 = -ONE

      alpha2 = ONE / TEN
      a2 = SIX / FIVE
      b2 = ZERO
      c2 = ZERO ! not used
      d2 = ZERO ! not used
      
    else  ! default 2nd CD
      alpha1 = ZERO
      a1 = TWO
      b1 = -FIVE
      c1 = FOUR
      d1 = -ONE

      alpha2 = ZERO
      a2 = ONE
      b2 = ZERO
      c2 = ZERO ! not used
      d2 = ZERO ! not used
      
    end if

    !C2C i = 1
    d2fC2C(1, 1, IBC_UDIRICHLET) = ZERO ! not used
    d2fC2C(1, 2, IBC_UDIRICHLET) = ONE
    d2fC2C(1, 3, IBC_UDIRICHLET) = alpha1

    d2rC2C(1, 1, IBC_UDIRICHLET) = a1
    d2rC2C(1, 2, IBC_UDIRICHLET) = b1
    d2rC2C(1, 3, IBC_UDIRICHLET) = c1
    d2rC2C(1, 4, IBC_UDIRICHLET) = d1

    !C2C i = 2
    d2fC2C(2, 1, IBC_UDIRICHLET) = alpha2
    d2fC2C(2, 2, IBC_UDIRICHLET) = ONE
    d2fC2C(2, 3, IBC_UDIRICHLET) = alpha2

    d2rC2C(2, 1, IBC_UDIRICHLET) = a2
    d2rC2C(2, 2, IBC_UDIRICHLET) = ZERO ! not used
    d2rC2C(2, 3, IBC_UDIRICHLET) = ZERO ! not used
    d2rC2C(2, 4, IBC_UDIRICHLET) = ZERO ! not used

    !C2C i = bulk
    d2fC2C(3, 1, IBC_UDIRICHLET) = alpha
    d2fC2C(3, 2, IBC_UDIRICHLET) = ONE
    d2fC2C(3, 3, IBC_UDIRICHLET) = alpha
    d2rC2C(3, 1, IBC_UDIRICHLET) = a / ONE
    d2rC2C(3, 2, IBC_UDIRICHLET) = b / FOUR
    d2rC2C(3, 3, IBC_UDIRICHLET) = ZERO ! not used
    d2rC2C(3, 4, IBC_UDIRICHLET) = ZERO ! not used

    !C2C i = n - 1
    d2fC2C(4, 1, IBC_UDIRICHLET) = alpha2
    d2fC2C(4, 2, IBC_UDIRICHLET) = ONE
    d2fC2C(4, 3, IBC_UDIRICHLET) = alpha2

    d2rC2C(4, 1, IBC_UDIRICHLET) = a2
    d2rC2C(4, 2, IBC_UDIRICHLET) = ZERO ! not used
    d2rC2C(4, 3, IBC_UDIRICHLET) = ZERO ! not used
    d2rC2C(4, 4, IBC_UDIRICHLET) = ZERO ! not used

    !C2C i = n
    d2fC2C(5, 1, IBC_UDIRICHLET) = alpha1
    d2fC2C(5, 2, IBC_UDIRICHLET) = ONE
    d2fC2C(5, 3, IBC_UDIRICHLET) = ZERO ! not used

    d2rC2C(5, 1, IBC_UDIRICHLET) = a1
    d2rC2C(5, 2, IBC_UDIRICHLET) = b1
    d2rC2C(5, 3, IBC_UDIRICHLET) = c1
    d2rC2C(5, 4, IBC_UDIRICHLET) = d1

    ! P2P
    d2fP2P(:, :, :) = d2fC2C(:, :, :)
    d2rP2P(:, :, :) = d2rC2C(:, :, :)
    call Print_debug_end_msg
    return
  end subroutine Assign_TDMA_coeffs
!===============================================================================
!===============================================================================
!> \brief Assigning the sparse matrix in the LHS of the compact scheme, and
!> calculating the geometry-only dependent variables for the TDMA scheme.
!>
!> This subroutine is called once locally.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n             the number of unknown array
!> \param[in]     bc            the boundary condition at two ends of the unknown
!> \param[in]     coeff         the basic TDMA coefficients defined above.
!> \param[out]    a             the coefficients for TDMA
!> \param[out]    b             a_i * x_(i-1) + b_i * x_(i) + c_i * x_(i+1)
!> \param[out]    c             = RHS
!> \param[out]    d             An assisting coeffients for the TDMA scheme.
!_______________________________________________________________________________
  subroutine Buildup_TDMA_LHS_array(n, bc, coeff, a, b, c, d)
!===============================================================================
! Module files
!===============================================================================
    use input_general_mod, only : IBC_PERIODIC
    use tridiagonal_matrix_algorithm
    implicit none
!===============================================================================
! Arguments
!===============================================================================
    integer(4), intent(in) :: n
    integer(4), intent(in) :: bc(2)
    real(WP), intent(in)   :: coeff(5, 3, 4)
    real(WP), intent(out)  :: a(n), b(n), c(n), d(n)
!===============================================================================
! Code
!===============================================================================
    a(1)         = coeff( 1, 1, bc(1) )
    a(2)         = coeff( 2, 1, bc(1) )
    a(3 : n - 2) = coeff( 3, 1, bc(1) )
    a(n - 1)     = coeff( 4, 1, bc(2) )
    a(n)         = coeff( 5, 1, bc(2) )

    b(1)         = coeff( 1, 2, bc(1) )
    b(2)         = coeff( 2, 2, bc(1) )
    b(3 : n - 2) = coeff( 3, 2, bc(1) )
    b(n - 1)     = coeff( 4, 2, bc(2) )
    b(n)         = coeff( 5, 2, bc(2) )

    c(1)         = coeff( 1, 3, bc(1) )
    c(2)         = coeff( 2, 3, bc(1) )
    c(3 : n - 2) = coeff( 3, 3, bc(1) )
    c(n - 1)     = coeff( 4, 3, bc(2) )
    c(n)         = coeff( 5, 3, bc(2) )

    if (bc(1) == IBC_PERIODIC) then
      call Preprocess_TDMA_coeffs(a(1:n-1), b(1:n-1), c(1:n-1), d(1:n-1), n-1)
    else 
      call Preprocess_TDMA_coeffs(a(:), b(:), c(:), d(:), n)
    end if 

    return
  end subroutine Buildup_TDMA_LHS_array
!===============================================================================
!===============================================================================
!> \brief Preparing the LHS matrix for the TDMA algorithm for compact scheme.
!>
!> This subroutine is called once locally.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!_______________________________________________________________________________
  subroutine Prepare_TDMA_LHS_matrix(d)
!===============================================================================
! Module files
!===============================================================================
    use tridiagonal_matrix_algorithm
    use input_general_mod, only : IBC_PERIODIC
    use udf_type_mod
    use parameters_constant_mod, only : ZERO
    implicit none

    type(t_domain), intent(in) :: d

    integer(4) :: i, nsz

    call Print_debug_start_msg &
         ("Preparing geometric matrix for the compact FD ...")
!-------------------------------------------------------------------------------
! 1st derivative in x direction
!-------------------------------------------------------------------------------
    i = 1
    ! 1st derivative in x direction with nc unknows
    nsz = d%nc(i)
    allocate (ad1x_C2C ( nsz ) ); ad1x_C2C(:) = ZERO
    allocate (bd1x_C2C ( nsz ) ); bd1x_C2C(:) = ZERO
    allocate (cd1x_C2C ( nsz ) ); cd1x_C2C(:) = ZERO
    allocate (dd1x_C2C ( nsz ) ); dd1x_C2C(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), d1fC2C, &
        ad1x_C2C, bd1x_C2C, cd1x_C2C, dd1x_C2C)

    allocate (ad1x_P2C ( nsz ) ); ad1x_P2C(:) = ZERO
    allocate (bd1x_P2C ( nsz ) ); bd1x_P2C(:) = ZERO
    allocate (cd1x_P2C ( nsz ) ); cd1x_P2C(:) = ZERO
    allocate (dd1x_P2C ( nsz ) ); dd1x_P2C(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), d1fP2C, &
        ad1x_P2C, bd1x_P2C, cd1x_P2C, dd1x_P2C)
    
    ! 1st derivative in x direction with np unknows
    nsz = d%np(i)    
    allocate (ad1x_P2P ( nsz ) ); ad1x_P2P(:) = ZERO
    allocate (bd1x_P2P ( nsz ) ); bd1x_P2P(:) = ZERO
    allocate (cd1x_P2P ( nsz ) ); cd1x_P2P(:) = ZERO
    allocate (dd1x_P2P ( nsz ) ); dd1x_P2P(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), d1fP2P, &
        ad1x_P2P, bd1x_P2P, cd1x_P2P, dd1x_P2P)

    allocate (ad1x_C2P ( nsz ) ); ad1x_C2P(:) = ZERO
    allocate (bd1x_C2P ( nsz ) ); bd1x_C2P(:) = ZERO
    allocate (cd1x_C2P ( nsz ) ); cd1x_C2P(:) = ZERO
    allocate (dd1x_C2P ( nsz ) ); dd1x_C2P(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), d1fC2P, &
        ad1x_C2P, bd1x_C2P, cd1x_C2P, dd1x_C2P)
!-------------------------------------------------------------------------------
! 1st derivative in y direction
!-------------------------------------------------------------------------------
    i = 2
    ! 1st derivative in y direction with nc unknows
    nsz = d%nc(i)
    allocate (ad1y_C2C ( nsz ) ); ad1y_C2C(:) = ZERO
    allocate (bd1y_C2C ( nsz ) ); bd1y_C2C(:) = ZERO
    allocate (cd1y_C2C ( nsz ) ); cd1y_C2C(:) = ZERO
    allocate (dd1y_C2C ( nsz ) ); dd1y_C2C(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), d1fC2C, &
        ad1y_C2C, bd1y_C2C, cd1y_C2C, dd1y_C2C)

    allocate (ad1y_P2C ( nsz ) ); ad1y_P2C(:) = ZERO
    allocate (bd1y_P2C ( nsz ) ); bd1y_P2C(:) = ZERO
    allocate (cd1y_P2C ( nsz ) ); cd1y_P2C(:) = ZERO
    allocate (dd1y_P2C ( nsz ) ); dd1y_P2C(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), d1fP2C, &
        ad1y_P2C, bd1y_P2C, cd1y_P2C, dd1y_P2C)

    ! 1st derivative in y direction with np unknows
    nsz = d%np(i)
    allocate (ad1y_P2P ( nsz ) ); ad1y_P2P(:) = ZERO
    allocate (bd1y_P2P ( nsz ) ); bd1y_P2P(:) = ZERO
    allocate (cd1y_P2P ( nsz ) ); cd1y_P2P(:) = ZERO
    allocate (dd1y_P2P ( nsz ) ); dd1y_P2P(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), d1fP2P, &
        ad1y_P2P, bd1y_P2P, cd1y_P2P, dd1y_P2P)

    allocate (ad1y_C2P ( nsz ) ); ad1y_C2P(:) = ZERO
    allocate (bd1y_C2P ( nsz ) ); bd1y_C2P(:) = ZERO
    allocate (cd1y_C2P ( nsz ) ); cd1y_C2P(:) = ZERO
    allocate (dd1y_C2P ( nsz ) ); dd1y_C2P(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), d1fC2P, &
        ad1y_C2P, bd1y_C2P, cd1y_C2P, dd1y_C2P) 
!-------------------------------------------------------------------------------
! 1st derivative in z direction
!-------------------------------------------------------------------------------
    i = 3
    ! 1st derivative in z direction with nc unknows
    nsz = d%nc(i)
    allocate (ad1z_C2C ( nsz ) ); ad1z_C2C(:) = ZERO
    allocate (bd1z_C2C ( nsz ) ); bd1z_C2C(:) = ZERO
    allocate (cd1z_C2C ( nsz ) ); cd1z_C2C(:) = ZERO
    allocate (dd1z_C2C ( nsz ) ); dd1z_C2C(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), d1fC2C, &
        ad1z_C2C, bd1z_C2C, cd1z_C2C, dd1z_C2C)

    allocate (ad1z_P2C ( nsz ) ); ad1z_P2C(:) = ZERO
    allocate (bd1z_P2C ( nsz ) ); bd1z_P2C(:) = ZERO
    allocate (cd1z_P2C ( nsz ) ); cd1z_P2C(:) = ZERO
    allocate (dd1z_P2C ( nsz ) ); dd1z_P2C(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), d1fP2C, &
        ad1z_P2C, bd1z_P2C, cd1z_P2C, dd1z_P2C)
    ! 1st derivative in z direction with np unknows
    nsz = d%np(i)
    allocate (ad1z_P2P ( nsz ) ); ad1z_P2P(:) = ZERO
    allocate (bd1z_P2P ( nsz ) ); bd1z_P2P(:) = ZERO
    allocate (cd1z_P2P ( nsz ) ); cd1z_P2P(:) = ZERO
    allocate (dd1z_P2P ( nsz ) ); dd1z_P2P(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), d1fP2P, &
        ad1z_P2P, bd1z_P2P, cd1z_P2P, dd1z_P2P)

    allocate (ad1z_C2P ( nsz ) ); ad1z_C2P(:) = ZERO
    allocate (bd1z_C2P ( nsz ) ); bd1z_C2P(:) = ZERO
    allocate (cd1z_C2P ( nsz ) ); cd1z_C2P(:) = ZERO
    allocate (dd1z_C2P ( nsz ) ); dd1z_C2P(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), d1fC2P, &
        ad1z_C2P, bd1z_C2P, cd1z_C2P, dd1z_C2P)

!-------------------------------------------------------------------------------
! mid-point interpolation in x direction with nc unknows
!-------------------------------------------------------------------------------
    i = 1
    nsz = d%nc(i)
    allocate (am1x_P2C ( nsz ) ); am1x_P2C(:) = ZERO
    allocate (bm1x_P2C ( nsz ) ); bm1x_P2C(:) = ZERO
    allocate (cm1x_P2C ( nsz ) ); cm1x_P2C(:) = ZERO
    allocate (dm1x_P2C ( nsz ) ); dm1x_P2C(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), m1fP2C, &
        am1x_P2C, bm1x_P2C, cm1x_P2C, dm1x_P2C)
!-------------------------------------------------------------------------------
! mid-point interpolation in x direction with np unknows
!-------------------------------------------------------------------------------
    nsz = d%np(i)
    allocate (am1x_C2P ( nsz ) ); am1x_C2P(:) = ZERO
    allocate (bm1x_C2P ( nsz ) ); bm1x_C2P(:) = ZERO
    allocate (cm1x_C2P ( nsz ) ); cm1x_C2P(:) = ZERO
    allocate (dm1x_C2P ( nsz ) ); dm1x_C2P(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), m1fC2P, &
        am1x_C2P, bm1x_C2P, cm1x_C2P, dm1x_C2P)
!-------------------------------------------------------------------------------
! mid-point interpolation in y direction with nc unknows
!-------------------------------------------------------------------------------
    i = 2
    nsz = d%nc(i)
    allocate (am1y_P2C ( nsz ) ); am1y_P2C(:) = ZERO
    allocate (bm1y_P2C ( nsz ) ); bm1y_P2C(:) = ZERO
    allocate (cm1y_P2C ( nsz ) ); cm1y_P2C(:) = ZERO
    allocate (dm1y_P2C ( nsz ) ); dm1y_P2C(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), m1fP2C, &
        am1y_P2C, bm1y_P2C, cm1y_P2C, dm1y_P2C)
!-------------------------------------------------------------------------------
! mid-point interpolation in y direction with np unknows
!-------------------------------------------------------------------------------
    nsz = d%np(i)
    allocate (am1y_C2P ( nsz ) ); am1y_C2P(:) = ZERO
    allocate (bm1y_C2P ( nsz ) ); bm1y_C2P(:) = ZERO
    allocate (cm1y_C2P ( nsz ) ); cm1y_C2P(:) = ZERO
    allocate (dm1y_C2P ( nsz ) ); dm1y_C2P(:) = ZERO
    call Buildup_TDMA_LHS_array(nsz, d%bc(:, i), m1fC2P, &
        am1y_C2P, bm1y_C2P, cm1y_C2P, dm1y_C2P)
!-------------------------------------------------------------------------------
! mid-point interpolation in z direction with nc unknows
!-------------------------------------------------------------------------------
    i = 3
    nsz = d%nc(i)
    allocate (am1z_P2C ( nsz ) ); am1z_P2C(:) = ZERO
    allocate (bm1z_P2C ( nsz ) ); bm1z_P2C(:) = ZERO
    allocate (cm1z_P2C ( nsz ) ); cm1z_P2C(:) = ZERO
    allocate (dm1z_P2C ( nsz ) ); dm1z_P2C(:) = ZERO
    call Buildup_TDMA_LHS_array( nsz, d%bc(:, i), m1fP2C, &
        am1z_P2C, bm1z_P2C, cm1z_P2C, dm1z_P2C)
!-------------------------------------------------------------------------------
! mid-point interpolation in z direction with np unknows
!-------------------------------------------------------------------------------
    nsz = d%np(i)
    allocate (am1z_C2P ( nsz ) ); am1z_C2P(:) = ZERO
    allocate (bm1z_C2P ( nsz ) ); bm1z_C2P(:) = ZERO
    allocate (cm1z_C2P ( nsz ) ); cm1z_C2P(:) = ZERO
    allocate (dm1z_C2P ( nsz ) ); dm1z_C2P(:) = ZERO
    call Buildup_TDMA_LHS_array( nsz, d%bc(:, i), m1fC2P, &
        am1z_C2P, bm1z_C2P, cm1z_C2P, dm1z_C2P)
!-------------------------------------------------------------------------------
! 2nd order deriviative in x direction with nc unknows
!-------------------------------------------------------------------------------
    i = 1
    nsz = d%nc(i)
    allocate (ad2x_C2C ( nsz ) ); ad2x_C2C(:) = ZERO
    allocate (bd2x_C2C ( nsz ) ); bd2x_C2C(:) = ZERO
    allocate (cd2x_C2C ( nsz ) ); cd2x_C2C(:) = ZERO
    allocate (dd2x_C2C ( nsz ) ); dd2x_C2C(:) = ZERO
    call Buildup_TDMA_LHS_array( nsz, d%bc(:, i), d2fC2C, &
        ad2x_C2C, bd2x_C2C, cd2x_C2C, dd2x_C2C)

!-------------------------------------------------------------------------------
! 2nd order deriviative in x direction with np unknows
!-------------------------------------------------------------------------------
    nsz = d%np(i)
    allocate (ad2x_P2P ( nsz ) ); ad2x_P2P(:) = ZERO
    allocate (bd2x_P2P ( nsz ) ); bd2x_P2P(:) = ZERO
    allocate (cd2x_P2P ( nsz ) ); cd2x_P2P(:) = ZERO
    allocate (dd2x_P2P ( nsz ) ); dd2x_P2P(:) = ZERO
    call Buildup_TDMA_LHS_array( nsz, d%bc(:, i), d2fP2P, &
        ad2x_P2P, bd2x_P2P, cd2x_P2P, dd2x_P2P)
!-------------------------------------------------------------------------------
! 2nd order deriviative in y direction with nc unknows
!-------------------------------------------------------------------------------
    i = 2
    nsz = d%nc(i)
    allocate (ad2y_C2C ( nsz ) ); ad2y_C2C(:) = ZERO
    allocate (bd2y_C2C ( nsz ) ); bd2y_C2C(:) = ZERO
    allocate (cd2y_C2C ( nsz ) ); cd2y_C2C(:) = ZERO
    allocate (dd2y_C2C ( nsz ) ); dd2y_C2C(:) = ZERO
    call Buildup_TDMA_LHS_array( nsz, d%bc(:, i), d2fC2C, &
        ad2y_C2C, bd2y_C2C, cd2y_C2C, dd2y_C2C)
!-------------------------------------------------------------------------------
! 2nd order deriviative in y direction with np unknows
!-------------------------------------------------------------------------------
    nsz = d%np(i)
    allocate (ad2y_P2P ( nsz ) ); ad2y_P2P(:) = ZERO
    allocate (bd2y_P2P ( nsz ) ); bd2y_P2P(:) = ZERO
    allocate (cd2y_P2P ( nsz ) ); cd2y_P2P(:) = ZERO
    allocate (dd2y_P2P ( nsz ) ); dd2y_P2P(:) = ZERO
    call Buildup_TDMA_LHS_array( nsz, d%bc(:, i), d2fP2P, &
        ad2y_P2P, bd2y_P2P, cd2y_P2P, dd2y_P2P)
!-------------------------------------------------------------------------------
! 2nd order deriviative in z direction with nc unknows
!-------------------------------------------------------------------------------
    i = 3
    nsz = d%nc(i)
    allocate (ad2z_C2C ( nsz ) ); ad2z_C2C(:) = ZERO
    allocate (bd2z_C2C ( nsz ) ); bd2z_C2C(:) = ZERO
    allocate (cd2z_C2C ( nsz ) ); cd2z_C2C(:) = ZERO
    allocate (dd2z_C2C ( nsz ) ); dd2z_C2C(:) = ZERO
    call Buildup_TDMA_LHS_array( nsz, d%bc(:, i), d2fC2C, &
        ad2z_C2C, bd2z_C2C, cd2z_C2C, dd2z_C2C)

!-------------------------------------------------------------------------------
! 2nd order deriviative in z direction with np unknows
!-------------------------------------------------------------------------------
    nsz = d%np(i)
    allocate (ad2z_P2P ( nsz ) ); ad2z_P2P(:) = ZERO
    allocate (bd2z_P2P ( nsz ) ); bd2z_P2P(:) = ZERO
    allocate (cd2z_P2P ( nsz ) ); cd2z_P2P(:) = ZERO
    allocate (dd2z_P2P ( nsz ) ); dd2z_P2P(:) = ZERO
    call Buildup_TDMA_LHS_array( nsz, d%bc(:, i), d2fP2P, &
        ad2z_P2P, bd2z_P2P, cd2z_P2P, dd2z_P2P)

    call Print_debug_end_msg
    return
  end subroutine Prepare_TDMA_LHS_matrix
!===============================================================================
!===============================================================================
  subroutine Prepare_coeffs_for_operations
    use input_general_mod, only : iAccuracy
    use geometry_mod,      only : domain
    implicit none

    call Assign_TDMA_coeffs (iAccuracy)
    call Prepare_TDMA_LHS_matrix (domain)
    return
  end subroutine Prepare_coeffs_for_operations
!===============================================================================
!===============================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for interpolation.
!>
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> for the interpolation.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     str           string to flag interpolation type, C2P or P2C
!> \param[in]     n             the number of unknowns
!> \param[in]     bc            the b.c. at two ends of the unknown array
!> \param[in]     inbr          the neibouring index of unknown index
!> \param[in]     coeff         the defined TDMA coefficients
!> \param[in]     fi            the input variable to build up the RHS array
!> \param[out]    fo            the output RHS array
!_______________________________________________________________________________
  subroutine Prepare_TDMA_interp_RHS_array(str, n, bc, inbr, coeff, fi, fo)
!===============================================================================
! Module files
!===============================================================================
    use parameters_constant_mod
    use input_general_mod
    implicit none
!===============================================================================
! Arguments
!===============================================================================
    character(3), intent(in) :: str
    integer(4),   intent(in) :: n ! unknow numbers
    integer(4),   intent(in) :: bc(2)
    integer(4),   intent(in) :: inbr(:, :)
    real(WP),     intent(in) :: coeff(5, 3, 4)
    real(WP),     intent(in) :: fi(:)
    real(WP),     intent(out):: fo(n)
!===============================================================================
! Local arguments
!===============================================================================
    integer(4) :: i
    integer(4) :: im2, im1, ip1, ip2
    logical :: fbc = .false.
    real(WP) :: fsign
!===============================================================================
! Code
!===============================================================================
    ! initilisation
    fo(:) = ZERO
!-------------------------------------------------------------------------------
! bulk body for i from 2 to n-1
! bulk body and b.c. for periodic b.c.
!-------------------------------------------------------------------------------
    do i = 1, n
      ! exclude non-periodic b.c. at both sides
      fbc = (i == 1 .or. i == 2 .or. i == n-1 .or. i==n)
      if( (.not. bc(1)==IBC_PERIODIC) .and. fbc) cycle

      im2 = inbr(1, i)
      im1 = inbr(2, i)
      ip1 = inbr(3, i)
      ip2 = inbr(4, i)

      if (str == 'P2C') then
        fo(i) = coeff( 3, 1, bc(1) ) * ( fi(ip1) + fi(i) ) + &
                coeff( 3, 2, bc(1) ) * ( fi(ip2) + fi(im1) )
      else if (str == 'C2P') then
        fo(i) = coeff( 3, 1, bc(1) ) * ( fi(i) + fi(im1) ) + &
                coeff( 3, 2, bc(1) ) * ( fi(ip1) + fi(im2) )
      else 
        call Print_error_msg("Error input in prepare_FD_TDMA_RHS in Subroutine: " // &
        "Prepare_TDMA_interp_RHS_array")
      end if

    end do
!-------------------------------------------------------------------------------
! boundary at the side of i = 1
!-------------------------------------------------------------------------------
    if (bc(1) == IBC_PERIODIC) then
      ! do nothing
    else if (bc(1) == IBC_ASYMMETRIC .or. bc(1) == IBC_SYMMETRIC) then
      if (bc(1) == IBC_ASYMMETRIC) then
        fsign = - ONE
      else 
        fsign = ONE
      end if

      if (str == 'P2C') then
        i = 1
        im2 = 3 ! -1'
        im1 = 2 ! 0'
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(1) ) * ( fi(ip1) +         fi(i) ) + &
                coeff( i, 2, bc(1) ) * ( fi(ip2) + fsign * fi(im1) )

        i = 2
        im2 = 2 ! 0'
        im1 = i - 1
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(1) ) * ( fi(ip1) + fi(i) ) + &
                coeff( i, 2, bc(1) ) * ( fi(ip2) + fi(im1) )

      else if (str == 'C2P') then

        i = 1
        im2 = 2 ! -1
        im1 = 1 ! 0
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(1) ) * ( fi(i)   + fsign * fi(im1) ) + &
                coeff( i, 2, bc(1) ) * ( fi(ip1) + fsign * fi(im2) )

        i = 2
        im2 = 1 ! 0
        im1 = i - 1
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(1) ) * ( fi(i)   +         fi(im1) ) + &
                coeff( i, 2, bc(1) ) * ( fi(ip1) + fsign * fi(im2) )

      else 
        call Print_error_msg("Error input in prepare_FD_TDMA_RHS in Subroutine: " // &
        "Prepare_TDMA_interp_RHS_array")
      end if

    else if (bc(1) == IBC_UDIRICHLET) then
      
      if (str == 'P2C') then

        fo(1) = coeff( 1, 1, bc(1) ) * fi(1) + &
                coeff( 1, 2, bc(1) ) * fi(2)  + &
                coeff( 1, 3, bc(1) ) * fi(3) 
        fo(2) = coeff( 2, 1, bc(1) ) * ( fi(3) + fi(2) )

      else if (str == 'C2P') then

        fo(1) = coeff( 1, 1, bc(1) ) * fi(1) + &
                coeff( 1, 2, bc(1) ) * fi(2)  + &
                coeff( 1, 3, bc(1) ) * fi(3) 
        fo(2) = coeff( 2, 1, bc(1) ) * ( fi(2) + fi(1) )
        
      else 
        call Print_error_msg("Error input in prepare_FD_TDMA_RHS in Subroutine: " // &
        "Prepare_TDMA_interp_RHS_array")
      end if

    else 
      call Print_error_msg("No Such Boundary Defined in Subroutine: " // &
      "Prepare_TDMA_interp_RHS_array")
    end if
!-------------------------------------------------------------------------------
! boundary at the side of i = n
!-------------------------------------------------------------------------------
    if (bc(2) == IBC_PERIODIC) then
      ! do nothing
    else if (bc(2) == IBC_ASYMMETRIC .or. bc(2) == IBC_SYMMETRIC) then
      if (bc(2) == IBC_ASYMMETRIC) then
        fsign = - ONE
      else 
        fsign = ONE
      end if

      if (str == 'P2C') then
        i = n
        im2 = i - 2
        im1 = i - 1
        ip1 = i + 1 
        ip2 = n     ! n' + 2
        fo(i) = coeff( i, 1, bc(2) ) * (         fi(ip1) + fi(i) ) + &
                coeff( i, 2, bc(2) ) * ( fsign * fi(ip2) + fi(im1) )

        i = n - 1
        im2 = i - 2
        im1 = i - 1
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(2) ) * ( fi(ip1) + fi(i) ) + &
                coeff( i, 2, bc(2) ) * ( fi(ip2) + fi(im1) )

      else if (str == 'C2P') then

        i = n
        im2 = i - 2
        im1 = i - 1
        ip1 = n     ! n + 1
        ip2 = n - 1 ! n + 2
        fo(i) = coeff( i, 1, bc(2) ) * (         fi(i)   + fi(im1) ) + &
                coeff( i, 2, bc(2) ) * ( fsign * fi(ip1) + fi(im2) )

        i = n - 1
        im2 = i - 2
        im1 = i - 1
        ip1 = i + 1
        ip2 = n ! n + 1
        fo(i) = coeff( i, 1, bc(2) ) * ( fi(i)   + fi(im1) ) + &
                coeff( i, 2, bc(2) ) * ( fi(ip1) + fi(im2) )

      else 
        call Print_error_msg("Error input in prepare_FD_TDMA_RHS in Subroutine: " // &
        "Prepare_TDMA_interp_RHS_array")
      end if

    else if (bc(2) == IBC_UDIRICHLET) then
      
      if (str == 'P2C') then

        fo(n) = coeff( 5, 1, bc(2) ) * fi(n + 1) + &
                coeff( 5, 2, bc(2) ) * fi(n    ) + &
                coeff( 5, 3, bc(2) ) * fi(n - 1) 
        fo(n - 1) = coeff( 4, 1, bc(2) ) * ( fi(n) + fi(n - 1) )

      else if (str == 'C2P') then

        fo(n) = coeff( 5, 1, bc(2) ) * fi(n - 1) + &
                coeff( 5, 2, bc(2) ) * fi(n - 2) + &
                coeff( 5, 3, bc(2) ) * fi(n - 3) 
        fo(n - 1) = coeff( 4, 1, bc(2) ) * ( fi(n - 1) + fi(n - 2) )
        
      else 
        call Print_error_msg("Error input in prepare_FD_TDMA_RHS in Subroutine: " // &
        "Prepare_TDMA_interp_RHS_array")
      end if

    else 
      call Print_error_msg("No Such Boundary Defined in Subroutine: " // &
      "Prepare_TDMA_interp_RHS_array")
    end if

    return
  end subroutine Prepare_TDMA_interp_RHS_array
!===============================================================================
!===============================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for 1st derivative.
!>
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> for the 1st derivative.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     str           string to flag interpolation type
!>                              C2C, P2P, C2P, P2C
!> \param[in]     n             the number of unknowns
!> \param[in]     bc            the b.c. at two ends of the unknown array
!> \param[in]     inbr          the neibouring index of unknown index
!> \param[in]     coeff         the defined TDMA coefficients
!> \param[in]     fi            the input variable to build up the RHS array
!> \param[out]    fo            the output RHS array
!_______________________________________________________________________________
  subroutine Prepare_TDMA_1deri_RHS_array(str, n, bc, inbr, dd, coeff, fi, fo)
!===============================================================================
! Module files
!===============================================================================
    use parameters_constant_mod
    use input_general_mod
    implicit none
!===============================================================================
! Arguments
!===============================================================================
    character(3), intent(in) :: str
    integer(4),   intent(in) :: n ! unknow numbers
    integer(4),   intent(in) :: bc(2)
    integer(4),   intent(in) :: inbr(:, :)
    real(WP),     intent(in) :: dd
    real(WP),     intent(in) :: coeff(5, 3, 4)
    real(WP),     intent(in) :: fi(:)
    real(WP),     intent(out):: fo(n)
!===============================================================================
! Local arguments
!===============================================================================    
    integer(4) :: i
    integer(4) :: im2, im1, ip1, ip2
    logical :: fbc = .false.
    real(WP) :: fsign
!===============================================================================
! Code
!===============================================================================
    ! initilisation
    fo(:) = ZERO
!-------------------------------------------------------------------------------
! bulk body for i from 2 to n-1
! bulk body and b.c. for periodic b.c.
!-------------------------------------------------------------------------------
    do i = 1, n

      ! exclude non-periodic b.c. at both sides
      fbc = (i == 1 .or. i == 2 .or. i == n-1 .or. i == n)
      if( (.not. bc(1)==IBC_PERIODIC) .and. fbc) cycle

      im2 = inbr(1, i)
      im1 = inbr(2, i)
      ip1 = inbr(3, i)
      ip2 = inbr(4, i)

      if (str == 'C2C' .or. str == 'P2P') then
        fo(i) = coeff( 3, 1, bc(1) ) * ( fi(ip1) - fi(im1) ) + &
                coeff( 3, 2, bc(1) ) * ( fi(ip2) - fi(im2) )
      else if (str == 'P2C') then
        fo(i) = coeff( 3, 1, bc(1) ) * ( fi(ip1) - fi(i  ) ) + &
                coeff( 3, 2, bc(1) ) * ( fi(ip2) - fi(im1) )
      else if (str == 'C2P') then
        fo(i) = coeff( 3, 1, bc(1) ) * ( fi(i  ) - fi(im1) ) + &
                coeff( 3, 2, bc(1) ) * ( fi(ip1) - fi(im2) )
      else 
        call Print_error_msg("101: Error input in prepare_FD_TDMA_RHS.")
      end if

    end do
!-------------------------------------------------------------------------------
! boundary at the side of i = 1
!-------------------------------------------------------------------------------
    if (bc(1) == IBC_PERIODIC) then
      ! do nothing
    else if (bc(1) == IBC_ASYMMETRIC .or. bc(1) == IBC_SYMMETRIC) then
      if (bc(1) == IBC_ASYMMETRIC) then
        fsign = - ONE
      else 
        fsign = ONE
      end if

      if (str == 'C2C') then

        i = 1
        im2 = 2 ! -1
        im1 = 1 ! 0
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(1) ) * ( fi(ip1) - fsign * fi(im1) ) + &
                coeff( i, 2, bc(1) ) * ( fi(ip2) - fsign * fi(im2) )

        i = 2
        im2 = 1 ! 0
        im1 = i - 1
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(1) ) * ( fi(ip1) -         fi(im1) ) + &
                coeff( i, 2, bc(1) ) * ( fi(ip2) - fsign * fi(im2) )

      else if (str == 'P2P') then

        i = 1
        im2 = 3 ! -1'
        im1 = 2 ! 0'
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(1) ) * ( fi(ip1) - fsign * fi(im1) ) + &
                coeff( i, 2, bc(1) ) * ( fi(ip2) - fsign * fi(im2) )

        i = 2
        im2 = 2 ! 0'
        im1 = i - 1
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(1) ) * ( fi(ip1) -         fi(im1) ) + &
                coeff( i, 2, bc(1) ) * ( fi(ip2) - fsign * fi(im2) )

      else if (str == 'P2C') then
        i = 1
        im2 = 3 ! -1'
        im1 = 2 ! 0'
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(1) ) * ( fi(ip1) -         fi(i) ) + &
                coeff( i, 2, bc(1) ) * ( fi(ip2) - fsign * fi(im1) )

        i = 2
        im2 = 2 ! 0'
        im1 = i - 1
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(1) ) * ( fi(ip1) - fi(i) ) + &
                coeff( i, 2, bc(1) ) * ( fi(ip2) - fi(im1) )

      else if (str == 'C2P') then

        i = 1
        im2 = 2 ! -1
        im1 = 1 ! 0
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(1) ) * ( fi(i)   - fsign * fi(im1) ) + &
                coeff( i, 2, bc(1) ) * ( fi(ip1) - fsign * fi(im2) )

        i = 2
        im2 = 1 ! 0
        im1 = i - 1
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(1) ) * ( fi(i) -           fi(im1) ) + &
                coeff( i, 2, bc(1) ) * ( fi(ip1) - fsign * fi(im2) )

      else 
        call Print_error_msg("102: Error input in prepare_FD_TDMA_RHS.")
      end if

    else if (bc(1) == IBC_UDIRICHLET) then
      
      if (str == 'C2C' .or. str == 'P2P') then

        fo(1) = coeff( 1, 1, bc(1) ) * fi(1) + &
                coeff( 1, 2, bc(1) ) * fi(2) + &
                coeff( 1, 3, bc(1) ) * fi(3) 
        fo(2) = coeff( 2, 1, bc(1) ) * ( fi(3) - fi(1) )

      else if (str == 'P2C') then

        fo(1) = coeff( 1, 1, bc(1) ) * fi(1) + &
                coeff( 1, 2, bc(1) ) * fi(2)  + &
                coeff( 1, 3, bc(1) ) * fi(3) 
        fo(2) = coeff( 2, 1, bc(1) ) * ( fi(3) - fi(2) )

      else if (str == 'C2P') then

        fo(1) = coeff( 1, 1, bc(1) ) * fi(1) + &
                coeff( 1, 2, bc(1) ) * fi(2)  + &
                coeff( 1, 3, bc(1) ) * fi(3) 
        fo(2) = coeff( 2, 1, bc(1) ) * ( fi(2) - fi(1) )
        
      else 
        call Print_error_msg("103: Error input in prepare_FD_TDMA_RHS.")
      end if

    else 
      call Print_error_msg("104: No Such Boundary Defined.")
    end if
!-------------------------------------------------------------------------------
! boundary at the side of i = n
!-------------------------------------------------------------------------------
    if (bc(2) == IBC_PERIODIC) then
      ! do nothing
    else if (bc(2) == IBC_ASYMMETRIC .or. bc(2) == IBC_SYMMETRIC) then
      if (bc(2) == IBC_ASYMMETRIC) then
        fsign = - ONE
      else 
        fsign = ONE
      end if

      if (str == 'C2C') then

        i = n
        im2 = i - 2
        im1 = i - 1
        ip1 = n     ! n + 1
        ip2 = n - 1 ! n + 2
        fo(i) = coeff( i, 1, bc(2) ) * ( fsign * fi(ip1) - fi(im1) ) + &
                coeff( i, 2, bc(2) ) * ( fsign * fi(ip2) - fi(im2) )

        i = n - 1
        im2 = i - 2
        im1 = i - 1
        ip1 = i + 1
        ip2 = n ! n + 1
        fo(i) = coeff( i, 1, bc(2) ) * (         fi(ip1) - fi(im1) ) + &
                coeff( i, 2, bc(2) ) * ( fsign * fi(ip2) - fi(im2) )

      else if (str == 'P2P') then

        i = n
        im2 = i - 2
        im1 = i - 1
        ip1 = i + 1 
        ip2 = n     ! n' + 2
        fo(i) = coeff( i, 1, bc(2) ) * (         fi(ip1) - fi(im1) ) + &
                coeff( i, 2, bc(2) ) * ( fsign * fi(ip2) - fi(im2) )

        i = n - 1
        im2 = i - 2
        im1 = i - 1
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(2) ) * ( fi(ip1) - fi(im1) ) + &
                coeff( i, 2, bc(2) ) * ( fi(ip2) - fi(im2) )

      else if (str == 'P2C') then
        i = n
        im2 = i - 2
        im1 = i - 1
        ip1 = i + 1 
        ip2 = n     ! n' + 2
        fo(i) = coeff( i, 1, bc(2) ) * (         fi(ip1) - fi(i) ) + &
                coeff( i, 2, bc(2) ) * ( fsign * fi(ip2) - fi(im1) )

        i = n - 1
        im2 = i - 2
        im1 = i - 1
        ip1 = i + 1
        ip2 = i + 2
        fo(i) = coeff( i, 1, bc(2) ) * ( fi(ip1) - fi(i) ) + &
                coeff( i, 2, bc(2) ) * ( fi(ip2) - fi(im1) )

      else if (str == 'C2P') then

        i = n
        im2 = i - 2
        im1 = i - 1
        ip1 = n     ! n + 1
        ip2 = n - 1 ! n + 2
        fo(i) = coeff( i, 1, bc(2) ) * (         fi(i)   - fi(im1) ) + &
                coeff( i, 2, bc(2) ) * ( fsign * fi(ip1) - fi(im2) )

        i = n - 1
        im2 = i - 2
        im1 = i - 1
        ip1 = i + 1
        ip2 = n ! n + 1
        fo(i) = coeff( i, 1, bc(2) ) * ( fi(i)   - fi(im1) ) + &
                coeff( i, 2, bc(2) ) * ( fi(ip1) - fi(im2) )

      else 
        call Print_error_msg("105: Error input in prepare_FD_TDMA_RHS.")
      end if

    else if (bc(2) == IBC_UDIRICHLET) then
      
      if (str == 'C2C' .or. str == 'P2P') then

        fo(n) = coeff( 5, 1, bc(2) ) * fi(n) + &
                coeff( 5, 2, bc(2) ) * fi(n - 1) + &
                coeff( 5, 3, bc(2) ) * fi(n - 2) 
        fo(n - 1) = coeff( 4, 1, bc(2) ) * ( fi(n) - fi(n - 2) )

      else if (str == 'P2C') then

        fo(n) = coeff( 5, 1, bc(2) ) * fi(n + 1) + &
                coeff( 5, 2, bc(2) ) * fi(n    ) + &
                coeff( 5, 3, bc(2) ) * fi(n - 1) 
        fo(n - 1) = coeff( 4, 1, bc(2) ) * ( fi(n) - fi(n - 1) )

      else if (str == 'C2P') then

        fo(n) = coeff( 5, 1, bc(2) ) * fi(n - 1) + &
                coeff( 5, 2, bc(2) ) * fi(n - 2)  + &
                coeff( 5, 3, bc(2) ) * fi(n - 3) 
        fo(n - 1) = coeff( 4, 1, bc(2) ) * ( fi(n - 1) - fi(n - 2) )
        
      else 
        call Print_error_msg("106: Error input in prepare_FD_TDMA_RHS.")
      end if

    else 
      call Print_error_msg("107: Error input in prepare_FD_TDMA_RHS.")
    end if

    fo(:) = fo(:) * dd

    return
  end subroutine Prepare_TDMA_1deri_RHS_array
!===============================================================================
!===============================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for 2nd derivative.
!>
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> for the 2nd derivative.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     str           string to flag interpolation type
!>                              C2C, P2P
!> \param[in]     n             the number of unknowns
!> \param[in]     bc            the b.c. at two ends of the unknown array
!> \param[in]     inbr          the neibouring index of unknown index
!> \param[in]     coeff         the defined TDMA coefficients
!> \param[in]     fi            the input variable to build up the RHS array
!> \param[out]    fo            the output RHS array
!_______________________________________________________________________________
  subroutine Prepare_TDMA_2deri_RHS_array(str, n, bc, inbr, dd, coeff, fi, fo)
!===============================================================================
! Module files
!===============================================================================
    use parameters_constant_mod
    use input_general_mod
    implicit none
!===============================================================================
! Arguments
!===============================================================================
    character(3), intent(in) :: str
    integer(4),   intent(in) :: n ! unknow numbers
    integer(4),   intent(in) :: bc(2)
    integer(4),   intent(in) :: inbr(:, :)
    real(WP),     intent(in) :: dd
    real(WP),     intent(in) :: coeff(5, 4, 4)
    real(WP),     intent(in) :: fi(:)
    real(WP),     intent(out):: fo(n)
!===============================================================================
! Local arguments
!===============================================================================    
    integer(4) :: i
    integer(4) :: im2, im1, ip1, ip2
    logical :: fbc = .false.
!===============================================================================
! Code
!===============================================================================
    ! initilisation
    fo(:) = ZERO
!-------------------------------------------------------------------------------
! bulk body for i from 2 to n-1
! bulk body and b.c. for periodic b.c.
!-------------------------------------------------------------------------------
    do i = 1, n

      ! exclude non-periodic b.c. at both sides
      fbc = (i == 1 .or. i == 2 .or. i == n-1 .or. i == n)
      if( (.not. bc(1)==IBC_PERIODIC) .and. fbc) cycle

      im2 = inbr(1, i)
      im1 = inbr(2, i)
      ip1 = inbr(3, i)
      ip2 = inbr(4, i)

      if (str == 'C2C' .or. str == 'P2P') then
        fo(i) = coeff( 3, 1, bc(1) ) * ( fi(ip1) - TWO * fi(i) + fi(im1) ) + &
                coeff( 3, 2, bc(1) ) * ( fi(ip2) - TWO * fi(i) + fi(im2) )
      else 
        call Print_error_msg("101: Error input in prepare_FD_TDMA_RHS.")
      end if

    end do
!-------------------------------------------------------------------------------
! boundary at the side of i = 1
!-------------------------------------------------------------------------------
    if (bc(1) == IBC_PERIODIC) then
      ! do nothing
    else if (bc(1) == IBC_UDIRICHLET) then
      
      if (str == 'C2C' .or. str == 'P2P') then

        fo(1) = coeff( 1, 1, bc(1) ) * fi(1) + &
                coeff( 1, 2, bc(1) ) * fi(2) + &
                coeff( 1, 3, bc(1) ) * fi(3) + &
                coeff( 1, 4, bc(1) ) * fi(4)

        fo(2) = coeff( 2, 1, bc(1) ) * ( fi(3) - TWO * fi(2) + fi(1) )
      else 
        call Print_error_msg("103: Error input in prepare_FD_TDMA_RHS.")
      end if

    else 
      call Print_error_msg("104: No Such Boundary Defined.")
    end if
!-------------------------------------------------------------------------------
! boundary at the side of i = n
!-------------------------------------------------------------------------------
    if (bc(2) == IBC_PERIODIC) then
      ! do nothing
    else if (bc(2) == IBC_UDIRICHLET) then
      
      if (str == 'C2C' .or. str == 'P2P') then

        fo(n) = coeff( 5, 1, bc(2) ) * fi(n) + &
                coeff( 5, 2, bc(2) ) * fi(n - 1) + &
                coeff( 5, 3, bc(2) ) * fi(n - 2) + &
                coeff( 5, 4, bc(2) ) * fi(n - 3) 
        fo(n - 1) = coeff( 4, 1, bc(2) ) * ( fi(n) - TWO * fi(n - 1) + fi(n - 2) )

      else 
        call Print_error_msg("106: Error input in prepare_FD_TDMA_RHS.")
      end if

    else 
      call Print_error_msg("107: Error input in prepare_FD_TDMA_RHS.")
    end if

    fo(:) = fo(:) * dd

    return
  end subroutine Prepare_TDMA_2deri_RHS_array
!===============================================================================
!===============================================================================
!> \brief To caculate the mid-point interpolation in 1D.
!>
!> This subroutine is called as required to get the mid-point interpolation.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     str1          string to flag which direction to impletment
!> \param[in]     str2          string to flag C2P or P2C
!> \param[in]     d             domain
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!_______________________________________________________________________________
  subroutine Get_midp_interpolation_1D(str1, str2, d, fi, fo)
!===============================================================================
! Module files
!===============================================================================
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
!===============================================================================
! Arguments
!===============================================================================
    type(t_domain), intent(in) :: d
    character(1),   intent(in) :: str1
    character(3),   intent(in) :: str2
    real(WP),       intent(in) :: fi(:)
    real(WP),       intent(out):: fo(:)
!===============================================================================
! Local arguments
!===============================================================================
    integer(4) :: i, nsz
!===============================================================================
! Code
!===============================================================================
    nsz = size(fo)

    if(str1=='x') then
      i = 1
      if (str2 == 'P2C') then
      
        call Prepare_TDMA_interp_RHS_array(str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
            m1rP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), am1x_P2C(:), bm1x_P2C(:), cm1x_P2C(:), dm1x_P2C(:), nsz)
        
      else if (str2 == 'C2P') then

        call Prepare_TDMA_interp_RHS_array(str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
            m1rC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), am1x_C2P(:), bm1x_C2P(:), cm1x_C2P(:), dm1x_C2P(:), nsz)

      else
        call Print_error_msg("108: Error input in prepare_FD_TDMA_RHS.")
      end if

    else if (str1 == 'y') then
      i = 2

      if (str2 == 'P2C') then
        call Prepare_TDMA_interp_RHS_array(str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
            m1rP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), am1y_P2C(:), bm1y_P2C(:), cm1y_P2C(:), dm1y_P2C(:), nsz)

      else if (str2 == 'C2P') then

        call Prepare_TDMA_interp_RHS_array(str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
            m1rC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), am1y_C2P(:), bm1y_C2P(:), cm1y_C2P(:), dm1y_C2P(:), nsz)

      else
        call Print_error_msg("109: Error input in prepare_FD_TDMA_RHS.")
      end if

    else if (str1 == 'z') then

      i = 3

      if (str2 == 'P2C') then
        call Prepare_TDMA_interp_RHS_array(str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
            m1rP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), am1z_P2C(:), bm1z_P2C(:), cm1z_P2C(:), dm1z_P2C(:), nsz)

      else if (str2 == 'C2P') then

        call Prepare_TDMA_interp_RHS_array(str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
            m1rC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), am1z_C2P(:), bm1z_C2P(:), cm1z_C2P(:), dm1z_C2P(:), nsz)

      else
        call Print_error_msg("110: Error input in prepare_FD_TDMA_RHS.")
      end if

    else
      call Print_error_msg("111: No such direction.")
    end if

    return 
  end subroutine Get_midp_interpolation_1D
!===============================================================================
!===============================================================================
!> \brief To caculate the 1st derivative in 1D.
!>
!> This subroutine is called as required to get the 1st derivative
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     str1          string to flag which direction to impletment
!> \param[in]     str2          string to flag C2P or P2C or C2C or P2P
!> \param[in]     d             domain
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!_______________________________________________________________________________
  subroutine Get_1st_derivative_1D(str1, str2, d, fi, fo)
  !===============================================================================
  ! Module files
  !===============================================================================
      use parameters_constant_mod
      use udf_type_mod, only : t_domain
      use tridiagonal_matrix_algorithm
      implicit none
  !===============================================================================
  ! Arguments
  !===============================================================================
      type(t_domain), intent(in) :: d
      character(1), intent(in) :: str1
      character(3), intent(in) :: str2
      real(WP), intent(in) :: fi(:)
      real(WP), intent(out) :: fo(:)
  !===============================================================================
  ! Local arguments
  !===============================================================================
      integer(4) :: i, nsz
  !===============================================================================
  ! Code
  !===============================================================================
      nsz = size(fo)
  
      if(str1=='x') then
        i = 1
  
        if (str2 == 'C2C') then
  
          call Prepare_TDMA_1deri_RHS_array( str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
                d%h1r(i), d1rC2C(:, :, :), fi(:), fo(:) )
          call Solve_TDMA( d%is_periodic(i), fo(:), ad1x_C2C(:), bd1x_C2C(:), cd1x_C2C(:), dd1x_C2C(:), nsz )
  
        else if (str2 == 'P2C') then
          call Prepare_TDMA_1deri_RHS_array(str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
              d%h1r(i), d1rP2C(:, :, :), fi(:), fo(:) )
          call Solve_TDMA(d%is_periodic(i), fo(:), ad1x_P2C(:), bd1x_P2C(:), cd1x_P2C(:), dd1x_P2C(:), nsz)
  
        else if (str2 == 'P2P') then
          
          call Prepare_TDMA_1deri_RHS_array(str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
              d%h1r(i), d1rP2P(:, :, :), fi(:), fo(:) )
              !write(*,'(A,7F8.4)') 'a', ad1x_P2P(:)
              !write(*,'(A,7F8.4)') 'b', bd1x_P2P(:)
              !write(*,'(A,7F8.4)') 'c', cd1x_P2P(:)
              !write(*,'(A,7F8.4)') 'd', dd1x_P2P(:)
              !write(*,'(A,7F8.4)') 'r', fo(:)
          call Solve_TDMA(d%is_periodic(i), fo(:), ad1x_P2P(:), bd1x_P2P(:), cd1x_P2P(:), dd1x_P2P(:), nsz)
          !write(*,'(A,7F8.4)') 'o', fo(:)
        else if (str2 == 'C2P') then
  
          call Prepare_TDMA_1deri_RHS_array(str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
              d%h1r(i), d1rC2P(:, :, :), fi(:), fo(:) )
          call Solve_TDMA(d%is_periodic(i), fo(:), ad1x_C2P(:), bd1x_C2P(:), cd1x_C2P(:), dd1x_C2P(:), nsz)
  
        else
          call Print_error_msg("112: No such staggered scheme defined")
        end if
  
      else if (str1 == 'y') then
        i = 2
  
        if (str2 == 'C2C') then
  
          call Prepare_TDMA_1deri_RHS_array( str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
                d%h1r(i), d1rC2C(:, :, :), fi(:), fo(:) )
          call Solve_TDMA( d%is_periodic(i), fo(:), ad1y_C2C(:), bd1y_C2C(:), cd1y_C2C(:), dd1y_C2C(:), nsz )
          if(d%is_stretching(2)) fo(:) = fo(:) * d%yMappingcc(:, 1)
        
        else if (str2 == 'P2C') then
          call Prepare_TDMA_1deri_RHS_array(str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
              d%h1r(i), d1rP2C(:, :, :), fi(:), fo(:) )
          call Solve_TDMA(d%is_periodic(i), fo(:), ad1y_P2C(:), bd1y_P2C(:), cd1y_P2C(:), dd1y_P2C(:), nsz)
          if(d%is_stretching(2)) fo(:) = fo(:) * d%yMappingcc(:, 1)
  
        else if (str2 == 'P2P') then
  
          call Prepare_TDMA_1deri_RHS_array(str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
              d%h1r(i), d1rP2P(:, :, :), fi(:), fo(:) )
          call Solve_TDMA(d%is_periodic(i), fo(:), ad1y_P2P(:), bd1y_P2P(:), cd1y_P2P(:), dd1y_P2P(:), nsz)
          if(d%is_stretching(2)) fo(:) = fo(:) * d%yMappingpt(1:nsz, 1)
  
        else if (str2 == 'C2P') then
  
          call Prepare_TDMA_1deri_RHS_array(str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
              d%h1r(i), d1rC2P(:, :, :), fi(:), fo(:) )
          call Solve_TDMA(d%is_periodic(i), fo(:), ad1y_C2P(:), bd1y_C2P(:), cd1y_C2P(:), dd1y_C2P(:), nsz)
          if(d%is_stretching(2)) fo(:) = fo(:) * d%yMappingpt(1:nsz, 1)
  
        else
          call Print_error_msg("113: No such staggered scheme defined")
        end if
  
      else if (str1 == 'z') then
  
        i = 3
  
        if (str2 == 'C2C') then
  
          call Prepare_TDMA_1deri_RHS_array( str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
                d%h1r(i), d1rC2C(:, :, :), fi(:), fo(:) )
          call Solve_TDMA( d%is_periodic(i), fo(:), ad1z_C2C(:), bd1z_C2C(:), cd1z_C2C(:), dd1z_C2C(:), nsz )
  
        else if (str2 == 'P2C') then
          call Prepare_TDMA_1deri_RHS_array(str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
              d%h1r(i), d1rP2C(:, :, :), fi(:), fo(:) )
          call Solve_TDMA(d%is_periodic(i), fo(:), ad1z_P2C(:), bd1z_P2C(:), cd1z_P2C(:), dd1z_P2C(:), nsz)
  
        else if (str2 == 'P2P') then
  
          call Prepare_TDMA_1deri_RHS_array(str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
              d%h1r(i), d1rP2P(:, :, :), fi(:), fo(:) )
          call Solve_TDMA(d%is_periodic(i), fo(:), ad1z_P2P(:), bd1z_P2P(:), cd1z_P2P(:), dd1z_P2P(:), nsz)
  
        else if (str2 == 'C2P') then
  
          call Prepare_TDMA_1deri_RHS_array(str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
              d%h1r(i), d1rC2P(:, :, :), fi(:), fo(:) )
          call Solve_TDMA(d%is_periodic(i), fo(:), ad1z_C2P(:), bd1z_C2P(:), cd1z_C2P(:), dd1z_C2P(:), nsz)
  
        else
          call Print_error_msg("114: No such staggered scheme defined")
        end if
  
      else
        call Print_error_msg("115: No such direction defined.")
      end if
  
      return 
    end subroutine Get_1st_derivative_1D


!===============================================================================
!===============================================================================
!> \brief To caculate the 2nd derivative in 1D.
!>
!> This subroutine is called as required to get the 2nd derivative
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     str1          string to flag which direction to impletment
!> \param[in]     str2          string to flag C2P or P2C or C2C or P2P
!> \param[in]     d             domain
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!_______________________________________________________________________________
  subroutine Get_2nd_derivative_1D(str1, str2, d, fi, fo)
  !===============================================================================
  ! Module files
  !===============================================================================
      use parameters_constant_mod
      use udf_type_mod
      use tridiagonal_matrix_algorithm
      implicit none
  !===============================================================================
  ! Arguments
  !===============================================================================
      type(t_domain), intent(in) :: d
      character(1), intent(in) :: str1
      character(3), intent(in) :: str2
      real(WP), intent(in) :: fi(:)
      real(WP), intent(out) :: fo(:)
  !===============================================================================
  ! Local arguments
  !===============================================================================
      integer(4) :: i, nsz
      real(WP), allocatable :: fo1(:)
  !===============================================================================
  ! Code
  !===============================================================================
      nsz = size(fo)
  
      if(str1=='x') then
        i = 1
  
        if (str2 == 'C2C') then
  
          call Prepare_TDMA_2deri_RHS_array(str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
                d%h2r(i), d2rC2C(:, :, :), fi(:), fo(:) )
          call Solve_TDMA(d%is_periodic(i), fo(:), ad2x_C2C(:), bd2x_C2C(:), cd2x_C2C(:), dd2x_C2C(:), nsz )
  
        else if (str2 == 'P2P') then
          
          call Prepare_TDMA_2deri_RHS_array(str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
                d%h2r(i), d2rP2P(:, :, :), fi(:), fo(:) )
          call Solve_TDMA(d%is_periodic(i), fo(:), ad2x_P2P(:), bd2x_P2P(:), cd2x_P2P(:), dd2x_P2P(:), nsz)

        else
          call Print_error_msg("112: No such staggered scheme defined")
        end if
  
      else if (str1 == 'y') then
        i = 2
  
        if (str2 == 'C2C') then
          allocate ( fo1(nsz) ); fo1(:) = ZERO
          ! 2nd, uniform
          call Prepare_TDMA_2deri_RHS_array( str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
                d%h2r(i), d2rC2C(:, :, :), fi(:), fo(:) )
          call Solve_TDMA( d%is_periodic(i), fo(:), ad2y_C2C(:), bd2y_C2C(:), cd2y_C2C(:), dd2y_C2C(:), nsz )
          
          if(d%is_stretching(2)) then
            ! 1st, uniform
            call Prepare_TDMA_1deri_RHS_array( str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
                d%h1r(i), d1rC2C(:, :, :), fi(:), fo1(:) )
            call Solve_TDMA( d%is_periodic(i), fo1(:), ad1y_C2C(:), bd1y_C2C(:), cd1y_C2C(:), dd1y_C2C(:), nsz )
            fo(:) = fo(:) * d%yMappingcc(:, 2) + fo1(:) * d%yMappingcc(:, 3)
          end if
          deallocate(fo1)

        else if (str2 == 'P2P') then
  
          call Prepare_TDMA_2deri_RHS_array(str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
              d%h2r(i), d2rP2P(:, :, :), fi(:), fo(:) )
          call Solve_TDMA(d%is_periodic(i), fo(:), ad2y_P2P(:), bd2y_P2P(:), cd2y_P2P(:), dd2y_P2P(:), nsz )

          if(d%is_stretching(2)) then
            allocate ( fo1(nsz) ); fo1(:) = ZERO
            call Prepare_TDMA_1deri_RHS_array(str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
              d%h1r(i), d1rP2P(:, :, :), fi(:), fo1(:) )
            call Solve_TDMA(d%is_periodic(i), fo1(:), ad1y_P2P(:), bd1y_P2P(:), cd1y_P2P(:), dd1y_P2P(:), nsz )
            fo(:) = fo(:) * d%yMappingpt(1:nsz, 2) + fo1(:) * d%yMappingpt(1:nsz, 3)
            deallocate(fo1)
          end if
  
        else
          call Print_error_msg("113: No such staggered scheme defined")
        end if
  
      else if (str1 == 'z') then
  
        i = 3
  
        if (str2 == 'C2C') then
  
          call Prepare_TDMA_2deri_RHS_array( str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
                d%h2r(i), d2rC2C(:, :, :), fi(:), fo(:) )
          call Solve_TDMA( d%is_periodic(i), fo(:), ad2z_C2C(:), bd2z_C2C(:), cd2z_C2C(:), dd2z_C2C(:), nsz )
  
        else if (str2 == 'P2P') then
  
          call Prepare_TDMA_2deri_RHS_array(str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
              d%h2r(i), d2rP2P(:, :, :), fi(:), fo(:) )
          call Solve_TDMA(d%is_periodic(i), fo(:), ad2z_P2P(:), bd2z_P2P(:), cd2z_P2P(:), dd2z_P2P(:), nsz )
  
        else
          call Print_error_msg("114: No such staggered scheme defined")
        end if
  
      else
        call Print_error_msg("115: No such direction defined.")
      end if
  
      return 
  end subroutine Get_2nd_derivative_1D
!===============================================================================
!===============================================================================
!> \brief To expend to 3D from 1D operation
!>  for interpolation
!-------------------------------------------------------------------------------
!===============================================================================
  subroutine Get_x_midp_C2P_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer(4) :: dim, nox
    integer(4) :: k, j

    dim = 1
    nox = size(fo3d, 1)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Prepare_TDMA_interp_RHS_array('C2P', nox, d%bc(:, dim), &
                d%iNeighb(:, :), m1rC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(dim), fo(:), am1x_C2P(:), bm1x_C2P(:), &
                cm1x_C2P(:), dm1x_C2P(:), nox)
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_midp_C2P_3dArray
!===============================================================================
  subroutine Get_x_midp_P2C_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer(4) :: dim, nox
    integer(4) :: k, j

    dim = 1
    nox = size(fo3d, 1)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Prepare_TDMA_interp_RHS_array('P2C', nox, d%bc(:, dim), &
                d%iNeighb(:, :), m1rP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(dim), fo(:), am1x_P2C(:), bm1x_P2C(:), &
                cm1x_P2C(:), dm1x_P2C(:), nox)
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_midp_P2C_3dArray
!===============================================================================
  subroutine Get_y_midp_C2P_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer(4) :: dim, noy
    integer(4) :: k, i

    dim = 2
    noy = size(fo3d, 2)

    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Prepare_TDMA_interp_RHS_array('C2P', noy, d%bc(:, dim), &
                d%jNeighb(:, :), m1rC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(dim), fo(:), am1y_C2P(:), bm1y_C2P(:), &
                cm1y_C2P(:), dm1y_C2P(:), noy)
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_midp_C2P_3dArray
!===============================================================================
  subroutine Get_y_midp_P2C_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer(4) :: dim, noy
    integer(4) :: k, i

    dim = 2
    noy = size(fo3d, 2)

    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Prepare_TDMA_interp_RHS_array('P2C', noy, d%bc(:, dim), &
                d%jNeighb(:, :), m1rP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(dim), fo(:), am1y_P2C(:), bm1y_P2C(:), &
                cm1y_P2C(:), dm1y_P2C(:), noy)
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_midp_P2C_3dArray
  !===============================================================================
  subroutine Get_z_midp_C2P_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer(4) :: dim, noz
    integer(4) :: j, i

    dim = 3
    noz = size(fo3d, 3)

    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Prepare_TDMA_interp_RHS_array('C2P', noz, d%bc(:, dim), &
                d%kNeighb(:, :), m1rC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(dim), fo(:), am1z_C2P(:), bm1z_C2P(:), &
                cm1z_C2P(:), dm1z_C2P(:), noz)
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_midp_C2P_3dArray
!===============================================================================
  subroutine Get_z_midp_P2C_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer(4) :: dim, noz
    integer(4) :: j, i

    dim = 3
    noz = size(fo3d, 3)

    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Prepare_TDMA_interp_RHS_array('P2C', noz, d%bc(:, dim), &
                d%kNeighb(:, :), m1rP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(dim), fo(:), am1z_P2C(:), bm1z_P2C(:), &
                cm1z_P2C(:), dm1z_P2C(:), noz)
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_midp_P2C_3dArray
!===============================================================================
!> \brief To expend to 3D from 1D operation
!>  for 1st-derivative
!-------------------------------------------------------------------------------
!===============================================================================
  subroutine Get_x_1st_derivative_C2C_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer(4) :: dim, nox
    integer(4) :: k, j

    dim = 1
    nox = size(fo3d, 1)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Prepare_TDMA_1deri_RHS_array( 'C2C', nox, d%bc(:, dim), &
                d%iNeighb(:, :), d%h1r(dim), d1rC2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(dim), fo(:), ad1x_C2C(:), bd1x_C2C(:), &
                cd1x_C2C(:), dd1x_C2C(:), nox )
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_1st_derivative_C2C_3dArray

!===============================================================================
  subroutine Get_x_1st_derivative_P2P_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer(4) :: dim, nox
    integer(4) :: k, j

    dim = 1
    nox = size(fo3d, 1)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Prepare_TDMA_1deri_RHS_array( 'P2P', nox, d%bc(:, dim), &
                d%iNeighb(:, :), d%h1r(dim), d1rP2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(dim), fo(:), ad1x_P2P(:), bd1x_P2P(:), &
                cd1x_P2P(:), dd1x_P2P(:), nox )
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_1st_derivative_P2P_3dArray

!===============================================================================
  subroutine Get_x_1st_derivative_C2P_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer(4) :: dim, nox
    integer(4) :: k, j

    dim = 1
    nox = size(fo3d, 1)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Prepare_TDMA_1deri_RHS_array( 'C2P', nox, d%bc(:, dim), &
                d%iNeighb(:, :), d%h1r(dim), d1rC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(dim), fo(:), ad1x_C2P(:), bd1x_C2P(:), &
                cd1x_C2P(:), dd1x_C2P(:), nox )
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_1st_derivative_C2P_3dArray

!===============================================================================
  subroutine Get_x_1st_derivative_P2C_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer(4) :: dim, nox
    integer(4) :: k, j

    dim = 1
    nox = size(fo3d, 1)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Prepare_TDMA_1deri_RHS_array( 'P2C', nox, d%bc(:, dim), &
                d%iNeighb(:, :), d%h1r(dim), d1rP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(dim), fo(:), ad1x_P2C(:), bd1x_P2C(:), &
                cd1x_P2C(:), dd1x_P2C(:), nox )
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_1st_derivative_P2C_3dArray

!===============================================================================
  subroutine Get_y_1st_derivative_C2C_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer(4) :: dim, noy
    integer(4) :: k, i

    dim = 2
    noy = size(fi3d, 2)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Prepare_TDMA_1deri_RHS_array( 'C2C', noy, d%bc(:, dim), &
                d%jNeighb(:, :), d%h1r(dim), d1rC2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(dim), fo(:), ad1y_C2C(:), bd1y_C2C(:), &
                cd1y_C2C(:), dd1y_C2C(:), noy )
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_1st_derivative_C2C_3dArray

!===============================================================================
  subroutine Get_y_1st_derivative_P2P_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer(4) :: dim, noy
    integer(4) :: k, i

    dim = 2
    noy = size(fi3d, 2)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Prepare_TDMA_1deri_RHS_array( 'P2P', noy, d%bc(:, dim), &
                d%jNeighb(:, :), d%h1r(dim), d1rP2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(dim), fo(:), ad1y_P2P(:), bd1y_P2P(:), &
                cd1y_P2P(:), dd1y_P2P(:), noy )
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_1st_derivative_P2P_3dArray

!===============================================================================
  subroutine Get_y_1st_derivative_C2P_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer(4) :: dim, noy
    integer(4) :: k, i

    dim = 2
    noy = size(fi3d, 2)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Prepare_TDMA_1deri_RHS_array( 'C2P', noy, d%bc(:, dim), &
                d%jNeighb(:, :), d%h1r(dim), d1rC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(dim), fo(:), ad1y_C2P(:), bd1y_C2P(:), &
                cd1y_C2P(:), dd1y_C2P(:), noy )
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_1st_derivative_C2P_3dArray

!===============================================================================
  subroutine Get_y_1st_derivative_P2C_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer(4) :: dim, noy
    integer(4) :: k, i

    dim = 2
    noy = size(fi3d, 2)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Prepare_TDMA_1deri_RHS_array( 'P2C', noy, d%bc(:, dim), &
                d%jNeighb(:, :), d%h1r(dim), d1rP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(dim), fo(:), ad1y_P2C(:), bd1y_P2C(:), &
                cd1y_P2C(:), dd1y_P2C(:), noy )
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_1st_derivative_P2C_3dArray

!===============================================================================
  subroutine Get_z_1st_derivative_C2C_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer(4) :: dim, noz
    integer(4) :: j, i

    dim = 3
    noz = size(fo3d, 3)

    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Prepare_TDMA_1deri_RHS_array( 'C2C', noz, d%bc(:, dim), &
                d%kNeighb(:, :), d%h1r(dim), d1rC2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(dim), fo(:), ad1z_C2C(:), bd1z_C2C(:), &
                cd1z_C2C(:), dd1z_C2C(:), noz )
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_1st_derivative_C2C_3dArray

!===============================================================================
  subroutine Get_z_1st_derivative_P2P_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer(4) :: dim, noz
    integer(4) :: j, i

    dim = 3
    noz = size(fo3d, 3)

    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Prepare_TDMA_1deri_RHS_array( 'P2P', noz, d%bc(:, dim), &
                d%kNeighb(:, :), d%h1r(dim), d1rP2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(dim), fo(:), ad1z_P2P(:), bd1z_P2P(:), &
                cd1z_P2P(:), dd1z_P2P(:), noz )
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_1st_derivative_P2P_3dArray

!===============================================================================
  subroutine Get_z_1st_derivative_C2P_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer(4) :: dim, noz
    integer(4) :: j, i

    dim = 3
    noz = size(fo3d, 3)

    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Prepare_TDMA_1deri_RHS_array( 'C2P', noz, d%bc(:, dim), &
                d%kNeighb(:, :), d%h1r(dim), d1rC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(dim), fo(:), ad1z_C2P(:), bd1z_C2P(:), &
                cd1z_C2P(:), dd1z_C2P(:), noz )
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_1st_derivative_C2P_3dArray

  !===============================================================================
  subroutine Get_z_1st_derivative_P2C_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer(4) :: dim, noz
    integer(4) :: j, i

    dim = 3
    noz = size(fo3d, 3)

    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Prepare_TDMA_1deri_RHS_array( 'P2C', noz, d%bc(:, dim), &
                d%kNeighb(:, :), d%h1r(dim), d1rP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(dim), fo(:), ad1z_P2C(:), bd1z_P2C(:), &
                cd1z_P2C(:), dd1z_P2C(:), noz )
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_1st_derivative_P2C_3dArray
!===============================================================================
!> \brief To expend to 3D from 1D operation
!>  for 2nd-derivative
!-------------------------------------------------------------------------------
!===============================================================================
  subroutine Get_x_2nd_derivative_C2C_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer(4) :: dim, nox
    integer(4) :: k, j

    dim = 1
    nox = size(fo3d, 1)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Prepare_TDMA_2deri_RHS_array( 'C2C', nox, d%bc(:, dim), &
                d%iNeighb(:, :), d%h2r(dim), d2rC2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(dim), fo(:), ad2x_C2C(:), bd2x_C2C(:), &
                cd2x_C2C(:), dd2x_C2C(:), nox )
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return
  end subroutine Get_x_2nd_derivative_C2C_3dArray

  !===============================================================================
  subroutine Get_x_2nd_derivative_P2P_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer(4) :: dim, nox
    integer(4) :: k, j

    dim = 1
    nox = size(fo3d, 1)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        call Prepare_TDMA_2deri_RHS_array( 'P2P', nox, d%bc(:, dim), &
                d%iNeighb(:, :), d%h2r(dim), d2rP2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(dim), fo(:), ad2x_P2P(:), bd2x_P2P(:), &
                cd2x_P2P(:), dd2x_P2P(:), nox )
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return
  end subroutine Get_x_2nd_derivative_P2P_3dArray
!===============================================================================
  subroutine Get_y_2nd_derivative_C2C_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer(4) :: dim, noy
    integer(4) :: k, i

    dim = 2
    noy = size(fo3d, 2)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Prepare_TDMA_2deri_RHS_array( 'C2C', noy, d%bc(:, dim), &
                d%jNeighb(:, :), d%h2r(dim), d2rC2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(dim), fo(:), ad2y_C2C(:), bd2y_C2C(:), &
                cd2y_C2C(:), dd2y_C2C(:), noy )
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return
  end subroutine Get_y_2nd_derivative_C2C_3dArray

!===============================================================================
  subroutine Get_y_2nd_derivative_P2P_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer(4) :: dim, noy
    integer(4) :: k, i

    dim = 2
    noy = size(fo3d, 2)
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        call Prepare_TDMA_2deri_RHS_array( 'P2P', noy, d%bc(:, dim), &
                d%jNeighb(:, :), d%h2r(dim), d2rP2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(dim), fo(:), ad2y_P2P(:), bd2y_P2P(:), &
                cd2y_P2P(:), dd2y_P2P(:), noy )
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return
  end subroutine Get_y_2nd_derivative_P2P_3dArray
!===============================================================================
  subroutine Get_z_2nd_derivative_C2C_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in)   :: d
    real(WP),       intent(in)   :: fi3d(:, :, :)
    real(WP),       intent(out)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer(4) :: dim, noz
    integer(4) :: j, i

    dim = 3
    noz = size(fo3d, 3)
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Prepare_TDMA_2deri_RHS_array( 'C2C', noz, d%bc(:, dim), &
                d%kNeighb(:, :), d%h2r(dim), d2rC2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(dim), fo(:), ad2z_C2C(:), bd2z_C2C(:), &
                cd2z_C2C(:), dd2z_C2C(:), noz )
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return
  end subroutine Get_z_2nd_derivative_C2C_3dArray

!===============================================================================
  subroutine Get_z_2nd_derivative_P2P_3dArray(d, fi3d, fo3d)
    use parameters_constant_mod, only : ZERO
    use udf_type_mod, only : t_domain
    use tridiagonal_matrix_algorithm, only : Solve_TDMA
    implicit none

    type(t_domain), intent(in   )  :: d
    real(WP),       intent(in   )  :: fi3d(:, :, :)
    real(WP),       intent(inout)  :: fo3d(:, :, :)

    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer(4) :: dim, noz
    integer(4) :: j, i

    dim = 3
    noz = size(fo3d, dim)
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        call Prepare_TDMA_2deri_RHS_array( 'P2P', noz, d%bc(:, dim), &
                d%kNeighb(:, :), d%h2r(dim), d2rP2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(dim), fo(:), ad2z_P2P(:), bd2z_P2P(:), &
                cd2z_P2P(:), dd2z_P2P(:), noz )
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return
  end subroutine Get_z_2nd_derivative_P2P_3dArray

!===============================================================================
!>\brief : to calculate:
!>         fo = \int_1^nx \int_
!===============================================================================
  subroutine Get_volumetric_average_3d(d, fi3d, fo, is_stored_nyp)
    ! how to get a high order bulk value?
    use parameters_constant_mod, only : ZERO, HALF
    use udf_type_mod,            only : t_domain
    implicit none
  
    type(t_domain), intent(in)  :: d
    real(WP),       intent(in)  :: fi3d(:, :, :)
    real(WP),       intent(out) :: fo
    logical,        intent(in)  :: is_stored_nyp
  
    real(WP), allocatable   :: fi3dy(:, :, :)
    real(WP)   :: vol
    integer(4) :: i, j, k
    integer(4) :: nix, niy, niz

    fo = ZERO
    vol = ZERO
    nix = size(fi3d, 1)
    niz = size(fi3d, 3)
    if(.not. is_stored_nyp) then
!-------------------------------------------------------------------------------
!   if variable is not stored in y-nodes, extends them to y-nodes
!-------------------------------------------------------------------------------
      if( d%is_periodic(2) ) then
        niy = size(fi3d, 2)
      else
        niy = size(fi3d, 2) + 1
      end if
      allocate( fi3dy(nix, niy, niz) )
      call Get_y_midp_C2P_3dArray ( d, fi3d, fi3dy)
      do k = 1, niz
        do i = 1, nix
          do j = 1, niy
            fo = fo + &
                ( fi3dy(i, d%jNeighb(3, j), k) + fi3d(i, j, k) ) * &
                ( d%yp(j + 1) - d%yc(j) ) * HALF + &
                ( fi3dy(i, j,               k) + fi3d(i, j, k) ) * &
                ( d%yc(j    ) - d%yp(j) ) * HALF
            vol = vol + ( d%yp(j + 1) - d%yp(j) )
          end do
        end do
      end do
      deallocate(fi3dy)
    else
!-------------------------------------------------------------------------------
!   if variable is stored in y-nodes, extend them to y-cell centres
!-------------------------------------------------------------------------------
      if( d%is_periodic(2) ) then
        niy = size(fi3d, 2)
      else
        niy = size(fi3d, 2) - 1
      end if
      allocate( fi3dy(nix, niy, niz) )
      call Get_y_midp_P2C_3dArray ( d, fi3d, fi3dy)
      do k = 1, niz
        do i = 1, nix
          do j = 1, niy
            fo = fo + &
                ( fi3d(i, d%jNeighb(3, j), k) + fi3dy(i, j, k) ) * &
                ( d%yp(j + 1) - d%yc(j) ) * HALF + &
                ( fi3d(i, j,               k) + fi3dy(i, j, k) ) * &
                ( d%yc(j    ) - d%yp(j) ) * HALF
            vol = vol + ( d%yp(j + 1) - d%yp(j) )
          end do
        end do
      end do
      deallocate(fi3dy)
    end if

    fo = fo / real(nix * niy, WP) / vol
    return 
  end subroutine Get_volumetric_average_3d

end module