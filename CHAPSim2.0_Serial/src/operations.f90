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
  use precision_mod
  implicit none
  !----------------------------------------------------------------
  ! coefficients for TDMA of 1st deriviative  
  ! to store coefficients for TDMA
  ! eg, cfpC2C(5, 3, 4)
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
  ! collocated C2C
  real(WP) :: cfpC2C(5, 3, 4)
  real(WP) :: cfrC2C(5, 3, 4)
  
  ! collocated P2P
  real(WP) :: cfpP2P(5, 3, 4)
  real(WP) :: cfrP2P(5, 3, 4)

  ! staggered C2P
  real(WP) :: cfpC2P(5, 3, 4)
  real(WP) :: cfrC2P(5, 3, 4)

  ! staggered P2C
  real(WP) :: cfpP2C(5, 3, 4)
  real(WP) :: cfrP2C(5, 3, 4)

  ! interpolation P2C
  real(WP) :: mfpP2C(5, 3, 4)
  real(WP) :: mfrP2C(5, 3, 4)

  ! interpolation C2P
  real(WP) :: mfpC2P(5, 3, 4)
  real(WP) :: mfrC2P(5, 3, 4)


  ! TDMA 1st LHS Matrix
  ! x
  real(WP), allocatable :: adx_P2P(:)
  real(WP), allocatable :: bdx_P2P(:)
  real(WP), allocatable :: cdx_P2P(:)
  real(WP), allocatable :: ddx_P2P(:)

  real(WP), allocatable :: adx_C2C(:)
  real(WP), allocatable :: bdx_C2C(:)
  real(WP), allocatable :: cdx_C2C(:)
  real(WP), allocatable :: ddx_C2C(:)

  real(WP), allocatable :: adx_P2C(:)
  real(WP), allocatable :: bdx_P2C(:)
  real(WP), allocatable :: cdx_P2C(:)
  real(WP), allocatable :: ddx_P2C(:)

  real(WP), allocatable :: adx_C2P(:)
  real(WP), allocatable :: bdx_C2P(:)
  real(WP), allocatable :: cdx_C2P(:)
  real(WP), allocatable :: ddx_C2P(:)

  real(WP), allocatable :: amx_P2C(:)
  real(WP), allocatable :: bmx_P2C(:)
  real(WP), allocatable :: cmx_P2C(:)
  real(WP), allocatable :: dmx_P2C(:)

  real(WP), allocatable :: amx_C2P(:)
  real(WP), allocatable :: bmx_C2P(:)
  real(WP), allocatable :: cmx_C2P(:)
  real(WP), allocatable :: dmx_C2P(:)

  ! y
  real(WP), allocatable :: ady_P2P(:)
  real(WP), allocatable :: bdy_P2P(:)
  real(WP), allocatable :: cdy_P2P(:)
  real(WP), allocatable :: ddy_P2P(:)

  real(WP), allocatable :: ady_C2C(:)
  real(WP), allocatable :: bdy_C2C(:)
  real(WP), allocatable :: cdy_C2C(:)
  real(WP), allocatable :: ddy_C2C(:)

  real(WP), allocatable :: ady_P2C(:)
  real(WP), allocatable :: bdy_P2C(:)
  real(WP), allocatable :: cdy_P2C(:)
  real(WP), allocatable :: ddy_P2C(:)

  real(WP), allocatable :: ady_C2P(:)
  real(WP), allocatable :: bdy_C2P(:)
  real(WP), allocatable :: cdy_C2P(:)
  real(WP), allocatable :: ddy_C2P(:)

  real(WP), allocatable :: amy_P2C(:)
  real(WP), allocatable :: bmy_P2C(:)
  real(WP), allocatable :: cmy_P2C(:)
  real(WP), allocatable :: dmy_P2C(:)

  real(WP), allocatable :: amy_C2P(:)
  real(WP), allocatable :: bmy_C2P(:)
  real(WP), allocatable :: cmy_C2P(:)
  real(WP), allocatable :: dmy_C2P(:)


  ! z
  real(WP), allocatable :: adz_P2P(:)
  real(WP), allocatable :: bdz_P2P(:)
  real(WP), allocatable :: cdz_P2P(:)
  real(WP), allocatable :: ddz_P2P(:)

  real(WP), allocatable :: adz_C2C(:)
  real(WP), allocatable :: bdz_C2C(:)
  real(WP), allocatable :: cdz_C2C(:)
  real(WP), allocatable :: ddz_C2C(:)

  real(WP), allocatable :: adz_P2C(:)
  real(WP), allocatable :: bdz_P2C(:)
  real(WP), allocatable :: cdz_P2C(:)
  real(WP), allocatable :: ddz_P2C(:)

  real(WP), allocatable :: adz_C2P(:)
  real(WP), allocatable :: bdz_C2P(:)
  real(WP), allocatable :: cdz_C2P(:)
  real(WP), allocatable :: ddz_C2P(:)

  real(WP), allocatable :: amz_P2C(:)
  real(WP), allocatable :: bmz_P2C(:)
  real(WP), allocatable :: cmz_P2C(:)
  real(WP), allocatable :: dmz_P2C(:)

  real(WP), allocatable :: amz_C2P(:)
  real(WP), allocatable :: bmz_C2P(:)
  real(WP), allocatable :: cmz_C2P(:)
  real(WP), allocatable :: dmz_C2P(:)

  private :: Assign_TDMA_coeffs
  private :: Buildup_TDMA_LHS_array
  private :: Prepare_TDMA_interp_RHS_array
  private :: Prepare_TDMA_deri_RHS_array

  public  :: Prepare_coeffs_for_operations
  public  :: Get_midp_interpolation
  public  :: Get_1st_derivative
  
  public  :: Test_interpolation

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

    real(WP) :: alpha, a, b, c
    real(WP) :: alpha1, a1, b1, c1
    real(WP) :: alpha2, a2, b2, c2
!______________________________________________________________________________!
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
    cfpC2C(1:5, 1, IBC_PERIODIC) = alpha
    cfpC2C(1:5, 2, IBC_PERIODIC) = ONE
    cfpC2C(1:5, 3, IBC_PERIODIC) = alpha
    cfrC2C(1:5, 1, IBC_PERIODIC) = a / TWO ! a/2
    cfrC2C(1:5, 2, IBC_PERIODIC) = b / FOUR ! b/4
    cfrC2C(1:5, 3, IBC_PERIODIC) = c ! not used

    !P2P for periodic b.c.
    cfpP2P(:, :, IBC_PERIODIC) = cfpC2C(:, :, IBC_PERIODIC)
    cfrP2P(:, :, IBC_PERIODIC) = cfrC2C(:, :, IBC_PERIODIC)

    !C2C for symmetric b.c.
    cfpC2C(1, 1, IBC_SYMMETRIC) = ZERO ! not used
    cfpC2C(1, 2, IBC_SYMMETRIC) = ONE - alpha
    cfpC2C(1, 3, IBC_SYMMETRIC) = alpha

    cfpC2C(2:4, 1, IBC_SYMMETRIC) = alpha
    cfpC2C(2:4, 2, IBC_SYMMETRIC) = ONE
    cfpC2C(2:4, 3, IBC_SYMMETRIC) = alpha

    cfpC2C(5, 1, IBC_SYMMETRIC) = alpha
    cfpC2C(5, 2, IBC_SYMMETRIC) = ONE - alpha
    cfpC2C(5, 3, IBC_SYMMETRIC) = ZERO ! not used

    cfrC2C(1:5, 1, IBC_SYMMETRIC) = a / TWO ! a/2
    cfrC2C(1:5, 2, IBC_SYMMETRIC) = b / FOUR ! b/4
    cfrC2C(1:5, 3, IBC_SYMMETRIC) = c ! not used

    !P2P for symmetric b.c.
    cfpP2P(1, 1, IBC_SYMMETRIC) = ZERO ! not used
    cfpP2P(1, 2, IBC_SYMMETRIC) = ONE
    cfpP2P(1, 3, IBC_SYMMETRIC) = ZERO

    cfpP2P(2:4, :, IBC_SYMMETRIC) = cfpC2C(2:4, :, IBC_SYMMETRIC)

    cfpP2P(5, 1, IBC_SYMMETRIC) = ZERO
    cfpP2P(5, 2, IBC_SYMMETRIC) = ONE
    cfpP2P(5, 3, IBC_SYMMETRIC) = ZERO ! not used

    cfrP2P(:, :, IBC_SYMMETRIC) = cfrC2C(:, :, IBC_SYMMETRIC)

    !C2C for asymmetric b.c.
    cfpC2C(:, :, IBC_ASYMMETRIC) = cfpC2C(:, :, IBC_SYMMETRIC)
    cfrC2C(:, :, IBC_ASYMMETRIC) = cfrC2C(:, :, IBC_SYMMETRIC)
    !P2P for asymmetric b.c.
    cfpP2P(:, :, IBC_ASYMMETRIC) = cfpP2P(:, :, IBC_SYMMETRIC)
    cfrP2P(:, :, IBC_ASYMMETRIC) = cfrP2P(:, :, IBC_SYMMETRIC)

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

    cfpC2C(1, 1, IBC_UDIRICHLET) = alpha1 ! not used
    cfpC2C(1, 2, IBC_UDIRICHLET) = ONE
    cfpC2C(1, 3, IBC_UDIRICHLET) = alpha1
    cfrC2C(1, 1, IBC_UDIRICHLET) = a1
    cfrC2C(1, 2, IBC_UDIRICHLET) = b1
    cfrC2C(1, 3, IBC_UDIRICHLET) = c1

    cfpC2C(2, 1, IBC_UDIRICHLET) = alpha2
    cfpC2C(2, 2, IBC_UDIRICHLET) = ONE
    cfpC2C(2, 3, IBC_UDIRICHLET) = alpha2
    cfrC2C(2, 1, IBC_UDIRICHLET) = a2 / TWO
    cfrC2C(2, 2, IBC_UDIRICHLET) = b2 / FOUR ! not used
    cfrC2C(2, 3, IBC_UDIRICHLET) = c2 ! not used

    cfpC2C(3, 1, IBC_UDIRICHLET) = alpha
    cfpC2C(3, 2, IBC_UDIRICHLET) = ONE
    cfpC2C(3, 3, IBC_UDIRICHLET) = alpha
    cfrC2C(3, 1, IBC_UDIRICHLET) = a / TWO ! a/2
    cfrC2C(3, 2, IBC_UDIRICHLET) = b / FOUR ! b/4
    cfrC2C(3, 3, IBC_UDIRICHLET) = c ! not used

    cfpC2C(4, 1, IBC_UDIRICHLET) = alpha2
    cfpC2C(4, 2, IBC_UDIRICHLET) = ONE
    cfpC2C(4, 3, IBC_UDIRICHLET) = alpha2
    cfrC2C(4, 1, IBC_UDIRICHLET) = a2 / TWO
    cfrC2C(4, 2, IBC_UDIRICHLET) = b2 / FOUR ! not used
    cfrC2C(4, 3, IBC_UDIRICHLET) = c2 ! not used

    cfpC2C(5, 1, IBC_UDIRICHLET) = alpha1
    cfpC2C(5, 2, IBC_UDIRICHLET) = ONE
    cfpC2C(5, 3, IBC_UDIRICHLET) = alpha1 ! not used
    cfrC2C(5, 1, IBC_UDIRICHLET) = -a1
    cfrC2C(5, 2, IBC_UDIRICHLET) = -b1
    cfrC2C(5, 3, IBC_UDIRICHLET) = -c1

    cfpP2P(:, :, IBC_UDIRICHLET) = cfpC2C(:, :, IBC_UDIRICHLET)
    cfrP2P(:, :, IBC_UDIRICHLET) = cfrC2C(:, :, IBC_UDIRICHLET)
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
    cfpC2P(1:5, 1, IBC_PERIODIC) = alpha
    cfpC2P(1:5, 2, IBC_PERIODIC) = ONE
    cfpC2P(1:5, 3, IBC_PERIODIC) = alpha
    cfrC2P(1:5, 1, IBC_PERIODIC) = a ! a
    cfrC2P(1:5, 2, IBC_PERIODIC) = b / THREE ! b/3
    cfrC2P(1:5, 3, IBC_PERIODIC) = c ! not used

    !P2C for periodic b.c.
    cfpP2C(:, :, IBC_PERIODIC) = cfpC2P(:, :, IBC_PERIODIC)
    cfrP2C(:, :, IBC_PERIODIC) = cfrC2P(:, :, IBC_PERIODIC)

    !C2P for symmetric 
    cfpC2P(:, :, IBC_SYMMETRIC) = cfpP2P(:, :, IBC_SYMMETRIC)
    cfrC2P(:, :, IBC_SYMMETRIC) = cfrC2P(:, :, IBC_PERIODIC)

    !P2C for symmetric 
    cfpP2C(:, :, IBC_SYMMETRIC) = cfpC2C(:, :, IBC_SYMMETRIC)
    cfrP2C(:, :, IBC_SYMMETRIC) = cfrC2P(:, :, IBC_SYMMETRIC)

    !C2P for asymmetric 
    cfpC2P(:, :, IBC_ASYMMETRIC) = cfpC2P(:, :, IBC_SYMMETRIC)
    cfrC2P(:, :, IBC_ASYMMETRIC) = cfrC2P(:, :, IBC_SYMMETRIC)

    !P2C for asymmetric 
    cfpP2C(:, :, IBC_ASYMMETRIC) = cfpP2C(:, :, IBC_SYMMETRIC)
    cfrP2C(:, :, IBC_ASYMMETRIC) = cfrP2C(:, :, IBC_SYMMETRIC)

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

    cfpP2C(1, 1, IBC_UDIRICHLET) = alpha1 ! not used
    cfpP2C(1, 2, IBC_UDIRICHLET) = ONE
    cfpP2C(1, 3, IBC_UDIRICHLET) = alpha1
    cfrP2C(1, 1, IBC_UDIRICHLET) = a1
    cfrP2C(1, 2, IBC_UDIRICHLET) = b1
    cfrP2C(1, 3, IBC_UDIRICHLET) = c1

    cfpP2C(2, 1, IBC_UDIRICHLET) = alpha2
    cfpP2C(2, 2, IBC_UDIRICHLET) = ONE
    cfpP2C(2, 3, IBC_UDIRICHLET) = alpha2
    cfrP2C(2, 1, IBC_UDIRICHLET) = a2
    cfrP2C(2, 2, IBC_UDIRICHLET) = b2 / THREE ! not used
    cfrP2C(2, 3, IBC_UDIRICHLET) = c2 ! not used

    cfpP2C(3, 1, IBC_UDIRICHLET) = alpha
    cfpP2C(3, 2, IBC_UDIRICHLET) = ONE
    cfpP2C(3, 3, IBC_UDIRICHLET) = alpha
    cfrP2C(3, 1, IBC_UDIRICHLET) = a
    cfrP2C(3, 2, IBC_UDIRICHLET) = b / THREE
    cfrP2C(3, 3, IBC_UDIRICHLET) = c ! not used

    cfpP2C(4, 1, IBC_UDIRICHLET) = alpha2
    cfpP2C(4, 2, IBC_UDIRICHLET) = ONE
    cfpP2C(4, 3, IBC_UDIRICHLET) = alpha2
    cfrP2C(4, 1, IBC_UDIRICHLET) = a2
    cfrP2C(4, 2, IBC_UDIRICHLET) = b2 / THREE ! not used
    cfrP2C(4, 3, IBC_UDIRICHLET) = c2 ! not used

    cfpP2C(5, 1, IBC_UDIRICHLET) = alpha1
    cfpP2C(5, 2, IBC_UDIRICHLET) = ONE
    cfpP2C(5, 3, IBC_UDIRICHLET) = alpha1 ! not used
    cfrP2C(5, 1, IBC_UDIRICHLET) = -a1
    cfrP2C(5, 2, IBC_UDIRICHLET) = -b1
    cfrP2C(5, 3, IBC_UDIRICHLET) = -c1

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

    cfpC2P(1, 1, IBC_UDIRICHLET) = alpha1 ! not used
    cfpC2P(1, 2, IBC_UDIRICHLET) = ONE
    cfpC2P(1, 3, IBC_UDIRICHLET) = alpha1
    cfrC2P(1, 1, IBC_UDIRICHLET) = a1
    cfrC2P(1, 2, IBC_UDIRICHLET) = b1
    cfrC2P(1, 3, IBC_UDIRICHLET) = c1

    cfpC2P(2, 1, IBC_UDIRICHLET) = alpha2
    cfpC2P(2, 2, IBC_UDIRICHLET) = ONE
    cfpC2P(2, 3, IBC_UDIRICHLET) = alpha2
    cfrC2P(2, 1, IBC_UDIRICHLET) = a2
    cfrC2P(2, 2, IBC_UDIRICHLET) = b2 / THREE ! not used
    cfrC2P(2, 3, IBC_UDIRICHLET) = c2 ! not used

    cfpC2P(3, 1, IBC_UDIRICHLET) = alpha
    cfpC2P(3, 2, IBC_UDIRICHLET) = ONE
    cfpC2P(3, 3, IBC_UDIRICHLET) = alpha
    cfrC2P(3, 1, IBC_UDIRICHLET) = a
    cfrC2P(3, 2, IBC_UDIRICHLET) = b / THREE
    cfrC2P(3, 3, IBC_UDIRICHLET) = c ! not used

    cfpC2P(4, 1, IBC_UDIRICHLET) = alpha2
    cfpC2P(4, 2, IBC_UDIRICHLET) = ONE
    cfpC2P(4, 3, IBC_UDIRICHLET) = alpha2
    cfrC2P(4, 1, IBC_UDIRICHLET) = a2
    cfrC2P(4, 2, IBC_UDIRICHLET) = b2 / THREE ! not used
    cfrC2P(4, 3, IBC_UDIRICHLET) = c2 ! not used

    cfpC2P(5, 1, IBC_UDIRICHLET) = alpha1
    cfpC2P(5, 2, IBC_UDIRICHLET) = ONE
    cfpC2P(5, 3, IBC_UDIRICHLET) = alpha1 ! not used
    cfrC2P(5, 1, IBC_UDIRICHLET) = -a1
    cfrC2P(5, 2, IBC_UDIRICHLET) = -b1
    cfrC2P(5, 3, IBC_UDIRICHLET) = -c1

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
    mfpC2P(1:5, 1, IBC_PERIODIC) = alpha
    mfpC2P(1:5, 2, IBC_PERIODIC) = ONE
    mfpC2P(1:5, 3, IBC_PERIODIC) = alpha
    mfrC2P(1:5, 1, IBC_PERIODIC) = a / TWO
    mfrC2P(1:5, 2, IBC_PERIODIC) = b / TWO
    mfrC2P(1:5, 3, IBC_PERIODIC) = c ! not used

    !P2C for periodic b.c.
    mfpP2C(:, :, IBC_PERIODIC) = mfpC2P(:, :, IBC_PERIODIC)
    mfrP2C(:, :, IBC_PERIODIC) = mfrC2P(:, :, IBC_PERIODIC)

    !C2P for symmetric, orthogonal, eg. u in y direction.
    mfpC2P(1, 1, IBC_SYMMETRIC) = ZERO ! not used
    mfpC2P(1, 2, IBC_SYMMETRIC) = ONE
    mfpC2P(1, 3, IBC_SYMMETRIC) = alpha + alpha

    mfpC2P(2:4, 1, IBC_SYMMETRIC) = alpha
    mfpC2P(2:4, 2, IBC_SYMMETRIC) = ONE
    mfpC2P(2:4, 3, IBC_SYMMETRIC) = alpha

    mfpC2P(5, 1, IBC_SYMMETRIC) = alpha + alpha
    mfpC2P(5, 2, IBC_SYMMETRIC) = ONE
    mfpC2P(5, 3, IBC_SYMMETRIC) = ZERO ! not used.

    mfrC2P(1:5, 1, IBC_SYMMETRIC) = a / TWO
    mfrC2P(1:5, 2, IBC_SYMMETRIC) = b / TWO
    mfrC2P(1:5, 3, IBC_SYMMETRIC) = ZERO ! not used. 

    !C2P for symmetric, parallel, eg. v in y direction.
    mfpC2P(1, 1, IBC_ASYMMETRIC) = ZERO ! not used
    mfpC2P(1, 2, IBC_ASYMMETRIC) = ONE
    mfpC2P(1, 3, IBC_ASYMMETRIC) = alpha - alpha

    mfpC2P(2:4, 1, IBC_ASYMMETRIC) = alpha
    mfpC2P(2:4, 2, IBC_ASYMMETRIC) = ONE
    mfpC2P(2:4, 3, IBC_ASYMMETRIC) = alpha

    mfpC2P(5, 1, IBC_ASYMMETRIC) = alpha - alpha
    mfpC2P(5, 2, IBC_ASYMMETRIC) = ONE
    mfpC2P(5, 3, IBC_ASYMMETRIC) = ZERO ! not used.

    mfrC2P(1:5, 1, IBC_ASYMMETRIC) = a / TWO
    mfrC2P(1:5, 2, IBC_ASYMMETRIC) = b / TWO
    mfrC2P(1:5, 3, IBC_ASYMMETRIC) = ZERO ! not used. 

    !P2C for symmetric, orthogonal, eg. u in y direction.
    mfpP2C(1, 1, IBC_SYMMETRIC) = ZERO ! not used
    mfpP2C(1, 2, IBC_SYMMETRIC) = ONE + alpha
    mfpP2C(1, 3, IBC_SYMMETRIC) = alpha

    mfpP2C(2:4, 1, IBC_SYMMETRIC) = alpha
    mfpP2C(2:4, 2, IBC_SYMMETRIC) = ONE
    mfpP2C(2:4, 3, IBC_SYMMETRIC) = alpha

    mfpP2C(5, 1, IBC_SYMMETRIC) = alpha
    mfpP2C(5, 2, IBC_SYMMETRIC) = ONE + alpha
    mfpP2C(5, 3, IBC_SYMMETRIC) = ZERO ! not used.

    mfrP2C(1:5, 1, IBC_SYMMETRIC) = a / TWO
    mfrP2C(1:5, 2, IBC_SYMMETRIC) = b / TWO
    mfrP2C(1:5, 3, IBC_SYMMETRIC) = ZERO ! not used. 

    !P2C for symmetric, parallel, eg. v in y direction.
    mfpP2C(1, 1, IBC_ASYMMETRIC) = ZERO ! not used
    mfpP2C(1, 2, IBC_ASYMMETRIC) = ONE - alpha
    mfpP2C(1, 3, IBC_ASYMMETRIC) = alpha

    mfpP2C(2:4, 1, IBC_ASYMMETRIC) = alpha
    mfpP2C(2:4, 2, IBC_ASYMMETRIC) = ONE
    mfpP2C(2:4, 3, IBC_ASYMMETRIC) = alpha

    mfpP2C(5, 1, IBC_ASYMMETRIC) = alpha
    mfpP2C(5, 2, IBC_ASYMMETRIC) = ONE - alpha
    mfpP2C(5, 3, IBC_ASYMMETRIC) = ZERO ! not used.

    mfrP2C(1:5, 1, IBC_ASYMMETRIC) = a / TWO
    mfrP2C(1:5, 2, IBC_ASYMMETRIC) = b / TWO
    mfrP2C(1:5, 3, IBC_ASYMMETRIC) = ZERO ! not used. 

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
    mfpP2C(1, 1, IBC_UDIRICHLET) = ZERO ! not used
    mfpP2C(1, 2, IBC_UDIRICHLET) = ONE
    mfpP2C(1, 3, IBC_UDIRICHLET) = alpha1
    mfrP2C(1, 1, IBC_UDIRICHLET) = a1
    mfrP2C(1, 2, IBC_UDIRICHLET) = b1
    mfrP2C(1, 3, IBC_UDIRICHLET) = c1

    mfpP2C(2, 1, IBC_UDIRICHLET) = alpha2
    mfpP2C(2, 2, IBC_UDIRICHLET) = ONE
    mfpP2C(2, 3, IBC_UDIRICHLET) = alpha2
    mfrP2C(2, 1, IBC_UDIRICHLET) = a2 / TWO
    mfrP2C(2, 2, IBC_UDIRICHLET) = ZERO ! not used
    mfrP2C(2, 3, IBC_UDIRICHLET) = ZERO ! not used

    mfpP2C(3, 1, IBC_UDIRICHLET) = alpha
    mfpP2C(3, 2, IBC_UDIRICHLET) = ONE
    mfpP2C(3, 3, IBC_UDIRICHLET) = alpha
    mfrP2C(3, 1, IBC_UDIRICHLET) = a / TWO
    mfrP2C(3, 2, IBC_UDIRICHLET) = b / TWO
    mfrP2C(3, 3, IBC_UDIRICHLET) = ZERO ! not used

    mfpP2C(4, 1, IBC_UDIRICHLET) = alpha2
    mfpP2C(4, 2, IBC_UDIRICHLET) = ONE
    mfpP2C(4, 3, IBC_UDIRICHLET) = alpha2
    mfrP2C(4, 1, IBC_UDIRICHLET) = a2 / TWO
    mfrP2C(4, 2, IBC_UDIRICHLET) = ZERO ! not used
    mfrP2C(4, 3, IBC_UDIRICHLET) = ZERO ! not used

    mfpP2C(5, 1, IBC_UDIRICHLET) = alpha1
    mfpP2C(5, 2, IBC_UDIRICHLET) = ONE
    mfpP2C(5, 3, IBC_UDIRICHLET) = ZERO ! not used
    mfrP2C(5, 1, IBC_UDIRICHLET) = a1
    mfrP2C(5, 2, IBC_UDIRICHLET) = b1
    mfrP2C(5, 3, IBC_UDIRICHLET) = c1

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

    mfpC2P(1, 1, IBC_UDIRICHLET) = ZERO ! not used
    mfpC2P(1, 2, IBC_UDIRICHLET) = ONE
    mfpC2P(1, 3, IBC_UDIRICHLET) = alpha1
    mfrC2P(1, 1, IBC_UDIRICHLET) = a1
    mfrC2P(1, 2, IBC_UDIRICHLET) = b1
    mfrC2P(1, 3, IBC_UDIRICHLET) = c1

    mfpC2P(2, 1, IBC_UDIRICHLET) = alpha2
    mfpC2P(2, 2, IBC_UDIRICHLET) = ONE
    mfpC2P(2, 3, IBC_UDIRICHLET) = alpha2
    mfrC2P(2, 1, IBC_UDIRICHLET) = a2 / TWO
    mfrC2P(2, 2, IBC_UDIRICHLET) = ZERO ! not used
    mfrC2P(2, 3, IBC_UDIRICHLET) = ZERO ! not used

    mfpC2P(3, 1, IBC_UDIRICHLET) = alpha
    mfpC2P(3, 2, IBC_UDIRICHLET) = ONE
    mfpC2P(3, 3, IBC_UDIRICHLET) = alpha
    mfrC2P(3, 1, IBC_UDIRICHLET) = a / TWO
    mfrC2P(3, 2, IBC_UDIRICHLET) = b / TWO
    mfrC2P(3, 3, IBC_UDIRICHLET) = ZERO ! not used

    mfpC2P(4, 1, IBC_UDIRICHLET) = alpha2
    mfpC2P(4, 2, IBC_UDIRICHLET) = ONE
    mfpC2P(4, 3, IBC_UDIRICHLET) = alpha2
    mfrC2P(4, 1, IBC_UDIRICHLET) = a2 / TWO
    mfrC2P(4, 2, IBC_UDIRICHLET) = ZERO ! not used
    mfrC2P(4, 3, IBC_UDIRICHLET) = ZERO ! not used

    mfpC2P(5, 1, IBC_UDIRICHLET) = alpha1
    mfpC2P(5, 2, IBC_UDIRICHLET) = ONE
    mfpC2P(5, 3, IBC_UDIRICHLET) = ZERO ! not used
    mfrC2P(5, 1, IBC_UDIRICHLET) = a1
    mfrC2P(5, 2, IBC_UDIRICHLET) = b1
    mfrC2P(5, 3, IBC_UDIRICHLET) = c1
    
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
    use input_general_mod, only: IBC_PERIODIC
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
    use input_general_mod, only: IBC_PERIODIC
    use udf_type_mod
    use parameters_constant_mod, only: ZERO
    implicit none

    type(t_domain), intent(in) :: d

    integer(4) :: i

!-------------------------------------------------------------------------------
! derivative in x direction with nc unknows
!-------------------------------------------------------------------------------
    i = 1
    ! c2c with nc unknows
    allocate (adx_C2C (d%nc(i) ) ); adx_C2C(:) = ZERO
    allocate (bdx_C2C (d%nc(i) ) ); bdx_C2C(:) = ZERO
    allocate (cdx_C2C (d%nc(i) ) ); cdx_C2C(:) = ZERO
    allocate (ddx_C2C (d%nc(i) ) ); ddx_C2C(:) = ZERO
    call Buildup_TDMA_LHS_array(d%nc(i), d%bc(:, i), cfpC2C, &
        adx_C2C, bdx_C2C, cdx_C2C, ddx_C2C)
    ! P2C with nc unknows
    allocate (adx_P2C (d%nc(i) ) ); adx_P2C(:) = ZERO
    allocate (bdx_P2C (d%nc(i) ) ); bdx_P2C(:) = ZERO
    allocate (cdx_P2C (d%nc(i) ) ); cdx_P2C(:) = ZERO
    allocate (ddx_P2C (d%nc(i) ) ); ddx_P2C(:) = ZERO
    call Buildup_TDMA_LHS_array(d%nc(i), d%bc(:, i), cfpP2C, &
        adx_P2C, bdx_P2C, cdx_P2C, ddx_P2C)
    
    ! derivative in x direction with np unknows
    allocate (adx_P2P (d%np(i) ) ); adx_P2P(:) = ZERO
    allocate (bdx_P2P (d%np(i) ) ); bdx_P2P(:) = ZERO
    allocate (cdx_P2P (d%np(i) ) ); cdx_P2P(:) = ZERO
    allocate (ddx_P2P (d%np(i) ) ); ddx_P2P(:) = ZERO
    call Buildup_TDMA_LHS_array(d%np(i), d%bc(:, i), cfpP2P, &
        adx_P2P, bdx_P2P, cdx_P2P, ddx_P2P)

    allocate (adx_C2P (d%np(i) ) ); adx_C2P(:) = ZERO
    allocate (bdx_C2P (d%np(i) ) ); bdx_C2P(:) = ZERO
    allocate (cdx_C2P (d%np(i) ) ); cdx_C2P(:) = ZERO
    allocate (ddx_C2P (d%np(i) ) ); ddx_C2P(:) = ZERO
    call Buildup_TDMA_LHS_array(d%np(i), d%bc(:, i), cfpC2P, &
        adx_C2P, bdx_C2P, cdx_C2P, ddx_C2P)

!-------------------------------------------------------------------------------
! derivative in y direction with nc unknows
!-------------------------------------------------------------------------------
    i = 2
    allocate (ady_C2C (d%nc(i) ) ); ady_C2C(:) = ZERO
    allocate (bdy_C2C (d%nc(i) ) ); bdy_C2C(:) = ZERO
    allocate (cdy_C2C (d%nc(i) ) ); cdy_C2C(:) = ZERO
    allocate (ddy_C2C (d%nc(i) ) ); ddy_C2C(:) = ZERO
    call Buildup_TDMA_LHS_array(d%nc(i), d%bc(:, i), cfpC2C, &
        ady_C2C, bdy_C2C, cdy_C2C, ddy_C2C)

    allocate (ady_P2C (d%nc(i) ) ); ady_P2C(:) = ZERO
    allocate (bdy_P2C (d%nc(i) ) ); bdy_P2C(:) = ZERO
    allocate (cdy_P2C (d%nc(i) ) ); cdy_P2C(:) = ZERO
    allocate (ddy_P2C (d%nc(i) ) ); ddy_P2C(:) = ZERO
    call Buildup_TDMA_LHS_array(d%nc(i), d%bc(:, i), cfpP2C, &
        ady_P2C, bdy_P2C, cdy_P2C, ddy_P2C)

!-------------------------------------------------------------------------------
! derivative in y direction with np unknows
!-------------------------------------------------------------------------------
    allocate (ady_P2P (d%np(i) ) ); ady_P2P(:) = ZERO
    allocate (bdy_P2P (d%np(i) ) ); bdy_P2P(:) = ZERO
    allocate (cdy_P2P (d%np(i) ) ); cdy_P2P(:) = ZERO
    allocate (ddy_P2P (d%np(i) ) ); ddy_P2P(:) = ZERO
    call Buildup_TDMA_LHS_array(d%np(i), d%bc(:, i), cfpP2P, &
        ady_P2P, bdy_P2P, cdy_P2P, ddy_P2P)

    allocate (ady_C2P (d%np(i) ) ); ady_C2P(:) = ZERO
    allocate (bdy_C2P (d%np(i) ) ); bdy_C2P(:) = ZERO
    allocate (cdy_C2P (d%np(i) ) ); cdy_C2P(:) = ZERO
    allocate (ddy_C2P (d%np(i) ) ); ddy_C2P(:) = ZERO
    call Buildup_TDMA_LHS_array(d%np(i), d%bc(:, i), cfpC2P, &
        ady_C2P, bdy_C2P, cdy_C2P, ddy_C2P)
    
!-------------------------------------------------------------------------------
! derivative in z direction with nc unknows
!-------------------------------------------------------------------------------
    i = 3
    allocate (adz_C2C (d%nc(i) ) ); adz_C2C(:) = ZERO
    allocate (bdz_C2C (d%nc(i) ) ); bdz_C2C(:) = ZERO
    allocate (cdz_C2C (d%nc(i) ) ); cdz_C2C(:) = ZERO
    allocate (ddz_C2C (d%nc(i) ) ); ddz_C2C(:) = ZERO
    call Buildup_TDMA_LHS_array(d%nc(i), d%bc(:, i), cfpC2C, &
        adz_C2C, bdz_C2C, cdz_C2C, ddz_C2C)

    allocate (adz_P2C (d%nc(i) ) ); adz_P2C(:) = ZERO
    allocate (bdz_P2C (d%nc(i) ) ); bdz_P2C(:) = ZERO
    allocate (cdz_P2C (d%nc(i) ) ); cdz_P2C(:) = ZERO
    allocate (ddz_P2C (d%nc(i) ) ); ddz_P2C(:) = ZERO
    call Buildup_TDMA_LHS_array(d%nc(i), d%bc(:, i), cfpP2C, &
        adz_P2C, bdz_P2C, cdz_P2C, ddz_P2C)

!-------------------------------------------------------------------------------
! derivative in z direction with np unknows
!-------------------------------------------------------------------------------
    allocate (adz_P2P (d%np(i) ) ); adz_P2P(:) = ZERO
    allocate (bdz_P2P (d%np(i) ) ); bdz_P2P(:) = ZERO
    allocate (cdz_P2P (d%np(i) ) ); cdz_P2P(:) = ZERO
    allocate (ddz_P2P (d%np(i) ) ); ddz_P2P(:) = ZERO
    call Buildup_TDMA_LHS_array(d%np(i), d%bc(:, i), cfpP2P, &
        adz_P2P, bdz_P2P, cdz_P2P, ddz_P2P)

    allocate (adz_C2P (d%np(i) ) ); adz_C2P(:) = ZERO
    allocate (bdz_C2P (d%np(i) ) ); bdz_C2P(:) = ZERO
    allocate (cdz_C2P (d%np(i) ) ); cdz_C2P(:) = ZERO
    allocate (ddz_C2P (d%np(i) ) ); ddz_C2P(:) = ZERO
    call Buildup_TDMA_LHS_array(d%np(i), d%bc(:, i), cfpC2P, &
        adz_C2P, bdz_C2P, cdz_C2P, ddz_C2P)

!-------------------------------------------------------------------------------
! mid-point interpolation in x direction with nc unknows
!-------------------------------------------------------------------------------
    i = 1
    allocate (amx_P2C (d%nc(i) ) ); amx_P2C(:) = ZERO
    allocate (bmx_P2C (d%nc(i) ) ); bmx_P2C(:) = ZERO
    allocate (cmx_P2C (d%nc(i) ) ); cmx_P2C(:) = ZERO
    allocate (dmx_P2C (d%nc(i) ) ); dmx_P2C(:) = ZERO
    call Buildup_TDMA_LHS_array(d%nc(i), d%bc(:, i), mfpP2C, &
        amx_P2C, bmx_P2C, cmx_P2C, dmx_P2C)

!-------------------------------------------------------------------------------
! mid-point interpolation in x direction with np unknows
!-------------------------------------------------------------------------------
    allocate (amx_C2P (d%np(i) ) ); amx_C2P(:) = ZERO
    allocate (bmx_C2P (d%np(i) ) ); bmx_C2P(:) = ZERO
    allocate (cmx_C2P (d%np(i) ) ); cmx_C2P(:) = ZERO
    allocate (dmx_C2P (d%np(i) ) ); dmx_C2P(:) = ZERO
    call Buildup_TDMA_LHS_array(d%np(i), d%bc(:, i), mfpC2P, &
        amx_C2P, bmx_C2P, cmx_C2P, dmx_C2P)

!-------------------------------------------------------------------------------
! mid-point interpolation in y direction with nc unknows
!-------------------------------------------------------------------------------
    i = 2
    allocate (amy_P2C (d%nc(i) ) ); amy_P2C(:) = ZERO
    allocate (bmy_P2C (d%nc(i) ) ); bmy_P2C(:) = ZERO
    allocate (cmy_P2C (d%nc(i) ) ); cmy_P2C(:) = ZERO
    allocate (dmy_P2C (d%nc(i) ) ); dmy_P2C(:) = ZERO
    call Buildup_TDMA_LHS_array(d%nc(i), d%bc(:, i), mfpP2C, &
        amy_P2C, bmy_P2C, cmy_P2C, dmy_P2C)

!-------------------------------------------------------------------------------
! mid-point interpolation in y direction with np unknows
!-------------------------------------------------------------------------------
    allocate (amy_C2P (d%np(i) ) ); amy_C2P(:) = ZERO
    allocate (bmy_C2P (d%np(i) ) ); bmy_C2P(:) = ZERO
    allocate (cmy_C2P (d%np(i) ) ); cmy_C2P(:) = ZERO
    allocate (dmy_C2P (d%np(i) ) ); dmy_C2P(:) = ZERO
    call Buildup_TDMA_LHS_array(d%np(i), d%bc(:, i), mfpC2P, &
        amy_C2P, bmy_C2P, cmy_C2P, dmy_C2P)

!-------------------------------------------------------------------------------
! mid-point interpolation in z direction with nc unknows
!-------------------------------------------------------------------------------
    i = 3
    allocate (amz_P2C (d%nc(i) ) ); amz_P2C(:) = ZERO
    allocate (bmz_P2C (d%nc(i) ) ); bmz_P2C(:) = ZERO
    allocate (cmz_P2C (d%nc(i) ) ); cmz_P2C(:) = ZERO
    allocate (dmz_P2C (d%nc(i) ) ); dmz_P2C(:) = ZERO
    call Buildup_TDMA_LHS_array(d%nc(i), d%bc(:, i), mfpP2C, &
        amz_P2C, bmz_P2C, cmz_P2C, dmz_P2C)

!-------------------------------------------------------------------------------
! mid-point interpolation in z direction with np unknows
!-------------------------------------------------------------------------------
    allocate (amz_C2P (d%np(i) ) ); amz_C2P(:) = ZERO
    allocate (bmz_C2P (d%np(i) ) ); bmz_C2P(:) = ZERO
    allocate (cmz_C2P (d%np(i) ) ); cmz_C2P(:) = ZERO
    allocate (dmz_C2P (d%np(i) ) ); dmz_C2P(:) = ZERO
    call Buildup_TDMA_LHS_array(d%np(i), d%bc(:, i), mfpC2P, &
        amz_C2P, bmz_C2P, cmz_C2P, dmz_C2P)

    return
  end subroutine Prepare_TDMA_LHS_matrix

  subroutine Prepare_coeffs_for_operations
    use input_general_mod, only: i1deriAccu, i2deriAccu
    use geometry_mod, only: domain

    call Assign_TDMA_coeffs(i1deriAccu)
    call Prepare_TDMA_LHS_matrix(domain)
    return
  end subroutine
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
  subroutine Prepare_TDMA_deri_RHS_array(str, n, bc, inbr, dd, coeff, fi, fo)
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

    fo(:) = fo(:) / dd

    return
  end subroutine Prepare_TDMA_deri_RHS_array
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
  subroutine Get_midp_interpolation(str1, str2, d, fi, fo)
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
            mfrP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), amx_P2C(:), bmx_P2C(:), cmx_P2C(:), dmx_P2C(:), d%nc(i))
        
      else if (str2 == 'C2P') then

        call Prepare_TDMA_interp_RHS_array(str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
            mfrC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), amx_C2P(:), bmx_C2P(:), cmx_C2P(:), dmx_C2P(:), d%np(i))

      else
        call Print_error_msg("108: Error input in prepare_FD_TDMA_RHS.")
      end if

    else if (str1 == 'y') then
      i = 2

      if (str2 == 'P2C') then
        call Prepare_TDMA_interp_RHS_array(str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
            mfrP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), amy_P2C(:), bmy_P2C(:), cmy_P2C(:), dmy_P2C(:), d%nc(i))

      else if (str2 == 'C2P') then

        call Prepare_TDMA_interp_RHS_array(str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
            mfrC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), amy_C2P(:), bmy_C2P(:), cmy_C2P(:), dmy_C2P(:), d%np(i))

      else
        call Print_error_msg("109: Error input in prepare_FD_TDMA_RHS.")
      end if

    else if (str1 == 'z') then

      i = 3

      if (str2 == 'P2C') then
        call Prepare_TDMA_interp_RHS_array(str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
            mfrP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), amz_P2C(:), bmz_P2C(:), cmz_P2C(:), dmz_P2C(:), d%nc(i))

      else if (str2 == 'C2P') then

        call Prepare_TDMA_interp_RHS_array(str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
            mfrC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), amz_C2P(:), bmz_C2P(:), cmz_C2P(:), dmz_C2P(:), d%np(i))

      else
        call Print_error_msg("110: Error input in prepare_FD_TDMA_RHS.")
      end if

    else
      call Print_error_msg("111: No such direction.")
    end if

    return 
  end subroutine Get_midp_interpolation
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
  subroutine Get_1st_derivative(str1, str2, d, fi, fo)
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
!===============================================================================
! Code
!===============================================================================
    nsz = size(fo)

    if(str1=='x') then
      i = 1

      if (str2 == 'C2C') then

        call Prepare_TDMA_deri_RHS_array( str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
              d%h(i), cfrC2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(i), fo(:), adx_C2C(:), bdx_C2C(:), cdx_C2C(:), ddx_C2C(:), d%nc(i) )

      else if (str2 == 'P2C') then
        call Prepare_TDMA_deri_RHS_array(str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
            d%h(i), cfrP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), adx_P2C(:), bdx_P2C(:), cdx_P2C(:), ddx_P2C(:), d%nc(i))

      else if (str2 == 'P2P') then
        
        call Prepare_TDMA_deri_RHS_array(str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
            d%h(i), cfrP2P(:, :, :), fi(:), fo(:) )
            !write(*,'(A,7F8.4)') 'a', adx_P2P(:)
            !write(*,'(A,7F8.4)') 'b', bdx_P2P(:)
            !write(*,'(A,7F8.4)') 'c', cdx_P2P(:)
            !write(*,'(A,7F8.4)') 'd', ddx_P2P(:)
            !write(*,'(A,7F8.4)') 'r', fo(:)
        call Solve_TDMA(d%is_periodic(i), fo(:), adx_P2P(:), bdx_P2P(:), cdx_P2P(:), ddx_P2P(:), d%np(i))
        !write(*,'(A,7F8.4)') 'o', fo(:)
      else if (str2 == 'C2P') then

        call Prepare_TDMA_deri_RHS_array(str2, nsz, d%bc(:, i), d%iNeighb(:, :), &
            d%h(i), cfrC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), adx_C2P(:), bdx_C2P(:), cdx_C2P(:), ddx_C2P(:), d%np(i))

      else
        call Print_error_msg("112: No such staggered scheme defined")
      end if

    else if (str1 == 'y') then
      i = 2

      if (str2 == 'C2C') then

        call Prepare_TDMA_deri_RHS_array( str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
              d%h(i), cfrC2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(i), fo(:), ady_C2C(:), bdy_C2C(:), cdy_C2C(:), ddy_C2C(:), d%nc(i) )
        if(d%is_stretching(2)) fo(:) = fo(:) * d%yMappingcc(:, 1)
      
      else if (str2 == 'P2C') then
        call Prepare_TDMA_deri_RHS_array(str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
            d%h(i), cfrP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), ady_P2C(:), bdy_P2C(:), cdy_P2C(:), ddy_P2C(:), d%nc(i))
        if(d%is_stretching(2)) fo(:) = fo(:) * d%yMappingcc(:, 1)

      else if (str2 == 'P2P') then

        call Prepare_TDMA_deri_RHS_array(str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
            d%h(i), cfrP2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), ady_P2P(:), bdy_P2P(:), cdy_P2P(:), ddy_P2P(:), d%np(i))
        if(d%is_stretching(2)) fo(:) = fo(:) * d%yMappingpt(:, 1)

      else if (str2 == 'C2P') then

        call Prepare_TDMA_deri_RHS_array(str2, nsz, d%bc(:, i), d%jNeighb(:, :), &
            d%h(i), cfrC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), ady_C2P(:), bdy_C2P(:), cdy_C2P(:), ddy_C2P(:), d%np(i))
        if(d%is_stretching(2)) fo(:) = fo(:) * d%yMappingpt(:, 1)

      else
        call Print_error_msg("113: No such staggered scheme defined")
      end if

    else if (str1 == 'z') then

      i = 3

      if (str2 == 'C2C') then

        call Prepare_TDMA_deri_RHS_array( str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
              d%h(i), cfrC2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA( d%is_periodic(i), fo(:), adz_C2C(:), bdz_C2C(:), cdz_C2C(:), ddz_C2C(:), d%nc(i) )

      else if (str2 == 'P2C') then
        call Prepare_TDMA_deri_RHS_array(str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
            d%h(i), cfrP2C(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), adz_P2C(:), bdz_P2C(:), cdz_P2C(:), ddz_P2C(:), d%nc(i))

      else if (str2 == 'P2P') then

        call Prepare_TDMA_deri_RHS_array(str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
            d%h(i), cfrP2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), adz_P2P(:), bdz_P2P(:), cdz_P2P(:), ddz_P2P(:), d%np(i))

      else if (str2 == 'C2P') then

        call Prepare_TDMA_deri_RHS_array(str2, nsz, d%bc(:, i), d%kNeighb(:, :), &
            d%h(i), cfrC2P(:, :, :), fi(:), fo(:) )
        call Solve_TDMA(d%is_periodic(i), fo(:), adz_C2P(:), bdz_C2P(:), cdz_C2P(:), ddz_C2P(:), d%np(i))

      else
        call Print_error_msg("114: No such staggered scheme defined")
      end if

    else
      call Print_error_msg("115: No such direction defined.")
    end if

    return 
  end subroutine Get_1st_derivative

!===============================================================================
!===============================================================================
!> \brief To test this subroutine for mid-point interpolation.
!>
!> This subroutine is called in \ref Test_schemes. Define the logicals to choose
!> which test section is required. 
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!_______________________________________________________________________________
  subroutine Test_interpolation(d)
    use flow_variables_mod
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in) :: d
    real(WP), allocatable :: fi(:), fo(:)
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    real(WP) :: ref
    integer(4) :: i, j, k
    real(WP) :: err(3), errmax
    logical :: dbg = .false.
    logical :: uix_p2c = .false.
    logical :: vix_c2p = .false.
    logical :: uiy_c2p = .true.
    logical :: viy_p2c = .true.

    if(uix_p2c) then
      !test interpolation. u in x, P2C
      !(i', j, k) --> (i, j, k)
      allocate ( fi( d%np(1) ) ); fi = ZERO
      allocate ( fo( d%nc(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# Test interp u in x P2C: kji, err1, err2, err3'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = qx(:, j, k)
          call Get_midp_interpolation('x', 'P2C', d, fi(:), fo(:))

          do i = 4, d%nc(1)-3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do

          do i = 1, 3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do

          do i = d%nc(1)-2, d%nc(1)
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do

        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3(","ES15.7))') err(1:3)
    end if

    if(uiy_c2p) then
      ! test interpolation. u in y, C2P 
      ! (i', j, k) --> (i', j', k)
      allocate ( fi( d%nc(2) ) ); fi = ZERO
      allocate ( fo( d%np(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# Test interp u in y C2P: kij, err1, err2, err3'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          fi(:) = qx(i, :, k)
          call Get_midp_interpolation('y', 'C2P', d, fi(:), fo(:))
          do j = 4, d%np(2)-3
            yp = d%yp(j)
            ref = sin_wp ( xp ) + sin_wp(yp) + sin_wp(zc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yp = d%yp(j)
            ref = sin_wp ( xp ) + sin_wp(yp) + sin_wp(zc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%np(2)-2, d%np(2)
            yp = d%yp(j)
            ref = sin_wp ( xp ) + sin_wp(yp) + sin_wp(zc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3(","ES15.7))') err(1:3)
    end if

    if(vix_c2p) then
      !test interpolation. v in x, C2P
      !(i', j, k) --> (i, j, k)
      allocate ( fi( d%nc(1) ) ); fi = ZERO
      allocate ( fo( d%np(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# Test interp v in x C2P: err1, err2, err3'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = qy(:, j, k)
          call Get_midp_interpolation('x', 'C2P', d, fi(:), fo(:))
          do i = 4, d%np(1)-3
            xp = d%h(1) * (real(i - 1, WP))
            ref = sin_wp ( xp ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do

          do i = 1, 3
            xp = d%h(1) * (real(i - 1, WP))
            ref = sin_wp ( xp ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do

          do i = d%np(1)-2, d%np(1)
            xp = d%h(1) * (real(i - 1, WP))
            ref = sin_wp ( xp ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do

        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3(","ES15.7))') err(1:3)
    end if

    if(viy_p2c) then
      ! test interpolation. v in y, P2C 
      ! (i', j, k) --> (i', j', k)
      allocate ( fi( d%np(2) ) ); fi = ZERO
      allocate ( fo( d%nc(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# Test interp v in y P2C: kij, err1, err2, err3'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%nc(1)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          fi(:) = qy(i, :, k)
          call Get_midp_interpolation('y', 'P2C', d, fi(:), fo(:))
          do j = 4, d%nc(2)-3
            yc = d%yc(j)
            ref = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yc = d%yc(j)
            ref = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%nc(2)-2, d%nc(2)
            yc = d%yc(j)
            ref = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3(","ES15.7))') err(1:3)
    end if

    return 
  end subroutine
!===============================================================================
!===============================================================================
!> \brief To test this subroutine for 1st derivative.
!>
!> This subroutine is called in \ref Test_schemes. Define the logicals to choose
!> which test section is required. 
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!_______________________________________________________________________________
  subroutine Test_1st_derivative(d)
    use flow_variables_mod
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in) :: d
    real(WP), allocatable :: fi(:), fo(:)
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    real(WP) :: ref
    integer(4) :: i, j, k
    real(WP) :: err(3), errmax
    logical :: dbg = .false.

    logical :: dudx_P2C = .false.
    logical :: dudx_P2P = .false.
    logical :: dvdx_C2P = .false.
    logical :: dvdx_C2C = .false.

    logical :: dudy_C2P = .true.
    logical :: dudy_C2C = .true.
    logical :: dvdy_P2P = .true.
    logical :: dvdy_P2C = .true.

    if(dudx_P2C) then
      ! du / dx, P2C
      ! (i', j, k) --> (i, j, k)
      allocate ( fi( d%np(1) ) ); fi = ZERO
      allocate ( fo( d%nc(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# Test 1stDeri u in x P2C: I, err1, err2, err3'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = qx(:, j, k)
          call Get_1st_derivative('x', 'P2C', d, fi(:), fo(:))
          do i = 4, d%nc(1)-3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = cos_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = 1, 3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = cos_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = d%nc(1)-2, d%nc(1)
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = cos_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3(","ES15.7))') err(1:3)
    end if

    if(dudx_P2P) then
    ! du / dx, P2P
    ! (i', j, k) --> (i', j, k)
      allocate ( fi( d%np(1) ) ); fi = ZERO
      allocate ( fo( d%np(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# Test 1stDeri u in x P2P: I, err1, err2, err3'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = qx(:, j, k)
          call Get_1st_derivative('x', 'P2P', d, fi(:), fo(:))
          do i = 4, d%np(1)-3
            xp = d%h(1) * real(i - 1, WP)
            ref = cos_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = 1, 3
            xp = d%h(1) * real(i - 1, WP)
            ref = cos_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = d%np(1)-2, d%np(1)
            xp = d%h(1) * real(i - 1, WP)
            ref = cos_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3(","ES15.7))') err(1:3)
    end if

    if(dudy_C2P) then
      ! du / dy, C2P
      ! (i', j, k) --> (i', j', k)
      allocate ( fi( d%nc(2) ) ); fi = ZERO
      allocate ( fo( d%np(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# Test 1stDeri u in y C2P: I, err1, err2, err3'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          fi(:) = qx(i, :, k)
          call Get_1st_derivative('y', 'C2P', d, fi(:), fo(:))
          do j = 4, d%np(2)-3
            yp = d%yp(j)
            ref = cos_wp(yp)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yp = d%yp(j)
            ref = cos_wp(yp)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%np(2)-2, d%np(2)
            yp = d%yp(j)
            ref = cos_wp(yp)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3(","ES15.7))') err(1:3)   
    end if

    if(dudy_C2C) then
      ! du / dy, C2C
      ! (i', j, k) --> (i, j, k)
      allocate ( fi( d%nc(2) ) ); fi = ZERO
      allocate ( fo( d%nc(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# Test 1stDeri u in y C2C: I, err1, err2, err3'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          fi(:) = qx(i, :, k)
          call Get_1st_derivative('y', 'C2C', d, fi(:), fo(:))
          do j = 4, d%nc(2)-3
            yc = d%yc(j)
            ref = cos_wp(yc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yc = d%yc(j)
            ref = cos_wp(yc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%nc(2)-2, d%nc(2)
            yc = d%yc(j)
            ref = cos_wp(yc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3(","ES15.7))') err(1:3)
    end if

    if(dvdy_P2C) then
      ! dv / dy, P2C
      ! (i, j', k) --> (i, j, k)
      allocate ( fi( d%np(2) ) ); fi = ZERO
      allocate ( fo( d%nc(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# Test 1stDeri v in y P2C: I, err1, err2, err3'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%nc(1)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          fi(:) = qy(i, :, k)
          call Get_1st_derivative('y', 'P2C', d, fi(:), fo(:))
          do j = 4, d%nc(2)-3
            yc = d%yc(j)
            ref = cos_wp ( yc )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yc = d%yc(j)
            ref = cos_wp ( yc )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%nc(2)-2, d%nc(2)
            yc = d%yc(j)
            ref = cos_wp ( yc )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3(","ES15.7))') err(1:3) 
    end if 
    
    if(dvdy_P2P) then
      ! dv / dy, P2P
      ! (i, j', k) --> (i, j', k)
      allocate ( fi( d%np(2) ) ); fi = ZERO
      allocate ( fo( d%np(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# Test 1stDeri v in y P2P: I, err1, err2, err3'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%nc(1)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          fi(:) = qy(i, :, k)
          call Get_1st_derivative('y', 'P2P', d, fi(:), fo(:))
          do j = 4, d%np(2)-3
            yp = d%yp(j)
            ref = cos_wp ( yp )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yp = d%yp(j)
            ref = cos_wp ( yp )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%np(2)-2, d%np(2)
            yp = d%yp(j)
            ref = cos_wp ( yp )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3(","ES15.7))') err(1:3) 
    end if

    if(dvdx_C2C) then
      ! du / dx, P2C
      ! (i', j, k) --> (i, j, k)
      allocate ( fi( d%nc(1) ) ); fi = ZERO
      allocate ( fo( d%nc(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# Test 1stDeri u in x P2C: I, err1, err2, err3'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = qy(:, j, k)
          call Get_1st_derivative('x', 'C2C', d, fi(:), fo(:))
          do i = 4, d%nc(1)-3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = cos_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = 1, 3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = cos_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = d%nc(1)-2, d%nc(1)
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = cos_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3(","ES15.7))') err(1:3)
    end if

    if(dvdx_C2P) then
    ! du / dx, P2P
    ! (i', j, k) --> (i', j, k)
      allocate ( fi( d%nc(1) ) ); fi = ZERO
      allocate ( fo( d%np(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# Test 1stDeri u in x P2P: I, err1, err2, err3'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = qx(:, j, k)
          call Get_1st_derivative('x', 'P2P', d, fi(:), fo(:))
          do i = 4, d%np(1)-3
            xp = d%h(1) * real(i - 1, WP)
            ref = cos_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = 1, 3
            xp = d%h(1) * real(i - 1, WP)
            ref = cos_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = d%np(1)-2, d%np(1)
            xp = d%h(1) * real(i - 1, WP)
            ref = cos_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3(","ES15.7))') err(1:3)
    end if
    
    
    return 
  end subroutine

end module