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
! This program is disatributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!----------------------------------------------------------------------------------------------------------
!==========================================================================================================
!> \file operations.f90
!>
!> \brief A general operation of derivative and interpolation in 1D.
!>
!==========================================================================================================
module operations
  use precision_mod, only : WP
  implicit none

  private
!----------------------------------------------------------------------------------------------------------
! basic coefficients for TDMA of 1st deriviative  
! to store coefficients for TDMA
!     d1fC2C vs d1rC2C :
!       f : coefficients in the LHS, unknown side.
!       r : coefficients in the RHS, known side. 
! eg, d1fC2C(5, 3, 5)
!     First column: 1:2 for L1,2, one side b.c.
!                   4:5 for L4,5, the other side b.c.
!                   3   for L-bulk, all other cells
!     Second column: 1 for coefficients of LHS f^(1)_{i-1}
!                    2 for coefficients of LHS f^(1)_{i}
!                    3 for coefficients of LHS f^(1)_{i+1}
!     Third column:  for b.c. types
!                    0 = IBC_INTERIOR  
!                    1 = IBC_PERIODIC  
!                    2 = IBC_SYMMETRIC 
!                    3 = IBC_ASYMMETRIC
!                    4 = IBC_DIRICHLET 
!                    5 = IBC_NEUMANN   
!                    6 = IBC_INTRPL    
!     Fourth column:  Accuracy 
!                    1 = IACCU_CD2
!                    2 = IACCU_CD4
!                    3 = IACCU_CP4
!                    4 = IACCU_CP6
!----------------------------------------------------------------------------------------------------------
  integer, parameter :: NL = 5 ! rows/line types
  integer, parameter :: NS = 3 ! how many coefficients
  integer, parameter :: NBCS = 0 ! bc index, start
  integer, parameter :: NBCE = 6 ! bc index, end
  integer, parameter :: NACC = 4 ! accuracy types
!----------------------------------------------------------------------------------------------------------
! for 1st derivative
!----------------------------------------------------------------------------------------------------------
  ! collocated C2C
  real(WP), save, public :: d1fC2C(NL, NS, NBCS:NBCE, NACC)
  real(WP), save, public :: d1rC2C(NL, NS, NBCS:NBCE, NACC)
  ! collocated P2P
  real(WP), save, public :: d1fP2P(NL, NS, NBCS:NBCE, NACC)
  real(WP), save, public :: d1rP2P(NL, NS, NBCS:NBCE, NACC)
  ! staggered C2P
  real(WP), save, public :: d1fC2P(NL, NS, NBCS:NBCE, NACC)
  real(WP), save, public :: d1rC2P(NL, NS, NBCS:NBCE, NACC)
  ! staggered P2C
  real(WP), save, public :: d1fP2C(NL, NS, NBCS:NBCE, NACC)
  real(WP), save, public :: d1rP2C(NL, NS, NBCS:NBCE, NACC)
!----------------------------------------------------------------------------------------------------------
! for iterpolation
!----------------------------------------------------------------------------------------------------------
  ! interpolation P2C
  real(WP), save, public :: m1fP2C(NL, NS, NBCS:NBCE, NACC)
  real(WP), save, public :: m1rP2C(NL, NS, NBCS:NBCE, NACC)
  ! interpolation C2P
  real(WP), save, public :: m1fC2P(NL, NS, NBCS:NBCE, NACC)
  real(WP), save, public :: m1rC2P(NL, NS, NBCS:NBCE, NACC)

!----------------------------------------------------------------------------------------------------------
! coefficients array for TDMA of 1st deriviative  
! to store coefficients array for TDMA
!----------------------------------------------------------------------------------------------------------
  type t_xtdma_lhs
!----------------------------------------------------------------------------------------------------------
!   x : pre-processed TDMA LHS Matrix for 1st deriviative
!----------------------------------------------------------------------------------------------------------
    real(WP), allocatable :: ad1x_P2P(:, :, :, :)
    real(WP), allocatable :: bd1x_P2P(:, :, :, :)
    real(WP), allocatable :: cd1x_P2P(:, :, :, :)
    real(WP), allocatable :: dd1x_P2P(:, :, :, :)

    real(WP), allocatable :: ad1x_C2C(:, :, :, :)
    real(WP), allocatable :: bd1x_C2C(:, :, :, :)
    real(WP), allocatable :: cd1x_C2C(:, :, :, :)
    real(WP), allocatable :: dd1x_C2C(:, :, :, :)

    real(WP), allocatable :: ad1x_P2C(:, :, :, :)
    real(WP), allocatable :: bd1x_P2C(:, :, :, :)
    real(WP), allocatable :: cd1x_P2C(:, :, :, :)
    real(WP), allocatable :: dd1x_P2C(:, :, :, :)

    real(WP), allocatable :: ad1x_C2P(:, :, :, :)
    real(WP), allocatable :: bd1x_C2P(:, :, :, :)
    real(WP), allocatable :: cd1x_C2P(:, :, :, :)
    real(WP), allocatable :: dd1x_C2P(:, :, :, :)
!----------------------------------------------------------------------------------------------------------
!   x : pre-processed TDMA LHS Matrix for mid-point interpolation
!----------------------------------------------------------------------------------------------------------
    real(WP), allocatable :: am1x_P2C(:, :, :, :)
    real(WP), allocatable :: bm1x_P2C(:, :, :, :)
    real(WP), allocatable :: cm1x_P2C(:, :, :, :)
    real(WP), allocatable :: dm1x_P2C(:, :, :, :)

    real(WP), allocatable :: am1x_C2P(:, :, :, :)
    real(WP), allocatable :: bm1x_C2P(:, :, :, :)
    real(WP), allocatable :: cm1x_C2P(:, :, :, :)
    real(WP), allocatable :: dm1x_C2P(:, :, :, :)
  end type t_xtdma_lhs

  type(t_xtdma_lhs), allocatable :: xtdma_lhs(:) 

!----------------------------------------------------------------------------------------------------------
! y : pre-processed TDMA LHS Matrix for 1st deriviative
!----------------------------------------------------------------------------------------------------------
  real(WP), allocatable :: ad1y_P2P(:, :, :, :)
  real(WP), allocatable :: bd1y_P2P(:, :, :, :)
  real(WP), allocatable :: cd1y_P2P(:, :, :, :)
  real(WP), allocatable :: dd1y_P2P(:, :, :, :)

  real(WP), allocatable :: ad1y_C2C(:, :, :, :)
  real(WP), allocatable :: bd1y_C2C(:, :, :, :)
  real(WP), allocatable :: cd1y_C2C(:, :, :, :)
  real(WP), allocatable :: dd1y_C2C(:, :, :, :)

  real(WP), allocatable :: ad1y_P2C(:, :, :, :)
  real(WP), allocatable :: bd1y_P2C(:, :, :, :)
  real(WP), allocatable :: cd1y_P2C(:, :, :, :)
  real(WP), allocatable :: dd1y_P2C(:, :, :, :)

  real(WP), allocatable :: ad1y_C2P(:, :, :, :)
  real(WP), allocatable :: bd1y_C2P(:, :, :, :)
  real(WP), allocatable :: cd1y_C2P(:, :, :, :)
  real(WP), allocatable :: dd1y_C2P(:, :, :, :)
!----------------------------------------------------------------------------------------------------------
! y : pre-processed TDMA LHS Matrix for mid-point interpolation
!----------------------------------------------------------------------------------------------------------
  real(WP), allocatable :: am1y_P2C(:, :, :, :)
  real(WP), allocatable :: bm1y_P2C(:, :, :, :)
  real(WP), allocatable :: cm1y_P2C(:, :, :, :)
  real(WP), allocatable :: dm1y_P2C(:, :, :, :)

  real(WP), allocatable :: am1y_C2P(:, :, :, :)
  real(WP), allocatable :: bm1y_C2P(:, :, :, :)
  real(WP), allocatable :: cm1y_C2P(:, :, :, :)
  real(WP), allocatable :: dm1y_C2P(:, :, :, :)
!----------------------------------------------------------------------------------------------------------
! z : pre-processed TDMA LHS Matrix for 1st deriviative
!----------------------------------------------------------------------------------------------------------
  real(WP), allocatable :: ad1z_P2P(:, :, :, :)
  real(WP), allocatable :: bd1z_P2P(:, :, :, :)
  real(WP), allocatable :: cd1z_P2P(:, :, :, :)
  real(WP), allocatable :: dd1z_P2P(:, :, :, :)

  real(WP), allocatable :: ad1z_C2C(:, :, :, :)
  real(WP), allocatable :: bd1z_C2C(:, :, :, :)
  real(WP), allocatable :: cd1z_C2C(:, :, :, :)
  real(WP), allocatable :: dd1z_C2C(:, :, :, :)

  real(WP), allocatable :: ad1z_P2C(:, :, :, :)
  real(WP), allocatable :: bd1z_P2C(:, :, :, :)
  real(WP), allocatable :: cd1z_P2C(:, :, :, :)
  real(WP), allocatable :: dd1z_P2C(:, :, :, :)

  real(WP), allocatable :: ad1z_C2P(:, :, :, :)
  real(WP), allocatable :: bd1z_C2P(:, :, :, :)
  real(WP), allocatable :: cd1z_C2P(:, :, :, :)
  real(WP), allocatable :: dd1z_C2P(:, :, :, :)
!----------------------------------------------------------------------------------------------------------
! z : pre-processed TDMA LHS Matrix for mid-point interpolation
!----------------------------------------------------------------------------------------------------------
  real(WP), allocatable :: am1z_P2C(:, :, :, :)
  real(WP), allocatable :: bm1z_P2C(:, :, :, :)
  real(WP), allocatable :: cm1z_P2C(:, :, :, :)
  real(WP), allocatable :: dm1z_P2C(:, :, :, :)

  real(WP), allocatable :: am1z_C2P(:, :, :, :)
  real(WP), allocatable :: bm1z_C2P(:, :, :, :)
  real(WP), allocatable :: cm1z_C2P(:, :, :, :)
  real(WP), allocatable :: dm1z_C2P(:, :, :, :)
!----------------------------------------------------------------------------------------------------------
! processures
!----------------------------------------------------------------------------------------------------------
  private :: Prepare_compact_coefficients
  private :: Buildup_TDMA_LHS_array
  public  :: Prepare_LHS_coeffs_for_operations

  private :: buildup_ghost_cells_P
  private :: buildup_ghost_cells_C

  private :: Prepare_TDMA_interp_P2C_RHS_array
  private :: Get_x_midp_P2C_1D
  private :: Get_y_midp_P2C_1D
  private :: Get_z_midp_P2C_1D
  public  :: Get_x_midp_P2C_3D
  public  :: Get_y_midp_P2C_3D 
  public  :: Get_z_midp_P2C_3D 

  private :: Prepare_TDMA_interp_C2P_RHS_array
  private :: Get_x_midp_C2P_1D
  private :: Get_y_midp_C2P_1D
  private :: Get_z_midp_C2P_1D
  public  :: Get_x_midp_C2P_3D
  public  :: Get_y_midp_C2P_3D
  public  :: Get_z_midp_C2P_3D

  private :: Prepare_TDMA_1deri_C2C_RHS_array
  private :: Get_x_1st_derivative_C2C_1D
  private :: Get_y_1st_derivative_C2C_1D
  private :: Get_z_1st_derivative_C2C_1D
  public  :: Get_x_1st_derivative_C2C_3D
  public  :: Get_y_1st_derivative_C2C_3D
  public  :: Get_z_1st_derivative_C2C_3D

  private :: Prepare_TDMA_1deri_P2P_RHS_array
  private :: Get_x_1st_derivative_P2P_1D
  private :: Get_y_1st_derivative_P2P_1D
  private :: Get_z_1st_derivative_P2P_1D
  public  :: Get_x_1st_derivative_P2P_3D
  public  :: Get_y_1st_derivative_P2P_3D
  public  :: Get_z_1st_derivative_P2P_3D

  private :: Prepare_TDMA_1deri_C2P_RHS_array
  private :: Get_x_1st_derivative_C2P_1D
  private :: Get_y_1st_derivative_C2P_1D
  private :: Get_z_1st_derivative_C2P_1D
  public  :: Get_x_1st_derivative_C2P_3D
  public  :: Get_y_1st_derivative_C2P_3D
  public  :: Get_z_1st_derivative_C2P_3D

  private :: Prepare_TDMA_1deri_P2C_RHS_array
  private :: Get_x_1st_derivative_P2C_1D
  private :: Get_y_1st_derivative_P2C_1D
  private :: Get_z_1st_derivative_P2C_1D
  public  :: Get_x_1st_derivative_P2C_3D
  public  :: Get_y_1st_derivative_P2C_3D
  public  :: Get_z_1st_derivative_P2C_3D

  public  :: Test_interpolation
  public  :: Test_1st_derivative

contains
!==========================================================================================================
!> \brief Assigned the cooefficients for the compact schemes     
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   private
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                           !
!----------------------------------------------------------------------------------------------------------
!> \param[in]     iaccu         the accuracy given by user
!==========================================================================================================
  subroutine Prepare_compact_coefficients
    use parameters_constant_mod
    use input_general_mod
    use mpi_mod
    implicit none

    real(WP) :: alpha (4),  a(4),  b(4),  c(4)
    real(WP) :: alpha1(4), a1(4), b1(4), c1(4)

    integer :: n

    if(nrank == 0) then
       call Print_debug_start_msg &
         ("Assigning coefficient matrix for the compact schemes ...")
       !write(*, *) "The given numerical accuracy =", iaccu
    end if

!----------------------------------------------------------------------------------------------------------
!   initialisation
!----------------------------------------------------------------------------------------------------------
    d1fC2C(:, :, :, :) = MAXP
    d1rC2C(:, :, :, :) = MAXP
    d1fP2P(:, :, :, :) = MAXP
    d1rP2P(:, :, :, :) = MAXP

    d1fC2P(:, :, :, :) = MAXP
    d1rC2P(:, :, :, :) = MAXP
    d1fP2C(:, :, :, :) = MAXP
    d1rP2C(:, :, :, :) = MAXP

    m1fC2P(:, :, :, :) = MAXP
    m1rC2P(:, :, :, :) = MAXP
    m1fP2C(:, :, :, :) = MAXP
    m1rP2C(:, :, :, :) = MAXP
!==========================================================================================================
! Set 1 : C2C, periodic & symmetric & asymmetric
!         1st derivative on collocated grids, C2C/P2P bulk coefficients
! d1fC2C : "d1"=first deriviative, "f"=f'  side, "C2C"= center 2 centre 
! d1rC2C : "d1"=first deriviative, "r"=rhs side, "C2C"= center 2 centre 
! alpha * f'_{i-1} + f'_i + alpha * f'_{i+1} = a/(2h) * ( f_{i+1} - f_{i-1} ) + &
!                                              b/(4h) * ( f_{i+2} - f_{i-2} )
! C2C unknows: 
! when i=1,    need: LHS: f'_0;      RHS: f_0, f_{-1}
! when i=2,    need:                 RHS: f_0
! when i=nc-1, need:                 RHS: f_{nc+1}
! when i=nc,   need: LHS: f'_{nc+1}; RHS: f_{nc+1}, f_{nc+2}
!==========================================================================================================
    alpha(IACCU_CD2) = ZERO
        a(IACCU_CD2) = ONE
        b(IACCU_CD2) = ZERO

    alpha(IACCU_CD4) = ZERO
        a(IACCU_CD4) = FOUR * ONE_THIRD
        b(IACCU_CD4) = - ONE_THIRD

    alpha(IACCU_CP4) = QUARTER
        a(IACCU_CP4) = ONEPFIVE
        b(IACCU_CP4) = ZERO

    alpha(IACCU_CP6) = ONE_THIRD
        a(IACCU_CP6) = FOURTEEN / NINE
        b(IACCU_CP6) = ONE / NINE
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2C, IBC_PERIODIC, unknowns from both rhs and lhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      d1fC2C(1:5, 1, IBC_PERIODIC, n) = alpha(n)
      d1fC2C(1:5, 2, IBC_PERIODIC, n) = ONE
      d1fC2C(1:5, 3, IBC_PERIODIC, n) = alpha(n)
      d1rC2C(1:5, 1, IBC_PERIODIC, n) = a(n) * HALF    ! a/2
      d1rC2C(1:5, 2, IBC_PERIODIC, n) = b(n) * QUARTER ! b/4
    end do 
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2C, IBC_SYMMETRIC, unknowns from both rhs and lhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      d1fC2C(1,   1, IBC_SYMMETRIC, n) = ZERO        ! not used
      d1fC2C(1,   2, IBC_SYMMETRIC, n) = ONE - alpha(n)
      d1fC2C(1,   3, IBC_SYMMETRIC, n) = alpha(n)
      d1fC2C(5,   1, IBC_SYMMETRIC, n) = d1fC2C(1,   3, IBC_SYMMETRIC, n)
      d1fC2C(5,   2, IBC_SYMMETRIC, n) = d1fC2C(1,   2, IBC_SYMMETRIC, n)
      d1fC2C(5,   3, IBC_SYMMETRIC, n) = d1fC2C(1,   1, IBC_SYMMETRIC, n)
      d1fC2C(2:4, :, IBC_SYMMETRIC, n) = d1fC2C(2:4, :, IBC_PERIODIC, n)
      d1rC2C(:,   :, IBC_SYMMETRIC, n) = d1rC2C(:,   :, IBC_PERIODIC, n)
    end do 
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2C, IBC_ASYMMETRIC, unknowns from both rhs and lhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      d1fC2C(1,   1, IBC_ASYMMETRIC, n) = ZERO        ! not used
      d1fC2C(1,   2, IBC_ASYMMETRIC, n) = ONE + alpha(n)
      d1fC2C(1,   3, IBC_ASYMMETRIC, n) = alpha(n)
      d1fC2C(5,   1, IBC_ASYMMETRIC, n) = d1fC2C(1,   3, IBC_ASYMMETRIC, n)
      d1fC2C(5,   2, IBC_ASYMMETRIC, n) = d1fC2C(1,   2, IBC_ASYMMETRIC, n)
      d1fC2C(5,   3, IBC_ASYMMETRIC, n) = d1fC2C(1,   1, IBC_ASYMMETRIC, n)
      d1fC2C(2:4, :, IBC_ASYMMETRIC, n) = d1fC2C(2:4, :, IBC_PERIODIC  , n)
      d1rC2C(:,   :, IBC_ASYMMETRIC, n) = d1rC2C(:,   :, IBC_PERIODIC  , n)
    end do 
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2C, IBC_INTERIOR, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    d1fC2C(:, :, IBC_INTERIOR, :) = d1fC2C(:, :, IBC_PERIODIC, :)
    d1rC2C(:, :, IBC_INTERIOR, :) = d1rC2C(:, :, IBC_PERIODIC, :)

    d1fC2C(1:2, :, IBC_INTERIOR, IACCU_CP4) = d1fC2C(1:2, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    d1rC2C(1:2, :, IBC_INTERIOR, IACCU_CP4) = d1rC2C(1:2, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    d1fC2C(4:5, :, IBC_INTERIOR, IACCU_CP4) = d1fC2C(4:5, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    d1rC2C(4:5, :, IBC_INTERIOR, IACCU_CP4) = d1rC2C(4:5, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD

    d1fC2C(1:2, :, IBC_INTERIOR, IACCU_CP6) = d1fC2C(1:2, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    d1rC2C(1:2, :, IBC_INTERIOR, IACCU_CP6) = d1rC2C(1:2, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    d1fC2C(4:5, :, IBC_INTERIOR, IACCU_CP6) = d1fC2C(4:5, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    d1rC2C(4:5, :, IBC_INTERIOR, IACCU_CP6) = d1rC2C(4:5, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2C IBC_DIRICHLET, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    d1fC2C(:, :, IBC_DIRICHLET, :) = d1fC2C(:, :, IBC_INTERIOR, :)
    d1rC2C(:, :, IBC_DIRICHLET, :) = d1rC2C(:, :, IBC_INTERIOR, :)
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2C IBC_NEUMANN, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    d1fC2C(:, :, IBC_NEUMANN,   :) = d1fC2C(:, :, IBC_INTERIOR, :)
    d1rC2C(:, :, IBC_NEUMANN,   :) = d1rC2C(:, :, IBC_INTERIOR, :)
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2C, IBC_INTRPL, no bc, no reconstuction. exterpolation only. 
!                    f'_1 + alpha * f'_{2}   = a1 * f_1 + b1 * f_2 + c1 * f-3
! alpha * f'_{1}   + f'_2 + alpha * f'_{3}   = a/(2h) * ( f_{i+1} - f_{i-1} ) 
! alpha * f'_{i-1} + f'_i + alpha * f'_{i+1} = a/(2h) * ( f_{i+1} - f_{i-1} ) + &
!                                              b/(4h) * ( f_{i+2} - f_{i-2} )
!----------------------------------------------------------------------------------------------------------
    alpha1(IACCU_CD2) = ZERO
        a1(IACCU_CD2) = - ONE
        b1(IACCU_CD2) = ONE 
        c1(IACCU_CD2) = ZERO

    alpha1(IACCU_CD4) = ZERO
        a1(IACCU_CD4) = -ONEPFIVE
        b1(IACCU_CD4) = TWO
        c1(IACCU_CD4) = -HALF

    alpha1(IACCU_CP4) = ONE
        a1(IACCU_CP4) = - TWO
        b1(IACCU_CP4) = TWO
        c1(IACCU_CP4) = ZERO

    alpha1(IACCU_CP6) = TWO
        a1(IACCU_CP6) = -TWOPFIVE
        b1(IACCU_CP6) = TWO
        c1(IACCU_CP6) = HALF

    do n = 1, NACC
      d1fC2C(1, 1,   IBC_INTRPL, n) = ZERO ! not used
      d1fC2C(1, 2,   IBC_INTRPL, n) = ONE
      d1fC2C(1, 3,   IBC_INTRPL, n) = alpha1(n)
      d1rC2C(1, 1,   IBC_INTRPL, n) = a1(n)
      d1rC2C(1, 2,   IBC_INTRPL, n) = b1(n)
      d1rC2C(1, 3,   IBC_INTRPL, n) = c1(n)

      d1fC2C(5, 1,   IBC_INTRPL, n) =   d1fC2C(1, 3, IBC_INTRPL, n)
      d1fC2C(5, 2,   IBC_INTRPL, n) =   d1fC2C(1, 2, IBC_INTRPL, n)
      d1fC2C(5, 3,   IBC_INTRPL, n) =   d1fC2C(1, 1, IBC_INTRPL, n)
      d1rC2C(5, 1,   IBC_INTRPL, n) = - d1rC2C(1, 1, IBC_INTRPL, n)
      d1rC2C(5, 2,   IBC_INTRPL, n) = - d1rC2C(1, 2, IBC_INTRPL, n)
      d1rC2C(5, 3,   IBC_INTRPL, n) = - d1rC2C(1, 3, IBC_INTRPL, n)

      d1fC2C(3, 1:3, IBC_INTRPL, n) = d1fC2C(3, 1:3, IBC_PERIODIC, n)
      d1rC2C(3, 1:3, IBC_INTRPL, n) = d1rC2C(3, 1:3, IBC_PERIODIC, n)
    end do
    ! check which below method works good! 
    do n = 1, 2
     !d1fC2C(2, 1:3, IBC_INTRPL, n) = d1fC2C(1, 1:3, IBC_INTRPL, n) ! exterpolation, following Line 1
     !d1rC2C(2, 1:3, IBC_INTRPL, n) = d1rC2C(1, 1:3, IBC_INTRPL, n) ! exterpolation, following Line 1
     !d1fC2C(4, 1:3, IBC_INTRPL, n) = d1fC2C(5, 1:3, IBC_INTRPL, n) ! exterpolation, following Line 1
     !d1rC2C(4, 1:3, IBC_INTRPL, n) = d1rC2C(5, 1:3, IBC_INTRPL, n) ! exterpolation, following Line 1

      d1fC2C(2, 1:3, IBC_INTRPL, n) = d1fC2C(2, 1:3, IBC_PERIODIC, IACCU_CD2) ! CD2/4 -> CD2
      d1rC2C(2, 1:3, IBC_INTRPL, n) = d1rC2C(2, 1:3, IBC_PERIODIC, IACCU_CD2) ! CD2/4 -> CD2
      d1fC2C(4, 1:3, IBC_INTRPL, n) = d1fC2C(4, 1:3, IBC_PERIODIC, IACCU_CD2) ! CD2/4 -> CD2
      d1rC2C(4, 1:3, IBC_INTRPL, n) = d1rC2C(4, 1:3, IBC_PERIODIC, IACCU_CD2) ! CD2/4 -> CD2
    end do

    do n = 3, 4
      !d1fC2C(2, 1:3, IBC_INTRPL, n) = d1fC2C(1, 1:3, IBC_INTRPL, n) ! exterpolation, following Line 1
      !d1rC2C(2, 1:3, IBC_INTRPL, n) = d1rC2C(1, 1:3, IBC_INTRPL, n) ! exterpolation, following Line 1
      !d1fC2C(4, 1:3, IBC_INTRPL, n) = d1fC2C(5, 1:3, IBC_INTRPL, n) ! exterpolation, following Line 1
      !d1rC2C(4, 1:3, IBC_INTRPL, n) = d1rC2C(5, 1:3, IBC_INTRPL, n) ! exterpolation, following Line 1

      d1fC2C(2, 1:3, IBC_INTRPL, n) = d1fC2C(2, 1:3, IBC_PERIODIC, IACCU_CP4) ! CP4/6 -> CP4
      d1rC2C(2, 1:3, IBC_INTRPL, n) = d1rC2C(2, 1:3, IBC_PERIODIC, IACCU_CP4) ! CP4/6 -> CP4
      d1fC2C(4, 1:3, IBC_INTRPL, n) = d1fC2C(2, 1:3, IBC_PERIODIC, IACCU_CP4) ! CP4/6 -> CP4
      d1rC2C(4, 1:3, IBC_INTRPL, n) = d1rC2C(2, 1:3, IBC_PERIODIC, IACCU_CP4) ! CP4/6 -> CP4
    end do
!==========================================================================================================
! 1st-derivative, P2P :
! d1fP2P : "d1"=first deriviative, "f"=f'  side, "P2P"= point(node) 2 point 
! d1rP2P : "d1"=first deriviative, "r"=rhs side, "P2P"= point(node) 2 point 
! alpha * f'_{i'-1} + f'_i' + alpha * f'_{i'+1} = a/(2h) * ( f_{i'+1} - f_{i'-1} ) + &
!                                                 b/(4h) * ( f_{i'+2} - f_{i'-2} )
! P2P unknows: 
! when i'=1,     need: LHS: f'_0';      RHS: f_0', f_{-1'}
! when i'=2,     need:                  RHS: f_0'
! when i'=np'-1, need:                  RHS: f_{np'+1}
! when i'=np',   need: LHS: f'_{np'+1}; RHS: f_{np'+1}, f_{np'+2}
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, P2P, IBC_PERIODIC, unknowns from both rhs and lhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    d1fP2P(:, :, IBC_PERIODIC,  :) = d1fC2C(:, :, IBC_PERIODIC,  :)
    d1rP2P(:, :, IBC_PERIODIC,  :) = d1rC2C(:, :, IBC_PERIODIC,  :)
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, P2P : IBC_SYMMETRIC, unknowns from both rhs and lhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      d1fP2P(1,   1, IBC_SYMMETRIC, n) = ZERO ! not used
      d1fP2P(1,   2, IBC_SYMMETRIC, n) = ONE
      d1fP2P(1,   3, IBC_SYMMETRIC, n) = ZERO
      d1fP2P(5,   1, IBC_SYMMETRIC, n) = d1fP2P(1,   3, IBC_SYMMETRIC, n)
      d1fP2P(5,   2, IBC_SYMMETRIC, n) = d1fP2P(1,   2, IBC_SYMMETRIC, n)
      d1fP2P(5,   3, IBC_SYMMETRIC, n) = d1fP2P(1,   1, IBC_SYMMETRIC, n)
      d1fP2P(2:4, :, IBC_SYMMETRIC, n) = d1fP2P(2:4, :, IBC_PERIODIC,  n)
      d1rP2P(:,   :, IBC_SYMMETRIC, n) = d1rP2P(:,   :, IBC_PERIODIC,  n)
    end do 
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, P2P : IBC_ASYMMETRIC, unknowns from both rhs and lhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      d1fP2P(1,   1, IBC_ASYMMETRIC, n) = ZERO ! not used
      d1fP2P(1,   2, IBC_ASYMMETRIC, n) = ONE
      d1fP2P(1,   3, IBC_ASYMMETRIC, n) = TWO * alpha(n)
      d1fP2P(5,   1, IBC_ASYMMETRIC, n) = d1fP2P(1,   3, IBC_ASYMMETRIC, n)
      d1fP2P(5,   2, IBC_ASYMMETRIC, n) = d1fP2P(1,   2, IBC_ASYMMETRIC, n)
      d1fP2P(5,   3, IBC_ASYMMETRIC, n) = d1fP2P(1,   1, IBC_ASYMMETRIC, n)
      d1fP2P(2:4, :, IBC_ASYMMETRIC, n) = d1fP2P(2:4, :, IBC_PERIODIC  , n)
      d1rP2P(:,   :, IBC_ASYMMETRIC, n) = d1rP2P(:,   :, IBC_PERIODIC  , n)
    end do
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, P2P : IBC_INTERIOR, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    d1fP2P(:, :, IBC_INTERIOR,  :) = d1fC2C(:, :, IBC_INTERIOR,  :)
    d1rP2P(:, :, IBC_INTERIOR,  :) = d1rC2C(:, :, IBC_INTERIOR,  :)
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, P2P : IBC_DIRICHLET, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    d1fP2P(:, :, IBC_DIRICHLET, :) = d1fC2C(:, :, IBC_DIRICHLET, :)
    d1rP2P(:, :, IBC_DIRICHLET, :) = d1rC2C(:, :, IBC_DIRICHLET, :)
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, P2P : NEUMANN, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      d1fP2P(1, 1,   IBC_NEUMANN, n) = ZERO ! not used
      d1fP2P(1, 2,   IBC_NEUMANN, n) = ONE
      d1fP2P(1, 3,   IBC_NEUMANN, n) = ZERO
      d1rP2P(1, 1:3, IBC_NEUMANN, n) = ZERO

      d1fP2P(5, 1,   IBC_NEUMANN, n) = ZERO
      d1fP2P(5, 2,   IBC_NEUMANN, n) = ONE
      d1fP2P(5, 3,   IBC_NEUMANN, n) = ZERO
      d1rP2P(5, 1:3, IBC_NEUMANN, n) = ZERO

      d1fP2P(2:4, 1:3, IBC_NEUMANN, n) = d1fP2P(2:4, 1:3, IBC_INTERIOR, n)
      d1rP2P(2:4, 1:3, IBC_NEUMANN, n) = d1rP2P(2:4, 1:3, IBC_INTERIOR, n)
    end do
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, P2P : exterpolation
! alpha * f'_{i'-1} + f'_i' + alpha * f'_{i'+1} = a/(2h) * ( f_{i'+1} - f_{i'-1} ) + &
!                                                 b/(4h) * ( f_{i'+2} - f_{i'-2} )
! 1st-derivative, C2C, IBC_INTRPL, no bc, no reconstuction. exterpolation only. 
!                     f'_1' + alpha * f'_{2'}   = a1 * f_1' + b1 * f_2' + c1 * f_3'
! alpha * f'_{1'}   + f'_2' + alpha * f'_{3'}   = a/(2h) * ( f_{i'+1} - f_{i'-1} ) 
! alpha * f'_{i'-1} + f'_i' + alpha * f'_{i'+1} = a/(2h) * ( f_{i'+1} - f_{i'-1} ) + &
!                                                 b/(4h) * ( f_{i'+2} - f_{i'-2} )
!----------------------------------------------------------------------------------------------------------
    d1fP2P(:, :, IBC_INTRPL,    :) = d1fC2C(:, :, IBC_INTRPL, :)
    d1rP2P(:, :, IBC_INTRPL,    :) = d1rC2C(:, :, IBC_INTRPL, :)
!==========================================================================================================
! 1st derivative on staggered grids C2P
! alpha * f'_{i'-1} + f'_i' + alpha * f'_{i'+1} = a/(h ) * ( f_{i}   - f_{i-1} ) + &
!                                                 b/(3h) * ( f_{i+1} - f_{i-2} )
! when i' = 1',    need: f'_0', f_0, f_{-1}
! when i' = 2',    need: f_0
! when i' = np-1', need: f_{np}
! when i' = np',   need: f'_{np+1'}, f_{np}, f_{np+1}
!==========================================================================================================
    alpha(IACCU_CD2) = ZERO
        a(IACCU_CD2) = ONE
        b(IACCU_CD2) = ZERO

    alpha(IACCU_CD4) = ZERO
        a(IACCU_CD4) = NINE * EIGHTH
        b(IACCU_CD4) = -ONE * EIGHTH

    alpha(IACCU_CP4) = ONE / TWENTYTWO
        a(IACCU_CP4) = TWELVE / ELEVEN
        b(IACCU_CP4) = ZERO

    alpha(IACCU_CP6) = NINE / SIXTYTWO
        a(IACCU_CP6) = SIXTYTHREE / SIXTYTWO
        b(IACCU_CP6) = SEVENTEEN / SIXTYTWO
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2P, IBC_PERIODIC, unknowns from both rhs and lhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      d1fC2P(1:5, 1, IBC_PERIODIC, n) = alpha(n)
      d1fC2P(1:5, 2, IBC_PERIODIC, n) = ONE
      d1fC2P(1:5, 3, IBC_PERIODIC, n) = alpha(n)
      d1rC2P(1:5, 1, IBC_PERIODIC, n) = a(n)             ! a
      d1rC2P(1:5, 2, IBC_PERIODIC, n) = b(n) * ONE_THIRD ! b/3
      d1rC2P(1:5, 3, IBC_PERIODIC, n) = ZERO             ! not used.
    end do
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2P, IBC_SYMMETRIC, unknowns from both rhs and lhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      d1fC2P(1,   1, IBC_SYMMETRIC, n) = ZERO ! not used
      d1fC2P(1,   2, IBC_SYMMETRIC, n) = ONE
      d1fC2P(1,   3, IBC_SYMMETRIC, n) = ZERO
      d1fC2P(5,   1, IBC_SYMMETRIC, n) = d1fC2P(1,   3, IBC_SYMMETRIC, n)
      d1fC2P(5,   2, IBC_SYMMETRIC, n) = d1fC2P(1,   2, IBC_SYMMETRIC, n)
      d1fC2P(5,   3, IBC_SYMMETRIC, n) = d1fC2P(1,   1, IBC_SYMMETRIC, n)
      d1fC2P(2:4, :, IBC_SYMMETRIC, n) = d1fC2P(2:4, :, IBC_PERIODIC,  n)
      d1rC2P(:,   :, IBC_SYMMETRIC, n) = d1rC2P(:,   :, IBC_PERIODIC,  n)
    end do
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2P : IBC_ASYMMETRIC, unknowns from both rhs and lhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      d1fC2P(1,   1, IBC_ASYMMETRIC, n) = ZERO ! not used
      d1fC2P(1,   2, IBC_ASYMMETRIC, n) = ONE
      d1fC2P(1,   3, IBC_ASYMMETRIC, n) = TWO * alpha(n)
      d1fC2P(5,   1, IBC_ASYMMETRIC, n) = d1fC2P(1,   3, IBC_ASYMMETRIC, n)
      d1fC2P(5,   2, IBC_ASYMMETRIC, n) = d1fC2P(1,   2, IBC_ASYMMETRIC, n)
      d1fC2P(5,   3, IBC_ASYMMETRIC, n) = d1fC2P(1,   1, IBC_ASYMMETRIC, n)
      d1fC2P(2:4, :, IBC_ASYMMETRIC, n) = d1fC2P(2:4, :, IBC_PERIODIC,   n)
      d1rC2P(:,   :, IBC_ASYMMETRIC, n) = d1rC2P(:,   :, IBC_PERIODIC,   n)
    end do
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2P, IBC_INTERIOR, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    d1fC2P(:, :, IBC_INTERIOR, :) = d1fC2P(:, :, IBC_PERIODIC, :)
    d1rC2P(:, :, IBC_INTERIOR, :) = d1rC2P(:, :, IBC_PERIODIC, :)

    d1fC2P(1:2, :, IBC_INTERIOR, IACCU_CP4) = d1fC2P(1:2, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    d1rC2P(1:2, :, IBC_INTERIOR, IACCU_CP4) = d1rC2P(1:2, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    d1fC2P(4:5, :, IBC_INTERIOR, IACCU_CP4) = d1fC2P(4:5, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    d1rC2P(4:5, :, IBC_INTERIOR, IACCU_CP4) = d1rC2P(4:5, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD

    d1fC2P(1:2, :, IBC_INTERIOR, IACCU_CP6) = d1fC2P(1:2, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    d1rC2P(1:2, :, IBC_INTERIOR, IACCU_CP6) = d1rC2P(1:2, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    d1fC2P(4:5, :, IBC_INTERIOR, IACCU_CP6) = d1fC2P(4:5, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    d1rC2P(4:5, :, IBC_INTERIOR, IACCU_CP6) = d1rC2P(4:5, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2P, IBC_DIRICHLET, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    d1fC2P(:, :, IBC_DIRICHLET, :) = d1fC2P(:, :, IBC_INTERIOR, :)
    d1rC2P(:, :, IBC_DIRICHLET, :) = d1rC2P(:, :, IBC_INTERIOR, :)
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, C2P : IBC_NEUMANN, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      d1fC2P(1, 1,   IBC_NEUMANN, n) = ZERO ! not used
      d1fC2P(1, 2,   IBC_NEUMANN, n) = ONE
      d1fC2P(1, 3,   IBC_NEUMANN, n) = ZERO
      d1rC2P(1, 1:3, IBC_NEUMANN, n) = ZERO ! not used
      d1fC2P(5, 1,   IBC_NEUMANN, n) = ZERO
      d1fC2P(5, 2,   IBC_NEUMANN, n) = ONE
      d1fC2P(5, 3,   IBC_NEUMANN, n) = ZERO
      d1rC2P(5, 1:3, IBC_NEUMANN, n) = ZERO
      d1fC2P(2:4, :, IBC_NEUMANN, n) = d1fC2P(2:4, :, IBC_INTERIOR, n)
      d1rC2P(2:4, :, IBC_NEUMANN, n) = d1rC2P(2:4, :, IBC_INTERIOR, n)
    end do
!----------------------------------------------------------------------------------------------------------
! 1st-derivative : C2P, IBC_INTRPL, no bc, no reconstuction. exterpolation only. 
! alpha * f'_{i'-1} + f'_i' + alpha * f'_{i'+1} = a/(h ) * ( f_{i}   - f_{i-1} ) + &
!                                                 b/(3h) * ( f_{i+1} - f_{i-2} )
! when i' = 1',    need: f'_0', f_0, f_{-1}
! when i' = 2',    need: f_0
! when i' = np-1', need: f_{np}
! when i' = np',   need: f'_{np+1'}, f_{np}, f_{np+1}
! [ 1     alpha1                            ][f'_1']=[a1 * f_{1}/h  + b1 * f_{2}/h + c1 * f_{3}/h  ]
! [alpha2 1      alpha2                     ][f'_2'] [a2 * (f_{2} - f_{1})/h  ]
! [       alpha  1      alpha               ][f'_i'] [a *  (f_{i} - f_{i-1})/h + b/3 * (f_{i+1} - f_{i-2})/h]
! [                     alpha2 1      alpha2][f'_4'] [a2 * (f_{n-1} - f_{n-2})/h]
! [                            alpha1 1     ][f'_5'] [-a1 * f_{n-1}/h  - b1 * f_{n-2}/h - c1 * f_{n-3}/h]
! tested: low accuracy at the line 2 and 4.
!----------------------------------------------------------------------------------------------------------
    alpha1(IACCU_CD2) = ZERO
        a1(IACCU_CD2) = -ONE
        b1(IACCU_CD2) = ONE
        c1(IACCU_CD2) = ZERO

    alpha1(IACCU_CD4) = alpha1(IACCU_CD2)
        a1(IACCU_CD4) =     a1(IACCU_CD2)
        b1(IACCU_CD4) =     b1(IACCU_CD2)
        c1(IACCU_CD4) =     c1(IACCU_CD2)

    alpha1(IACCU_CP4) = alpha1(IACCU_CD2)!TWENTYTHREE
        a1(IACCU_CP4) =     a1(IACCU_CD2)!-TWENTYFIVE
        b1(IACCU_CP4) =     b1(IACCU_CD2)!TWENTYSIX
        c1(IACCU_CP4) =     c1(IACCU_CD2)!-ONE, need to be zero to keep compact scheme

    alpha1(IACCU_CP6) = alpha1(IACCU_CD2)
        a1(IACCU_CP6) =     a1(IACCU_CD2) 
        b1(IACCU_CP6) =     b1(IACCU_CD2) 
        c1(IACCU_CP6) =     c1(IACCU_CD2) 

    do n = 1, NACC
      d1fC2P(1, 1,   IBC_INTRPL, n) = ZERO ! not used
      d1fC2P(1, 2,   IBC_INTRPL, n) = ONE
      d1fC2P(1, 3,   IBC_INTRPL, n) = alpha1(n)
      d1rC2P(1, 1,   IBC_INTRPL, n) = a1(n)
      d1rC2P(1, 2,   IBC_INTRPL, n) = b1(n)
      d1rC2P(1, 3,   IBC_INTRPL, n) = c1(n)

      d1fC2P(5, 1,   IBC_INTRPL, n) =   d1fC2P(1, 3, IBC_INTRPL, n)
      d1fC2P(5, 2,   IBC_INTRPL, n) =   d1fC2P(1, 2, IBC_INTRPL, n)
      d1fC2P(5, 3,   IBC_INTRPL, n) =   d1fC2P(1, 1, IBC_INTRPL, n)
      d1rC2P(5, 1,   IBC_INTRPL, n) = - d1rC2P(1, 1, IBC_INTRPL, n)
      d1rC2P(5, 2,   IBC_INTRPL, n) = - d1rC2P(1, 2, IBC_INTRPL, n)
      d1rC2P(5, 3,   IBC_INTRPL, n) = - d1rC2P(1, 3, IBC_INTRPL, n)

      d1fC2P(3, 1:3, IBC_INTRPL, n) = d1fC2P(3, 1:3, IBC_PERIODIC, n)
      d1rC2P(3, 1:3, IBC_INTRPL, n) = d1rC2P(3, 1:3, IBC_PERIODIC, n)
    end do
    ! check which below method works good! 
    do n = 1, 2
      ! 
      ! method 1: exterpolation, following Line 1
      !d1fC2P(2, 1:3, IBC_INTRPL, n) = d1fC2P(1, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      !d1rC2P(2, 1:3, IBC_INTRPL, n) = d1rC2P(1, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      !d1fC2P(4, 1:3, IBC_INTRPL, n) = d1fC2P(5, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      !d1rC2P(4, 1:3, IBC_INTRPL, n) = d1rC2P(5, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      ! method 2: CD2/4 -> CD2
      d1fC2P(2, 1:3, IBC_INTRPL, n) = d1fC2P(2, 1:3, IBC_PERIODIC, IACCU_CD2) ! CD2/4 -> CD2
      d1rC2P(2, 1:3, IBC_INTRPL, n) = d1rC2P(2, 1:3, IBC_PERIODIC, IACCU_CD2) ! CD2/4 -> CD2
      d1fC2P(4, 1:3, IBC_INTRPL, n) = d1fC2P(4, 1:3, IBC_PERIODIC, IACCU_CD2) ! CD2/4 -> CD2
      d1rC2P(4, 1:3, IBC_INTRPL, n) = d1rC2P(4, 1:3, IBC_PERIODIC, IACCU_CD2) ! CD2/4 -> CD2
    end do

    do n = 3, 4
      ! method 1: exterpolation, following Line 1
      !d1fC2P(2, 1:3, IBC_INTRPL, n) = d1fC2P(1, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      !d1rC2P(2, 1:3, IBC_INTRPL, n) = d1rC2P(1, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      !d1fC2P(4, 1:3, IBC_INTRPL, n) = d1fC2P(5, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      !d1rC2P(4, 1:3, IBC_INTRPL, n) = d1rC2P(5, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      ! method 2: CP4/6 -> CP4
      d1fC2P(2, 1:3, IBC_INTRPL, n) = d1fC2P(2, 1:3, IBC_PERIODIC, IACCU_CP4) ! CP4/6 -> CP4
      d1rC2P(2, 1:3, IBC_INTRPL, n) = d1rC2P(2, 1:3, IBC_PERIODIC, IACCU_CP4) ! CP4/6 -> CP4
      d1fC2P(4, 1:3, IBC_INTRPL, n) = d1fC2P(2, 1:3, IBC_PERIODIC, IACCU_CP4) ! CP4/6 -> CP4
      d1rC2P(4, 1:3, IBC_INTRPL, n) = d1rC2P(2, 1:3, IBC_PERIODIC, IACCU_CP4) ! CP4/6 -> CP4
    end do
!==========================================================================================================
! 1st derivative on staggered grids P2C
! alpha * f'_{i-1} +  f'_i +  alpha * f'_{i+1}  = a/(h ) * ( f_{i'+1} - f_{i'} ) + &
!                                                 b/(3h) * ( f_{i'+2} - f_{i'-1} )
! when i = 1,    need: f'_0, f_0'
! when i = 2,    need: nothing
! when i = nc-1, need: nothing
! when i = nc,   need: f'_{nc+1'}, f_{np'}, f_{np'+1}
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, P2C, IBC_PERIODIC, unknowns from both rhs and lhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    d1fP2C(:, :, IBC_PERIODIC, :) = d1fC2P(:, :, IBC_PERIODIC, :)
    d1rP2C(:, :, IBC_PERIODIC, :) = d1rC2P(:, :, IBC_PERIODIC, :)
!----------------------------------------------------------------------------------------------------------
! 1st-derivative : P2C : IBC_SYMMETRIC, unknowns from both rhs and lhs could be reconstructed from bc.
! alpha * f'_{i-1} +  f'_i +  alpha * f'_{i+1}  = a/(h ) * ( f_{i'+1} - f_{i'} ) + &
!                                                 b/(3h) * ( f_{i'+2} - f_{i'-1} )
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      d1fP2C(1,   1, IBC_SYMMETRIC, n) = ZERO ! not used
      d1fP2C(1,   2, IBC_SYMMETRIC, n) = ONE - alpha(n)
      d1fP2C(1,   3, IBC_SYMMETRIC, n) = alpha(n)
      d1fP2C(5,   1, IBC_SYMMETRIC, n) = d1fP2C(1,   3, IBC_SYMMETRIC, n)
      d1fP2C(5,   2, IBC_SYMMETRIC, n) = d1fP2C(1,   2, IBC_SYMMETRIC, n)
      d1fP2C(5,   3, IBC_SYMMETRIC, n) = d1fP2C(1,   1, IBC_SYMMETRIC, n)
      d1fP2C(2:4, :, IBC_SYMMETRIC, n) = d1fP2C(2:4, :, IBC_PERIODIC , n)
      d1rP2C(:,   :, IBC_SYMMETRIC, n) = d1rP2C(:,   :, IBC_PERIODIC , n)
    end do
!----------------------------------------------------------------------------------------------------------
! 1st-derivative : IBC_ASYMMETRIC, unknowns from both rhs and lhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      d1fP2C(1,   1, IBC_ASYMMETRIC, n) = ZERO ! not used
      d1fP2C(1,   2, IBC_ASYMMETRIC, n) = ONE + alpha(n)
      d1fP2C(1,   3, IBC_ASYMMETRIC, n) = alpha(n)
      d1fP2C(5,   1, IBC_ASYMMETRIC, n) = d1fP2C(1,   3, IBC_ASYMMETRIC, n)
      d1fP2C(5,   2, IBC_ASYMMETRIC, n) = d1fP2C(1,   2, IBC_ASYMMETRIC, n)
      d1fP2C(5,   3, IBC_ASYMMETRIC, n) = d1fP2C(1,   1, IBC_ASYMMETRIC, n)
      d1fP2C(2:4, :, IBC_ASYMMETRIC, n) = d1fP2C(2:4, :, IBC_PERIODIC  , n)
      d1rP2C(:,   :, IBC_ASYMMETRIC, n) = d1rP2C(:,   :, IBC_PERIODIC  , n)
    end do
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, P2C, IBC_INTERIOR, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    d1fP2C(:, :, IBC_INTERIOR, :) = d1fP2C(:, :, IBC_PERIODIC, :)
    d1rP2C(:, :, IBC_INTERIOR, :) = d1rP2C(:, :, IBC_PERIODIC, :)

    d1fP2C(1, :, IBC_INTERIOR, IACCU_CP4) = d1fP2C(1, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    d1rP2C(1, :, IBC_INTERIOR, IACCU_CP4) = d1rP2C(1, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    d1fP2C(5, :, IBC_INTERIOR, IACCU_CP4) = d1fP2C(5, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    d1rP2C(5, :, IBC_INTERIOR, IACCU_CP4) = d1rP2C(5, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD

    d1fP2C(1, :, IBC_INTERIOR, IACCU_CP6) = d1fP2C(1, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    d1rP2C(1, :, IBC_INTERIOR, IACCU_CP6) = d1rP2C(1, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    d1fP2C(5, :, IBC_INTERIOR, IACCU_CP6) = d1fP2C(4, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    d1rP2C(5, :, IBC_INTERIOR, IACCU_CP6) = d1rP2C(4, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, P2C, Dirichlet, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    d1fP2C(:, :, IBC_DIRICHLET, :) = d1fP2C(:, :, IBC_INTERIOR, :)
    d1rP2C(:, :, IBC_DIRICHLET, :) = d1rP2C(:, :, IBC_INTERIOR, :)
!----------------------------------------------------------------------------------------------------------
! 1st-derivative, P2C, Neumann, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    d1fP2C(:, :, IBC_NEUMANN, :) = d1fP2C(:, :, IBC_INTERIOR, :)
    d1rP2C(:, :, IBC_NEUMANN, :) = d1rP2C(:, :, IBC_INTERIOR, :)

!----------------------------------------------------------------------------------------------------------
! 1st-derivative : ! P2C : exterpolation, no bc, no reconstuction. exterpolation only. 
! [ 1     alpha1                            ][f'_1]=[a1 * f_{1'}/h  + b1 * f_{2'}/h + c1 * f_{3'}/h  ]
! [alpha2 1      alpha2                     ][f'_2] [a2 * (f_{3'} - f_{2'})/h  ]
! [       alpha  1      alpha               ][f'_i] [a *  (f_{i'+1} - f_{i'})/h + b/3 * (f_{i'+2} - f_{i'-1})/h]
! [                     alpha2 1      alpha2][f'_4] [a2 * (f_{n'} - f_{n'-1})/h]
! [                            alpha1 1     ][f'_5] [-a1 * f_{n'+1}/h  - b1 * f_{n'}/h - c1 * f_{n'-1}/h]
!----------------------------------------------------------------------------------------------------------
    alpha1(IACCU_CD2) = ZERO
        a1(IACCU_CD2) = -ONE
        b1(IACCU_CD2) = ONE
        c1(IACCU_CD2) = ZERO
    alpha1(IACCU_CD4) = alpha1(IACCU_CD2)
        a1(IACCU_CD4) =     a1(IACCU_CD2)
        b1(IACCU_CD4) =     b1(IACCU_CD2)
        c1(IACCU_CD4) =     c1(IACCU_CD2)
    alpha1(IACCU_CP4) = alpha1(IACCU_CD2) !ZERO
        a1(IACCU_CP4) =     a1(IACCU_CD2) !-ONE
        b1(IACCU_CP4) =     b1(IACCU_CD2) !ONE
        c1(IACCU_CP4) =     c1(IACCU_CD2) !ZERO
    alpha1(IACCU_CP6) = -ONE
        a1(IACCU_CP6) = -ONE
        b1(IACCU_CP6) = TWO
        c1(IACCU_CP6) = -ONE

    do n = 1, NACC
      d1fP2C(1, 1,   IBC_INTRPL, n) = ZERO ! not used
      d1fP2C(1, 2,   IBC_INTRPL, n) = ONE
      d1fP2C(1, 3,   IBC_INTRPL, n) = alpha1(n)
      d1rP2C(1, 1,   IBC_INTRPL, n) = a1(n)
      d1rP2C(1, 2,   IBC_INTRPL, n) = b1(n)
      d1rP2C(1, 3,   IBC_INTRPL, n) = c1(n)

      d1fP2C(5, 1,   IBC_INTRPL, n) =   d1fP2C(1, 3, IBC_INTRPL, n)
      d1fP2C(5, 2,   IBC_INTRPL, n) =   d1fP2C(1, 2, IBC_INTRPL, n)
      d1fP2C(5, 3,   IBC_INTRPL, n) =   d1fP2C(1, 1, IBC_INTRPL, n)
      d1rP2C(5, 1,   IBC_INTRPL, n) = - d1rP2C(1, 1, IBC_INTRPL, n)
      d1rP2C(5, 2,   IBC_INTRPL, n) = - d1rP2C(1, 2, IBC_INTRPL, n)
      d1rP2C(5, 3,   IBC_INTRPL, n) = - d1rP2C(1, 3, IBC_INTRPL, n)

      d1fP2C(2:4, 1:3, IBC_INTRPL, n) = d1fP2C(2:4, 1:3, IBC_PERIODIC, n)
      d1rP2C(2:4, 1:3, IBC_INTRPL, n) = d1rP2C(2:4, 1:3, IBC_PERIODIC, n)
    end do
!==========================================================================================================
!interpolation. C2P 
! alpha * f_{i'-1} + f_i' + alpha * f_{i'+1} = a/2 * ( f_{i}   + f_{i-1} ) + &
!                                              b/2 * ( f_{i+1} + f_{i-2} )
! when i' = 1,    need: f_{0'}, f_{0}, f_{-1}
! when i' = 2,    need:         f_{0}
! when i' = np-1, need:         f_{np'}
! when i' = np,   need: f_{np'+1}, f_{np'}, f_{np'+1}
!==========================================================================================================
      alpha(IACCU_CD2) = ZERO
          a(IACCU_CD2) = ONE
          b(IACCU_CD2) = ZERO
      alpha(IACCU_CD4) = ZERO
          a(IACCU_CD4) = NINE * EIGHTH
          b(IACCU_CD4) = -ONE * EIGHTH
      alpha(IACCU_CP4) = ONE / SIX
          a(IACCU_CP4) = FOUR * ONE_THIRD
          b(IACCU_CP4) = ZERO
      alpha(IACCU_CP6) = THREE * ZPONE
          a(IACCU_CP6) = ONEPFIVE
          b(IACCU_CP6) = ONE * ZPONE
!----------------------------------------------------------------------------------------------------------
!interpolation : C2P for IBC_PERIODIC, unknowns from both lhs and rhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      m1fC2P(1:5, 1, IBC_PERIODIC, n) = alpha(n)
      m1fC2P(1:5, 2, IBC_PERIODIC, n) = ONE
      m1fC2P(1:5, 3, IBC_PERIODIC, n) = alpha(n)
      m1rC2P(1:5, 1, IBC_PERIODIC, n) = a(n) * HALF
      m1rC2P(1:5, 2, IBC_PERIODIC, n) = b(n) * HALF
      m1rC2P(1:5, 3, IBC_PERIODIC, n) = ZERO ! not used
    end do 
!----------------------------------------------------------------------------------------------------------
!interpolation. C2P for IBC_SYMMETRIC, unknowns from both lhs and rhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      m1fC2P(1,   1, IBC_SYMMETRIC, n) = ZERO ! not used
      m1fC2P(1,   2, IBC_SYMMETRIC, n) = ONE
      m1fC2P(1,   3, IBC_SYMMETRIC, n) = alpha(n) + alpha(n)
      m1fC2P(5,   1, IBC_SYMMETRIC, n) = m1fC2P(1,   3, IBC_SYMMETRIC, n)
      m1fC2P(5,   2, IBC_SYMMETRIC, n) = m1fC2P(1,   2, IBC_SYMMETRIC, n)
      m1fC2P(5,   3, IBC_SYMMETRIC, n) = m1fC2P(1,   1, IBC_SYMMETRIC, n)
      m1fC2P(2:4, :, IBC_SYMMETRIC, n) = m1fC2P(2:4, :, IBC_PERIODIC , n)
      m1rC2P(:,   :, IBC_SYMMETRIC, n) = m1rC2P(:,   :, IBC_PERIODIC , n)
    end do
!----------------------------------------------------------------------------------------------------------
!interpolation. C2P for IBC_ASYMMETRIC, unknowns from both lhs and rhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      m1fC2P(1,   1, IBC_ASYMMETRIC, n) = ZERO ! not used
      m1fC2P(1,   2, IBC_ASYMMETRIC, n) = ONE
      m1fC2P(1,   3, IBC_ASYMMETRIC, n) = ZERO
      m1fC2P(5,   1, IBC_ASYMMETRIC, n) = m1fC2P(1,   3, IBC_ASYMMETRIC, n)
      m1fC2P(5,   2, IBC_ASYMMETRIC, n) = m1fC2P(1,   2, IBC_ASYMMETRIC, n)
      m1fC2P(5,   3, IBC_ASYMMETRIC, n) = m1fC2P(1,   1, IBC_ASYMMETRIC, n)
      m1fC2P(2:4, :, IBC_ASYMMETRIC, n) = m1fC2P(2:4, :, IBC_PERIODIC  , n)
      m1rC2P(:,   :, IBC_ASYMMETRIC, n) = m1rC2P(:,   :, IBC_PERIODIC  , n)
    end do
!----------------------------------------------------------------------------------------------------------
! interpolation : C2P for IBC_INTERIOR, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    m1fC2P(:,   :, IBC_INTERIOR, :) = m1fC2P(:,   :, IBC_PERIODIC, :)
    m1rC2P(:,   :, IBC_INTERIOR, :) = m1rC2P(:,   :, IBC_PERIODIC, :)

    m1fC2P(1:2, :, IBC_INTERIOR, IACCU_CP4) = m1fC2P(1:2, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    m1rC2P(1:2, :, IBC_INTERIOR, IACCU_CP4) = m1rC2P(1:2, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    m1fC2P(4:5, :, IBC_INTERIOR, IACCU_CP4) = m1fC2P(4:5, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    m1rC2P(4:5, :, IBC_INTERIOR, IACCU_CP4) = m1rC2P(4:5, :, IBC_PERIODIC, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD

    m1fC2P(1:2, :, IBC_INTERIOR, IACCU_CP6) = m1fC2P(1:2, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    m1rC2P(1:2, :, IBC_INTERIOR, IACCU_CP6) = m1rC2P(1:2, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    m1fC2P(4:5, :, IBC_INTERIOR, IACCU_CP6) = m1fC2P(4:5, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    m1rC2P(4:5, :, IBC_INTERIOR, IACCU_CP6) = m1rC2P(4:5, :, IBC_PERIODIC, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
!----------------------------------------------------------------------------------------------------------
! interpolation : C2P for IBC_NEUMANN, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    m1fC2P(:, :, IBC_NEUMANN, :) = m1fC2P(:, :, IBC_INTERIOR, :)
    m1rC2P(:, :, IBC_NEUMANN, :) = m1rC2P(:, :, IBC_INTERIOR, :)
!----------------------------------------------------------------------------------------------------------
! interpolation : C2P for IBC_DIRICHLET, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      m1fC2P(1, 1,   IBC_DIRICHLET, n) = ZERO ! not used
      m1fC2P(1, 2,   IBC_DIRICHLET, n) = ONE
      m1fC2P(1, 3,   IBC_DIRICHLET, n) = ZERO
      m1rC2P(1, 1:3, IBC_DIRICHLET, n) = ZERO ! not used
      m1fC2P(5, 1,   IBC_DIRICHLET, n) = ZERO
      m1fC2P(5, 2,   IBC_DIRICHLET, n) = ONE
      m1fC2P(5, 3,   IBC_DIRICHLET, n) = ZERO
      m1rC2P(5, 1:3, IBC_DIRICHLET, n) = ZERO
      m1fC2P(2:4, :, IBC_DIRICHLET, n) = m1fC2P(2:4, :, IBC_INTERIOR, n)
      m1rC2P(2:4, :, IBC_DIRICHLET, n) = m1rC2P(2:4, :, IBC_INTERIOR, n)
    end do
!----------------------------------------------------------------------------------------------------------
!interpolation. C2P, exterpolation
!-----------------------------------------------------------------------------
    alpha1(IACCU_CD2) = ZERO
        a1(IACCU_CD2) = ONEPFIVE
        b1(IACCU_CD2) = -HALF
        c1(IACCU_CD2) = ZERO
    alpha1(IACCU_CD4) = ZERO
        a1(IACCU_CD4) = FIFTEEN * EIGHTH
        b1(IACCU_CD4) = - FIVE * QUARTER
        c1(IACCU_CD4) = THREE * EIGHTH
    alpha1(IACCU_CP4) = THREE
        a1(IACCU_CP4) = THREE
        b1(IACCU_CP4) = ONE
        c1(IACCU_CP4) = ZERO
    alpha1(IACCU_CP6) = FIVE 
        a1(IACCU_CP6) = FIFTEEN * QUARTER
        b1(IACCU_CP6) = TWOPFIVE
        c1(IACCU_CP6) = -QUARTER

    do n = 1, NACC
      m1fC2P(1, 1,   IBC_INTRPL, n) = ZERO ! not used
      m1fC2P(1, 2,   IBC_INTRPL, n) = ONE
      m1fC2P(1, 3,   IBC_INTRPL, n) = alpha1(n)
      m1rC2P(1, 1,   IBC_INTRPL, n) = a1(n)
      m1rC2P(1, 2,   IBC_INTRPL, n) = b1(n)
      m1rC2P(1, 3,   IBC_INTRPL, n) = c1(n)

      m1fC2P(5, 1,   IBC_INTRPL, n) = m1fC2P(1, 3, IBC_INTRPL, n)
      m1fC2P(5, 2,   IBC_INTRPL, n) = m1fC2P(1, 2, IBC_INTRPL, n)
      m1fC2P(5, 3,   IBC_INTRPL, n) = m1fC2P(1, 1, IBC_INTRPL, n)
      m1rC2P(5, 1,   IBC_INTRPL, n) = m1rC2P(1, 1, IBC_INTRPL, n)
      m1rC2P(5, 2,   IBC_INTRPL, n) = m1rC2P(1, 2, IBC_INTRPL, n)
      m1rC2P(5, 3,   IBC_INTRPL, n) = m1rC2P(1, 3, IBC_INTRPL, n)

      m1fC2P(3, 1:3, IBC_INTRPL, n) = m1fC2P(3, 1:3, IBC_PERIODIC, n)
      m1rC2P(3, 1:3, IBC_INTRPL, n) = m1rC2P(3, 1:3, IBC_PERIODIC, n)
    end do
    ! check which below method works good! 
    do n = 1, 2
      ! method 1: exterpolation, following Line 1
      !m1fC2P(2, 1:3, IBC_INTRPL, n) = m1fC2P(1, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      !m1rC2P(2, 1:3, IBC_INTRPL, n) = m1rC2P(1, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      !m1fC2P(4, 1:3, IBC_INTRPL, n) = m1fC2P(5, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      !m1rC2P(4, 1:3, IBC_INTRPL, n) = m1rC2P(5, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      ! method 2: CD2/4 -> CD2
      m1fC2P(2, 1:3, IBC_INTRPL, n) = m1fC2P(2, 1:3, IBC_PERIODIC, IACCU_CD2) ! CD2/4 -> CD2
      m1rC2P(2, 1:3, IBC_INTRPL, n) = m1rC2P(2, 1:3, IBC_PERIODIC, IACCU_CD2) ! CD2/4 -> CD2
      m1fC2P(4, 1:3, IBC_INTRPL, n) = m1fC2P(4, 1:3, IBC_PERIODIC, IACCU_CD2) ! CD2/4 -> CD2
      m1rC2P(4, 1:3, IBC_INTRPL, n) = m1rC2P(4, 1:3, IBC_PERIODIC, IACCU_CD2) ! CD2/4 -> CD2
    end do

    do n = 3, 4
      ! method 1: exterpolation, following Line 1
      !m1fC2P(2, 1:3, IBC_INTRPL, n) = m1fC2P(1, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      !m1rC2P(2, 1:3, IBC_INTRPL, n) = m1rC2P(1, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      !m1fC2P(4, 1:3, IBC_INTRPL, n) = m1fC2P(5, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      !m1rC2P(4, 1:3, IBC_INTRPL, n) = m1rC2P(5, 1:3, IBC_INTRPL,   n) ! exterpolation, following Line 1
      ! method 2: CP4/6 -> CP4
      m1fC2P(2, 1:3, IBC_INTRPL, n) = m1fC2P(2, 1:3, IBC_PERIODIC, IACCU_CP4) ! CP4/6 -> CP4
      m1rC2P(2, 1:3, IBC_INTRPL, n) = m1rC2P(2, 1:3, IBC_PERIODIC, IACCU_CP4) ! CP4/6 -> CP4
      m1fC2P(4, 1:3, IBC_INTRPL, n) = m1fC2P(2, 1:3, IBC_PERIODIC, IACCU_CP4) ! CP4/6 -> CP4
      m1rC2P(4, 1:3, IBC_INTRPL, n) = m1rC2P(2, 1:3, IBC_PERIODIC, IACCU_CP4) ! CP4/6 -> CP4
    end do
!==========================================================================================================
!interpolation. P2C 
! P2C : i_max = nc
! alpha * f_{i-1} + f_i + alpha * f_{i+1} =    a/2 * ( f_{i'}   + f_{i'+1} ) + &
!                                              b/2 * ( f_{i'+2} + f_{i'-1} )
! when i = 1,    need: LHS: f_{0}, RHS: f_{0'}
! when i = 2,    need: nothing
! when i = nc-1, need: nothing
! when i = nc,   need: LHS: f_{np}, RHS: f_{np'+1}
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
!interpolation : P2C for IBC_PERIODIC, unknowns from both lhs and rhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    m1fP2C(:, :, IBC_PERIODIC, :) = m1fC2P(:, :, IBC_PERIODIC, :)
    m1rP2C(:, :, IBC_PERIODIC, :) = m1rC2P(:, :, IBC_PERIODIC, :)
!----------------------------------------------------------------------------------------------------------
!interpolation. P2C. IBC_SYMMETRIC, unknowns from both lhs and rhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      m1fP2C(1,   1, IBC_SYMMETRIC, n) = ZERO ! not used
      m1fP2C(1,   2, IBC_SYMMETRIC, n) = ONE + alpha(n)
      m1fP2C(1,   3, IBC_SYMMETRIC, n) = alpha(n)
      m1fP2C(5,   1, IBC_SYMMETRIC, n) = m1fP2C(1,   3, IBC_SYMMETRIC, n)
      m1fP2C(5,   2, IBC_SYMMETRIC, n) = m1fP2C(1,   2, IBC_SYMMETRIC, n)
      m1fP2C(5,   3, IBC_SYMMETRIC, n) = m1fP2C(1,   1, IBC_SYMMETRIC, n)
      m1fP2C(2:4, :, IBC_SYMMETRIC, n) = m1fP2C(2:4, :, IBC_PERIODIC , n)
      m1rP2C(:,   :, IBC_SYMMETRIC, n) = m1rP2C(:,   :, IBC_PERIODIC , n)
    end do
!----------------------------------------------------------------------------------------------------------
!interpolation. P2C. IBC_ASYMMETRIC, unknowns from both lhs and rhs could be reconstructed from bc.
!----------------------------------------------------------------------------------------------------------
    do n = 1, NACC
      m1fP2C(1,   1, IBC_ASYMMETRIC, n) = ZERO ! not used
      m1fP2C(1,   2, IBC_ASYMMETRIC, n) = ONE - alpha(n)
      m1fP2C(1,   3, IBC_ASYMMETRIC, n) = alpha(n)
      m1fP2C(5,   1, IBC_ASYMMETRIC, n) = m1fP2C(1,   3, IBC_ASYMMETRIC, n)
      m1fP2C(5,   2, IBC_ASYMMETRIC, n) = m1fP2C(1,   2, IBC_ASYMMETRIC, n)
      m1fP2C(5,   3, IBC_ASYMMETRIC, n) = m1fP2C(1,   1, IBC_ASYMMETRIC, n)
      m1fP2C(2:4, :, IBC_ASYMMETRIC, n) = m1fP2C(2:4, :, IBC_PERIODIC  , n)
      m1rP2C(:,   :, IBC_ASYMMETRIC, n) = m1rP2C(:,   :, IBC_PERIODIC  , n)
    end do
!----------------------------------------------------------------------------------------------------------
! interpolation : P2C, IBC_INTERIOR, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    m1fP2C(1:NL, 1:NS, IBC_INTERIOR, 1:NACC) = m1fP2C(1:NL, 1:NS, IBC_PERIODIC, 1:NACC)
    m1rP2C(1:NL, 1:NS, IBC_INTERIOR, 1:NACC) = m1rP2C(1:NL, 1:NS, IBC_PERIODIC, 1:NACC)

    m1fP2C(1, :, IBC_INTERIOR, IACCU_CP4) = m1fP2C(1, :, IBC_INTERIOR, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    m1rP2C(1, :, IBC_INTERIOR, IACCU_CP4) = m1rP2C(1, :, IBC_INTERIOR, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    m1fP2C(5, :, IBC_INTERIOR, IACCU_CP4) = m1fP2C(5, :, IBC_INTERIOR, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD
    m1rP2C(5, :, IBC_INTERIOR, IACCU_CP4) = m1rP2C(5, :, IBC_INTERIOR, IACCU_CD2) ! 3 cell stencil, 4th CP --> 2nd CD

    m1fP2C(1, :, IBC_INTERIOR, IACCU_CP6) = m1fP2C(1, :, IBC_INTERIOR, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    m1rP2C(1, :, IBC_INTERIOR, IACCU_CP6) = m1rP2C(1, :, IBC_INTERIOR, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    m1fP2C(5, :, IBC_INTERIOR, IACCU_CP6) = m1fP2C(5, :, IBC_INTERIOR, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
    m1rP2C(5, :, IBC_INTERIOR, IACCU_CP6) = m1rP2C(5, :, IBC_INTERIOR, IACCU_CD4) ! 5 cell stencil, 6th CP --> 4th CD
!----------------------------------------------------------------------------------------------------------
!interpolation. P2C. IBC_DIRICHLET, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    m1fP2C(1:NL, 1:NS, IBC_DIRICHLET, 1:NACC) = m1fP2C(1:NL, 1:NS, IBC_INTERIOR, 1:NACC)
    m1rP2C(1:NL, 1:NS, IBC_DIRICHLET, 1:NACC) = m1rP2C(1:NL, 1:NS, IBC_INTERIOR, 1:NACC)
!----------------------------------------------------------------------------------------------------------
!interpolation. P2C. IBC_NEUMANN, unknowns only from only rhs could be reconstructed from bc, thus explicit
!----------------------------------------------------------------------------------------------------------
    m1fP2C(:, :, IBC_NEUMANN, :) = m1fP2C(:, :, IBC_INTERIOR, :)
    m1rP2C(:, :, IBC_NEUMANN, :) = m1rP2C(:, :, IBC_INTERIOR, :)
!----------------------------------------------------------------------------------------------------------
! interpolation. P2C: exterpolation
! [ 1    alpha1                          ][f_1']=[a1 * f_{1'} + b1 * f_{2'} + c1 * f_{3'}  ]
! [      alpha2 1     alpha2             ][f_2'] [a2/2 * (f_{2'}   + f_{3'})]
! [             alpha 1      alpha       ][f_i'] [a/2  * (f_{i'}   + f_{i'+1}) + b/2 * (f_{i'+2} + f_{i'-1})]
! [                   alpha2 1     alpha2][f_4'] [a2/2 * (f_{n'-1}   + f_{n'})]
! [                          alpha1 1    ][f_5'] [a1   * f_{n'+1} + b1 * f_{n'} + c1 * f_{n'-1}]
!-----------------------------------------------------------------------------
    alpha1(IACCU_CD2) = ZERO
        a1(IACCU_CD2) = HALF
        b1(IACCU_CD2) = HALF
        c1(IACCU_CD2) = ZERO
    alpha1(IACCU_CD4) = ZERO
        a1(IACCU_CD4) = THREE * EIGHTH
        b1(IACCU_CD4) = THREE * QUARTER
        c1(IACCU_CD4) = -ONE * EIGHTH
    alpha1(IACCU_CP4) = ONE_THIRD
        a1(IACCU_CP4) = ONE_THIRD
        b1(IACCU_CP4) = ONE
        c1(IACCU_CP4) = ZERO
    alpha1(IACCU_CP6) = ONE
        a1(IACCU_CP6) = QUARTER
        b1(IACCU_CP6) = ONEPFIVE
        c1(IACCU_CP6) = QUARTER
    do n = 1, NACC
      m1fP2C(1, 1,   IBC_INTRPL, n) = ZERO ! not used
      m1fP2C(1, 2,   IBC_INTRPL, n) = ONE
      m1fP2C(1, 3,   IBC_INTRPL, n) = alpha1(n)
      m1rP2C(1, 1,   IBC_INTRPL, n) = a1(n)
      m1rP2C(1, 2,   IBC_INTRPL, n) = b1(n)
      m1rP2C(1, 3,   IBC_INTRPL, n) = c1(n)

      m1fP2C(5, 1,   IBC_INTRPL, n) = m1fP2C(1, 3, IBC_INTRPL, n)
      m1fP2C(5, 2,   IBC_INTRPL, n) = m1fP2C(1, 2, IBC_INTRPL, n)
      m1fP2C(5, 3,   IBC_INTRPL, n) = m1fP2C(1, 1, IBC_INTRPL, n)
      m1rP2C(5, 1,   IBC_INTRPL, n) = m1rP2C(1, 1, IBC_INTRPL, n)
      m1rP2C(5, 2,   IBC_INTRPL, n) = m1rP2C(1, 2, IBC_INTRPL, n)
      m1rP2C(5, 3,   IBC_INTRPL, n) = m1rP2C(1, 3, IBC_INTRPL, n)

      m1fP2C(2:4, 1:3, IBC_INTRPL, n) = m1fP2C(2:4, 1:3, IBC_PERIODIC, n)
      m1rP2C(2:4, 1:3, IBC_INTRPL, n) = m1rP2C(2:4, 1:3, IBC_PERIODIC, n)
    end do
    
    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine Prepare_compact_coefficients
!==========================================================================================================
!> \brief Assigning the sparse matrix in the LHS of the compact scheme, and
!> calculating the geometry-only dependent variables for the TDMA scheme.
!>
!> This subroutine is called once locally.
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                           !
!----------------------------------------------------------------------------------------------------------
!> \param[in]     n             the number of unknown array
!> \param[in]     bc            the boundary condition at two ends of the unknown
!> \param[in]     coeff         the basic TDMA coefficients defined above.
!> \param[out]    a             the coefficients for TDMA
!> \param[out]    b             a_i * x_(i-1) + b_i * x_(i) + c_i * x_(i+1)
!> \param[out]    c             = RHS
!> \param[out]    d             An assisting coeffients for the TDMA scheme.
!----------------------------------------------------------------------------------------------------------
  subroutine Buildup_TDMA_LHS_array(n, coeff, a, b, c, d)
    use tridiagonal_matrix_algorithm
    use parameters_constant_mod
    implicit none

    integer, intent(in) :: n
    real(WP), intent(in)   :: coeff(NL, NS, NBCS:NBCE, NACC)
    real(WP), intent(out)  :: a(n, NBCS:NBCE, NBCS:NBCE, NACC), &
                              b(n, NBCS:NBCE, NBCS:NBCE, NACC), &
                              c(n, NBCS:NBCE, NBCS:NBCE, NACC), &
                              d(n, NBCS:NBCE, NBCS:NBCE, NACC)

    integer :: i, j, m, k, s

    a(:, :, :, :) =  ZERO
    b(:, :, :, :) =  ZERO
    c(:, :, :, :) =  ZERO
    d(:, :, :, :) =  ZERO
    k = IBC_PERIODIC
    do m = 1, NACC
      do j = NBCS, NBCE
        do i = NBCS, NBCE

          a(1,         i, j, m) = coeff( 1, 1, i, m)
          a(2,         i, j, m) = coeff( 2, 1, i, m)
          a(3 : n - 2, i, j, m) = coeff( 3, 1, k, m)
          a(n - 1,     i, j, m) = coeff( 4, 1, j, m)
          a(n,         i, j, m) = coeff( 5, 1, j, m)

          b(1,         i, j, m) = coeff( 1, 2, i, m)
          b(2,         i, j, m) = coeff( 2, 2, i, m)
          b(3 : n - 2, i, j, m) = coeff( 3, 2, k, m)
          b(n - 1,     i, j, m) = coeff( 4, 2, j, m)
          b(n,         i, j, m) = coeff( 5, 2, j, m)

          c(1,         i, j, m) = coeff( 1, 3, i, m)
          c(2,         i, j, m) = coeff( 2, 3, i, m)
          c(3 : n - 2, i, j, m) = coeff( 3, 3, k, m)
          c(n - 1,     i, j, m) = coeff( 4, 3, j, m)
          c(n,         i, j, m) = coeff( 5, 3, j, m)

          if (j == k .and. i == k) then
            call Preprocess_TDMA_coeffs( a(1:n-1, i, j, m), &
                                         b(1:n-1, i, j, m), &
                                         c(1:n-1, i, j, m), &
                                         d(1:n-1, i, j, m), &
                                         n-1)
          else 
            call Preprocess_TDMA_coeffs( a(:, i, j, m), &
                                         b(:, i, j, m), &
                                         c(:, i, j, m), &
                                         d(:, i, j, m), &
                                         n)
          end if 
        end do
      end do
    end do

    return
  end subroutine Buildup_TDMA_LHS_array
!==========================================================================================================
!> \brief Preparing coefficients for TDMA calculation.
!----------------------------------------------------------------------------------------------------------
!> Scope:  mpi    called-freq    xdomain
!>         all    once           all
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                           !
!----------------------------------------------------------------------------------------------------------
!==========================================================================================================
  subroutine Prepare_LHS_coeffs_for_operations
    use vars_df_mod, only : domain
    use mpi_mod
    use parameters_constant_mod
    implicit none
    integer :: i, nsz

!==========================================================================================================
!   building up the basic lhs coeffients for compact schemes
!==========================================================================================================
    call Prepare_compact_coefficients
!==========================================================================================================
!   building up the full size lhs coeffients for compact schemes
!==========================================================================================================
!----------------------------------------------------------------------------------------------------------
! y-direction, with nc unknows
!----------------------------------------------------------------------------------------------------------
    i = 2
    nsz = domain(1)%nc(i)

    allocate (ad1y_C2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); ad1y_C2C = ZERO
    allocate (bd1y_C2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); bd1y_C2C = ZERO
    allocate (cd1y_C2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); cd1y_C2C = ZERO
    allocate (dd1y_C2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); dd1y_C2C = ZERO
    call Buildup_TDMA_LHS_array(nsz, d1fC2C, &
          ad1y_C2C, bd1y_C2C, cd1y_C2C, dd1y_C2C)

    allocate (ad1y_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); ad1y_P2C = ZERO
    allocate (bd1y_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); bd1y_P2C = ZERO
    allocate (cd1y_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); cd1y_P2C = ZERO
    allocate (dd1y_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); dd1y_P2C = ZERO
    call Buildup_TDMA_LHS_array(nsz, d1fP2C, &
          ad1y_P2C, bd1y_P2C, cd1y_P2C, dd1y_P2C)

    allocate (am1y_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); am1y_P2C = ZERO
    allocate (bm1y_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); bm1y_P2C = ZERO
    allocate (cm1y_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); cm1y_P2C = ZERO
    allocate (dm1y_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); dm1y_P2C = ZERO
    call Buildup_TDMA_LHS_array(nsz, m1fP2C, &
          am1y_P2C, bm1y_P2C, cm1y_P2C, dm1y_P2C)

!----------------------------------------------------------------------------------------------------------
! y-direction, with np unknows
!----------------------------------------------------------------------------------------------------------
    nsz = domain(1)%np(i)

    allocate (ad1y_P2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); ad1y_P2P = ZERO
    allocate (bd1y_P2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); bd1y_P2P = ZERO
    allocate (cd1y_P2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); cd1y_P2P = ZERO
    allocate (dd1y_P2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); dd1y_P2P = ZERO
    call Buildup_TDMA_LHS_array(nsz, d1fP2P, &
          ad1y_P2P, bd1y_P2P, cd1y_P2P, dd1y_P2P)

    allocate (ad1y_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); ad1y_C2P = ZERO
    allocate (bd1y_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); bd1y_C2P = ZERO
    allocate (cd1y_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); cd1y_C2P = ZERO
    allocate (dd1y_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); dd1y_C2P = ZERO
    call Buildup_TDMA_LHS_array(nsz, d1fC2P, &
          ad1y_C2P, bd1y_C2P, cd1y_C2P, dd1y_C2P) 

    allocate (am1y_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); am1y_C2P = ZERO
    allocate (bm1y_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); bm1y_C2P = ZERO
    allocate (cm1y_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); cm1y_C2P = ZERO
    allocate (dm1y_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); dm1y_C2P = ZERO
    call Buildup_TDMA_LHS_array(nsz, m1fC2P, &
          am1y_C2P, bm1y_C2P, cm1y_C2P, dm1y_C2P)

!----------------------------------------------------------------------------------------------------------
! z-direction, with nc unknows
!----------------------------------------------------------------------------------------------------------
    i = 3
    nsz = domain(1)%nc(i)

    allocate (ad1z_C2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); ad1z_C2C = ZERO
    allocate (bd1z_C2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); bd1z_C2C = ZERO
    allocate (cd1z_C2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); cd1z_C2C = ZERO
    allocate (dd1z_C2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); dd1z_C2C = ZERO
    call Buildup_TDMA_LHS_array(nsz, d1fC2C, &
          ad1z_C2C, bd1z_C2C, cd1z_C2C, dd1z_C2C)

    allocate (ad1z_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); ad1z_P2C = ZERO
    allocate (bd1z_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); bd1z_P2C = ZERO
    allocate (cd1z_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); cd1z_P2C = ZERO
    allocate (dd1z_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); dd1z_P2C = ZERO
    call Buildup_TDMA_LHS_array(nsz, d1fP2C, &
          ad1z_P2C, bd1z_P2C, cd1z_P2C, dd1z_P2C)

    allocate (am1z_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); am1z_P2C = ZERO
    allocate (bm1z_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); bm1z_P2C = ZERO
    allocate (cm1z_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); cm1z_P2C = ZERO
    allocate (dm1z_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); dm1z_P2C = ZERO
    call Buildup_TDMA_LHS_array(nsz, m1fP2C, &
          am1z_P2C, bm1z_P2C, cm1z_P2C, dm1z_P2C)

!----------------------------------------------------------------------------------------------------------
! z-direction, with np unknows
!----------------------------------------------------------------------------------------------------------
    nsz = domain(1)%np(i)
!----------------------------------------------------------------------------------------------------------
! 1st derivative in z direction with np unknows
!----------------------------------------------------------------------------------------------------------
    allocate (ad1z_P2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); ad1z_P2P = ZERO
    allocate (bd1z_P2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); bd1z_P2P = ZERO
    allocate (cd1z_P2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); cd1z_P2P = ZERO
    allocate (dd1z_P2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); dd1z_P2P = ZERO
    call Buildup_TDMA_LHS_array(nsz, d1fP2P, &
          ad1z_P2P, bd1z_P2P, cd1z_P2P, dd1z_P2P)

    allocate (ad1z_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); ad1z_C2P = ZERO
    allocate (bd1z_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); bd1z_C2P = ZERO
    allocate (cd1z_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); cd1z_C2P = ZERO
    allocate (dd1z_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); dd1z_C2P = ZERO
    call Buildup_TDMA_LHS_array(nsz, d1fC2P, &
          ad1z_C2P, bd1z_C2P, cd1z_C2P, dd1z_C2P)

    allocate (am1z_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); am1z_C2P = ZERO
    allocate (bm1z_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); bm1z_C2P = ZERO
    allocate (cm1z_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); cm1z_C2P = ZERO
    allocate (dm1z_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); dm1z_C2P = ZERO
    call Buildup_TDMA_LHS_array(nsz, m1fC2P, &
          am1z_C2P, bm1z_C2P, cm1z_C2P, dm1z_C2P)

!----------------------------------------------------------------------------------------------------------
! x-direction
!----------------------------------------------------------------------------------------------------------
    allocate ( xtdma_lhs (nxdomain) )
    do i = 1, nxdomain
!----------------------------------------------------------------------------------------------------------
! x-direction, with nc unknows
!----------------------------------------------------------------------------------------------------------
      nsz = domain(i)%nc(1)

      allocate (xtdma_lhs(i)%ad1x_C2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%ad1x_C2C = ZERO
      allocate (xtdma_lhs(i)%bd1x_C2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%bd1x_C2C = ZERO
      allocate (xtdma_lhs(i)%cd1x_C2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%cd1x_C2C = ZERO
      allocate (xtdma_lhs(i)%dd1x_C2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%dd1x_C2C = ZERO
      call Buildup_TDMA_LHS_array(nsz , d1fC2C, &
            xtdma_lhs(i)%ad1x_C2C, &
            xtdma_lhs(i)%bd1x_C2C, &
            xtdma_lhs(i)%cd1x_C2C, &
            xtdma_lhs(i)%dd1x_C2C)

      allocate (xtdma_lhs(i)%ad1x_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%ad1x_P2C = ZERO
      allocate (xtdma_lhs(i)%bd1x_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%bd1x_P2C = ZERO
      allocate (xtdma_lhs(i)%cd1x_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%cd1x_P2C = ZERO
      allocate (xtdma_lhs(i)%dd1x_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%dd1x_P2C = ZERO
      call Buildup_TDMA_LHS_array(nsz , d1fP2C, &
            xtdma_lhs(i)%ad1x_P2C, &
            xtdma_lhs(i)%bd1x_P2C, &
            xtdma_lhs(i)%cd1x_P2C, &
            xtdma_lhs(i)%dd1x_P2C)

      allocate (xtdma_lhs(i)%am1x_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%am1x_P2C = ZERO
      allocate (xtdma_lhs(i)%bm1x_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%bm1x_P2C = ZERO
      allocate (xtdma_lhs(i)%cm1x_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%cm1x_P2C = ZERO
      allocate (xtdma_lhs(i)%dm1x_P2C ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%dm1x_P2C = ZERO
      call Buildup_TDMA_LHS_array(nsz , m1fP2C, &
          xtdma_lhs(i)%am1x_P2C, &
          xtdma_lhs(i)%bm1x_P2C, &
          xtdma_lhs(i)%cm1x_P2C, &
          xtdma_lhs(i)%dm1x_P2C)      

!----------------------------------------------------------------------------------------------------------
! x-direction, with np unknows
!----------------------------------------------------------------------------------------------------------
      nsz = domain(i)%np(1)

      allocate (xtdma_lhs(i)%ad1x_P2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%ad1x_P2P = ZERO
      allocate (xtdma_lhs(i)%bd1x_P2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%bd1x_P2P = ZERO
      allocate (xtdma_lhs(i)%cd1x_P2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%cd1x_P2P = ZERO
      allocate (xtdma_lhs(i)%dd1x_P2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%dd1x_P2P = ZERO
      call Buildup_TDMA_LHS_array(nsz , d1fP2P, &
            xtdma_lhs(i)%ad1x_P2P, &
            xtdma_lhs(i)%bd1x_P2P, &
            xtdma_lhs(i)%cd1x_P2P, &
            xtdma_lhs(i)%dd1x_P2P)

      allocate (xtdma_lhs(i)%ad1x_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%ad1x_C2P = ZERO
      allocate (xtdma_lhs(i)%bd1x_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%bd1x_C2P = ZERO
      allocate (xtdma_lhs(i)%cd1x_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%cd1x_C2P = ZERO
      allocate (xtdma_lhs(i)%dd1x_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%dd1x_C2P = ZERO
      call Buildup_TDMA_LHS_array(nsz , d1fC2P, &
            xtdma_lhs(i)%ad1x_C2P, &
            xtdma_lhs(i)%bd1x_C2P, &
            xtdma_lhs(i)%cd1x_C2P, &
            xtdma_lhs(i)%dd1x_C2P)

      allocate (xtdma_lhs(i)%am1x_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%am1x_C2P = ZERO
      allocate (xtdma_lhs(i)%bm1x_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%bm1x_C2P = ZERO
      allocate (xtdma_lhs(i)%cm1x_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%cm1x_C2P = ZERO
      allocate (xtdma_lhs(i)%dm1x_C2P ( nsz, NBCS:NBCE, NBCS:NBCE, NACC ) ); xtdma_lhs(i)%dm1x_C2P = ZERO
      call Buildup_TDMA_LHS_array(nsz , m1fC2P, &
          xtdma_lhs(i)%am1x_C2P, &
          xtdma_lhs(i)%bm1x_C2P, &
          xtdma_lhs(i)%cm1x_C2P, &
          xtdma_lhs(i)%dm1x_C2P)      

    end do

    return
  end subroutine Prepare_LHS_coeffs_for_operations
!==========================================================================================================
!==========================================================================================================
  subroutine buildup_ghost_cells_C(fi, d1, ibc, fc, fbc)
    use parameters_constant_mod
    implicit none
    real(WP), intent(in ) :: fi(:)
    integer,  intent(in ) :: ibc(2)
    real(WP), intent(in ) :: d1(4)
    real(WP), intent(inout ) :: fc(-1:2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nc ! cell number

    nc = size(fi)
!----------------------------------------------------------------------------------------------------------
!>                   BC                                BC
!>      _|__.__|__.__||__.__|__.__|__...___|__.__|__.__||__.__|__.__|__.__
!>         -1     0      1     2            nc-1   nc    nc+1   nc+2 
!----------------------------------------------------------------------------------------------------------
    if ( ibc(1) == IBC_INTERIOR) then
      if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_INTERIOR @ buildup_ghost_cells_C')
      fc(0 ) = fbc(1)
      fc(-1) = fbc(3)
    else if ( ibc(1) == IBC_PERIODIC) then
      fc(0 ) = fi(nc   )
      fc(-1) = fi(nc -1)
    else if (ibc(1) == IBC_SYMMETRIC ) then
      fc(0 ) = fi(1)
      fc(-1) = fi(2)
    else if (ibc(1) == IBC_ASYMMETRIC) then
      fc(0 ) = -fi(1)
      fc(-1) = -fi(2)
    else if (ibc(1) == IBC_DIRICHLET) then
      if(.not. present(fbc)) call Print_warning_msg('Lack of fbc info for IBC_DIRICHLET @ buildup_ghost_cells_C')
      fc(0 ) = TWO * fbc(1) - fi(1)
      fc(-1) = TWO * fbc(1) - fi(2)
    else if (ibc(1) == IBC_NEUMANN) then
      if(.not. present(fbc)) call Print_warning_msg('Lack of fbc info for IBC_NEUMANN @ buildup_ghost_cells_C')
      fc(0 ) = fi(1) - fbc(1) * d1(1)
      fc(-1) = fi(2) - fbc(1) * ( TWO * d1(3) + d1(1) )
    else
      fc(0 ) = MAXP
      fc(-1) = MAXP
    end if
!----------------------------------------------------------------------------------------------------------
!>                   BC                                BC
!>      _|__.__|__.__||__.__|__.__|__...___|__.__|__.__||__.__|__.__|__.__
!>         -1     0      1     2            nc-1   nc    nc+1   nc+2 
!----------------------------------------------------------------------------------------------------------
    if ( ibc(2) == IBC_INTERIOR) then
      if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_INTERIOR @ buildup_ghost_cells_C2C')
      fc(1) = fbc(2)
      fc(2) = fbc(4)
    else if ( ibc(2) == IBC_PERIODIC) then
      fc(1) = fi(1)
      fc(2) = fi(2)
    else if (ibc(2) == IBC_SYMMETRIC ) then
      fc(1) = fi(nc   )
      fc(2) = fi(nc -1)
    else if (ibc(2) == IBC_ASYMMETRIC) then
      fc(1) = -fi(nc   )
      fc(2) = -fi(nc -1)
    else if (ibc(2) == IBC_DIRICHLET) then
      if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_DIRICHLET @ buildup_ghost_cells_C2C')
      fc(1) = TWO * fbc(2) - fi(nc)
      fc(2) = TWO * fbc(2) - fi(nc-1)
    else if (ibc(2) == IBC_NEUMANN) then
      if(.not. present(fbc)) call Print_warning_msg('Lack of fbc info for IBC_NEUMANN @ buildup_ghost_cells_C2C')
      fc(1) = fi(nc)     + fbc(2) * d1(2)
      fc(2) = fi(nc - 1) + fbc(2) * ( TWO * d1(4) + d1(2) )
    else
      fc(1) = MAXP
      fc(2) = MAXP
    end if

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine buildup_ghost_cells_P(fi, d1, ibc, fp, fbc)
    use parameters_constant_mod
    implicit none
    real(WP), intent(in ) :: fi(:)
    integer,  intent(in ) :: ibc(2)
    real(WP), intent(in ) :: d1(4)
    real(WP), intent(inout ) :: fp(-1:2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: np

    np = size(fi)
!----------------------------------------------------------------------------------------------------------
!>                   BC                                BC
!>      _|__.__|__.__||__.__|__.__|__...___|__.__|__.__||__.__|__.__|__.__
!>      -1     0      1     2     3      np-2   np-1  np    np+1  np+2 (non-periodic)
!>      -1     0      1     2     3      np-1   np    np+1  np+2  np+3 (periodic)
!----------------------------------------------------------------------------------------------------------
    if ( ibc(1) == IBC_INTERIOR) then
      if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_INTERIOR @ buildup_ghost_cells_P2P')
      fp( 0) = fbc(1)
      fp(-1) = fbc(3)
    else if ( ibc(1) == IBC_PERIODIC) then
      fp( 0) = fi(np   )
      fp(-1) = fi(np -1)
    else if (ibc(1) == IBC_SYMMETRIC ) then
      fp( 0) = fi(2)
      fp(-1) = fi(3)
    else if (ibc(1) == IBC_ASYMMETRIC) then
      fp( 0) = -fi(2)
      fp(-1) = -fi(3)
    else if (ibc(1) == IBC_DIRICHLET) then
      if(present(fbc)) then
        fp( 0) = TWO * fbc(1) - fi(2)
        fp(-1) = TWO * fbc(1) - fi(3)
      else
        fp( 0) = TWO * fi(1) - fi(2)
        fp(-1) = TWO * fi(1) - fi(3)
      end if
    else if (ibc(1) == IBC_NEUMANN) then
      if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_NEUMANN @ buildup_ghost_cells_P2P')
      fp( 0) = fi(2) - fbc(1) * TWO * d1(1)
      fp(-1) = fi(3) - fbc(1) * TWO * ( d1(1) + d1(3) ) 
    else
      fp( 0) = MAXP
      fp(-1) = MAXP
    end if
!----------------------------------------------------------------------------------------------------------
!>                   BC                                BC
!>      _|__.__|__.__||__.__|__.__|__...___|__.__|__.__||__.__|__.__|__.__
!>      -1     0      1     2     3      np-2   np-1  np    np+1  np+2 (non-periodic)
!>      -1     0      1     2     3      np-1   np    np+1  np+2  np+3 (periodic)
!----------------------------------------------------------------------------------------------------------
    if ( ibc(2) == IBC_INTERIOR) then
      if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_INTERIOR @ buildup_ghost_cells_P2P')
      fp(1) = fbc(2)
      fp(2) = fbc(4)
    else if ( ibc(2) == IBC_PERIODIC) then
      fp(1) = fi(1)
      fp(2) = fi(2)
    else if (ibc(2) == IBC_SYMMETRIC ) then
      fp(1) = fi(np - 1)
      fp(2) = fi(np - 2)
    else if (ibc(2) == IBC_ASYMMETRIC) then
      fp(1) = - fi(np - 1)
      fp(2) = - fi(np - 2)
    else if (ibc(2) == IBC_NEUMANN) then
      if(.not. present(fbc)) call Print_error_msg('Lack of fbc info for IBC_NEUMANN @ buildup_ghost_cells_P2P')
      fp(1) = fi(np - 1) + fbc(2) * TWO * d1(2)
      fp(2) = fi(np - 2) + fbc(2) * TWO * ( d1(2) + d1(4) ) 
    else if (ibc(2) == IBC_DIRICHLET) then
      if(present(fbc)) then
        fp(1) = TWO * fbc(2) - fi(np - 1)
        fp(2) = TWO * fbc(2) - fi(np - 2)
      else
        fp(1) = TWO * fi(np) - fi(np - 1)
        fp(2) = TWO * fi(np) - fi(np - 2)
      end if
    else
      fp(1) = MAXP
      fp(2) = MAXP
    end if

    return
  end subroutine
!==========================================================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for interpolation.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> \param[in]   np      the number of unknowns, here is np
!> \param[in]   ibc     the b.c. type at two ends of the unknown array
!> \param[in]   fbc     the b.c. values for the given ibc
!> \param[in]   coeff   the defined TDMA coefficients
!> \param[in]   d1      spacing
!> \param[in]   fi      the input variable to build up the RHS array
!> \param[out]  fo      the output RHS array
!==========================================================================================================
  subroutine Prepare_TDMA_interp_P2C_RHS_array(fi, fo, nc, coeff, d1, ibc, fbc)
    use parameters_constant_mod
    implicit none
    real(WP), intent(in ) :: fi(:)
    integer,  intent(in ) :: nc ! unknow numbers, nc
    real(WP), intent(out) :: fo(nc)
    real(WP), intent(in ) :: coeff(1:NL, 1:NS, NBCS:NBCE)
    integer,  intent(in ) :: ibc(2)
    real(WP), intent(in ) :: d1(4)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: i!, m, l
    real(WP) :: fp(-1:2)
    logical :: is_bc_main(2), is_bc_extd(2)

    fo(:) = ZERO
!----------------------------------------------------------------------------------------------------------
!   i = bulk
!----------------------------------------------------------------------------------------------------------
    do i = 2, nc - 2
      fo(i) = coeff( 3, 1, IBC_PERIODIC ) * ( fi(i    ) + fi(i + 1) ) + &
              coeff( 3, 2, IBC_PERIODIC ) * ( fi(i - 1) + fi(i + 2) )
    end do
!----------------------------------------------------------------------------------------------------------
!>                   BC                                BC
!                        1     2             nc-1   nc    nc+1  nc+2
!>      _|__.__|__.__||__.__|__.__|__...___|__.__|__.__||__.__|__.__|__.__
!>      -1     0      1     2     3      np-2   np-1  np    np+1  np+2 (non-periodic)
!>      -1     0      1     2     3      np-1   np    np+1  np+2  np+3 (periodic)
!----------------------------------------------------------------------------------------------------------
    call buildup_ghost_cells_P(fi(:), d1(:), ibc(:), fp(-1:2), fbc(:))
    do i = 1, 2
      is_bc_main(i) = (ibc(i) == IBC_INTERIOR   .or. &
                   ibc(i) == IBC_SYMMETRIC  .or. &
                   ibc(i) == IBC_ASYMMETRIC .or. &
                   ibc(i) == IBC_DIRICHLET  .or. &
                   ibc(i) == IBC_NEUMANN )
      is_bc_extd(i) = (is_bc_main(i) .or. &
                   ibc(i) == IBC_PERIODIC )
    end do 
!----------------------------------------------------------------------------------------------------------
    i = 1
    if(is_bc_extd(1)) then
      fo(i) = coeff( 1, 1, ibc(1) ) * ( fi(i) + fi(i + 1) ) + &
              coeff( 1, 2, ibc(1) ) * ( fp(0) + fi(i + 2) )
    else
      fo(i) = coeff( 1, 1, IBC_INTRPL) * fi(i    ) + &
              coeff( 1, 2, IBC_INTRPL) * fi(i + 1) + &
              coeff( 1, 3, IBC_INTRPL) * fi(i + 2) 
    end if
!----------------------------------------------------------------------------------------------------------
    i = nc
    if(is_bc_main(2)) then
      fo(i) = coeff( 5, 1, ibc(2) ) * ( fi(i    ) + fi(i + 1) ) + &
              coeff( 5, 2, ibc(2) ) * ( fi(i - 1) + fp(1) )
    else if( ibc(2) == IBC_PERIODIC ) then
      fo(i) = coeff( 5, 1, IBC_PERIODIC ) * ( fi(i    ) + fp(1) ) + &
              coeff( 5, 2, IBC_PERIODIC ) * ( fi(i - 1) + fp(2) )
    else
      fo(nc) = coeff( 5, 1, IBC_INTRPL) * fi(i + 1) + &
               coeff( 5, 2, IBC_INTRPL) * fi(i    ) + &
               coeff( 5, 3, IBC_INTRPL) * fi(i - 1) 
    end if
!----------------------------------------------------------------------------------------------------------
    i = nc - 1
    if( ibc(2) == IBC_PERIODIC ) then
      fo(i) = coeff( 4, 1, IBC_PERIODIC ) * ( fi(i    ) + fi(i + 1) ) + &
              coeff( 4, 2, IBC_PERIODIC ) * ( fi(i - 1) + fp(1)     )
    else
      fo(i) = coeff( 4, 1, ibc(2) ) * ( fi(i    ) + fi(i + 1) ) + &
              coeff( 4, 2, ibc(2) ) * ( fi(i - 1) + fi(i + 2) )
    end if
!----------------------------------------------------------------------------------------------------------
!   mesh-based scaling
!----------------------------------------------------------------------------------------------------------
    ! nothing.
    return
  end subroutine Prepare_TDMA_interp_P2C_RHS_array

!==========================================================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for interpolation.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> \param[in]   np      the number of unknowns, here is np
!> \param[in]   ibc     the b.c. type at two ends of the unknown array
!> \param[in]   fbc     the b.c. values for the given ibc
!> \param[in]   coeff   the defined TDMA coefficients
!> \param[in]   d1      spacing
!> \param[in]   fi      the input variable to build up the RHS array
!> \param[out]  fo      the output RHS array
!==========================================================================================================
  subroutine Prepare_TDMA_interp_C2P_RHS_array(fi, fo, np, coeff, d1, ibc, fbc)
    use parameters_constant_mod
    implicit none

    real(WP),           intent(in ) :: fi(:)
    integer,            intent(in ) :: np ! unknow numbers, np
    real(WP),           intent(out) :: fo(np)
    real(WP),           intent(in ) :: coeff(1:NL, 1:NS, NBCS:NBCE)
    real(WP),           intent(in ) :: d1(4)
    integer,            intent(in ) :: ibc(2)
    real(WP), optional, intent(in)  :: fbc(4) ! used for Dirichlet B.C. (1||2) & interior (3, 1,|| 2, 4)

    integer :: i
    real(WP) :: fc(-1:2)
    logical :: is_bc_main(2), is_bc_extd(2)
!==========================================================================================================
!interpolation. C2P 
! alpha * f_{i'-1} + f_i' + alpha * f_{i'+1} = a/2 * ( f_{i}   + f_{i-1} ) + &
!                                              b/2 * ( f_{i+1} + f_{i-2} )
! when i' = 1,    need: f_{0'}, f_{0}, f_{-1}
! when i' = 2,    need:         f_{0}
! when i' = np-1, need:         f_{np'}
! when i' = np,   need: f_{np'+1}, f_{np'}, f_{np'+1}
!==========================================================================================================
    fo(:) = ZERO
!----------------------------------------------------------------------------------------------------------
!   i = bulk
!----------------------------------------------------------------------------------------------------------
    do i = 3, np - 2
      fo(i) = coeff( 3, 1, IBC_PERIODIC ) * ( fi(i    ) + fi(i - 1) ) + &
              coeff( 3, 2, IBC_PERIODIC ) * ( fi(i + 1) + fi(i - 2) )
    end do
!----------------------------------------------------------------------------------------------------------
!>                   BC                                BC
!>      _|__.__|__.__||__.__|__.__|__...___|__.__|__.__||__.__|__.__|__.__
!>         -1     0      1     2            nc-1   nc    nc+1   nc+2 
!----------------------------------------------------------------------------------------------------------
    call buildup_ghost_cells_C(fi(:), d1(:), ibc(:), fc(-1:2), fbc(:))
    is_bc_main(1) = (ibc(1) == IBC_INTERIOR   .or. &
                     ibc(1) == IBC_PERIODIC   .or. &
                     ibc(1) == IBC_SYMMETRIC  .or. &
                     ibc(1) == IBC_ASYMMETRIC .or. &
                     ibc(1) == IBC_NEUMANN )
    is_bc_extd(1) = (is_bc_main(1) .or. &
                     ibc(1) == IBC_DIRICHLET)
    is_bc_main(2) = (ibc(2) == IBC_INTERIOR   .or. &
                     ibc(2) == IBC_SYMMETRIC  .or. &
                     ibc(2) == IBC_ASYMMETRIC .or. &
                     ibc(2) == IBC_NEUMANN )
    is_bc_extd(2) = (is_bc_main(2) .or. &
                     ibc(2) == IBC_DIRICHLET)
!----------------------------------------------------------------------------------------------------------
    i = 1    
    if(is_bc_main(1)) then
      fo(i) = coeff( 1, 1, ibc(1)) * ( fi(i    ) + fc( 0) )+ &
              coeff( 1, 2, ibc(1)) * ( fi(i + 1) + fc(-1) )
    else if (ibc(1) == IBC_DIRICHLET) then
      fo(i) = fbc(1)
    else
      fo(i) = coeff( 1, 1, IBC_INTRPL) * fi(i    ) + &
              coeff( 1, 2, IBC_INTRPL) * fi(i + 1) + &
              coeff( 1, 3, IBC_INTRPL) * fi(i + 2) 
    end if
!----------------------------------------------------------------------------------------------------------
    i = 2 
    if(is_bc_main(1)) then
      fo(i) = coeff( 2, 1, ibc(1)) * ( fi(i    ) + fi(i - 1) ) + &
              coeff( 2, 2, ibc(1)) * ( fi(i + 1) + fc(0)     )
    else if (ibc(1) == IBC_DIRICHLET) then
      fo(i) = coeff( 2, 1, IBC_DIRICHLET) * ( fi(i    ) + fi(i - 1) ) + &
              coeff( 2, 2, IBC_DIRICHLET) * ( fi(i + 1) + fc(0)     )
    else
      fo(i) = coeff( 2, 1, IBC_INTRPL ) * ( fi(i    ) + fi(i - 1) )
    end if
!----------------------------------------------------------------------------------------------------------
    i = np
    if(is_bc_main(2)) then
      fo(i) = coeff( 5, 1, ibc(2) ) * ( fc(1) + fi(i - 1) ) + &
              coeff( 5, 2, ibc(2) ) * ( fc(2) + fi(i - 2) )
    else if (ibc(2) == IBC_PERIODIC) then
      fo(i) = coeff( 5, 1, IBC_PERIODIC ) * ( fi(i) + fi(i - 1) ) + &
              coeff( 5, 2, IBC_PERIODIC ) * ( fc(1) + fi(i - 2) )
    else if (ibc(2) == IBC_DIRICHLET) then
      fo(i) = fbc(2)
    else
      fo(i)   = coeff( 5, 1, IBC_INTRPL) * fi(i - 1) + &
                coeff( 5, 2, IBC_INTRPL) * fi(i - 2) + &
                coeff( 5, 3, IBC_INTRPL) * fi(i - 3) 
    end if
!----------------------------------------------------------------------------------------------------------
    i = np - 1
    if(is_bc_main(2)) then
      fo(i) = coeff( 4, 1, ibc(2) ) * ( fi(i) + fi(i - 1) ) + &
              coeff( 4, 2, ibc(2) ) * ( fc(1) + fi(i - 2) )
    else if (ibc(2) == IBC_PERIODIC) then
      fo(i) = coeff( 4, 1, IBC_PERIODIC ) * ( fi(i    ) + fi(i - 1) ) + &
              coeff( 4, 2, IBC_PERIODIC ) * ( fi(i + 1) + fi(i - 2) )
    else if (ibc(2) == IBC_DIRICHLET) then
      fo(i) = coeff( 4, 1, IBC_DIRICHLET ) * ( fi(i) + fi(i - 1) ) + &
              coeff( 4, 2, IBC_DIRICHLET ) * ( fc(1) + fi(i - 2) )
    else
      fo(i) = coeff( 4, 1, IBC_INTRPL ) * ( fi(i) + fi(i - 1) )
    end if
!----------------------------------------------------------------------------------------------------------
!   mesh-based scaling
!----------------------------------------------------------------------------------------------------------
    ! nothing
    return
  end subroutine Prepare_TDMA_interp_C2P_RHS_array
!==========================================================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for 1st derivative.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> for the 1st derivative.
!----------------------------------------------------------------------------------------------------------
! 1st derivative on collocated grids, C2C/P2P coefficients : Periodic or Symmetric B.C.
! alpha * f'_{i-1} + f'_i + alpha * f'_{i+1} = a/(2h) * ( f_{i+1} - f_{i-1} ) + &
!                                              b/(4h) * ( f_{i+2} - f_{i-2} )
!----------------------------------------------------------------------------------------------------------
!> \param[in]   nc      the number of unknowns, here is nc
!> \param[in]   ibc     the b.c. type at two ends of the unknown array
!> \param[in]   fbc     the b.c. values for the given ibc
!> \param[in]   coeff   the defined TDMA coefficients
!> \param[in]   dd      1/spacing, ie. 1/dx, 1/dy, 1/dz
!> \param[in]   fi      the input variable to build up the RHS array
!> \param[out]  fo      the output RHS array
!==========================================================================================================
  subroutine Prepare_TDMA_1deri_C2C_RHS_array(fi, fo, nc, coeff, d1, dd, ibc, fbc)
    use parameters_constant_mod
    implicit none
    real(WP), intent(in ) :: fi(:)
    integer,  intent(in ) :: nc ! unknow numbers
    real(WP), intent(out) :: fo(nc)
    real(WP), intent(in ) :: coeff(1:NL, 1:NS, NBCS:NBCE)
    real(WP), intent(in ) :: d1(4)
    real(WP), intent(in ) :: dd
    integer,  intent(in ) :: ibc(2)
    real(WP), optional, intent(in)  :: fbc(4) ! used for Dirichlet B.C. or interior

    integer :: i, m, l
    real(WP) :: fc(-1:2)
    logical :: is_bc(2)

    fo(:) = ZERO
!----------------------------------------------------------------------------------------------------------
!   i = bulk
!----------------------------------------------------------------------------------------------------------
    l = 3
    do i = 3, nc - 2
      fo(i) = coeff( l, 1, IBC_PERIODIC ) * ( fi(i + 1) - fi(i - 1) ) + &
              coeff( l, 2, IBC_PERIODIC ) * ( fi(i + 2) - fi(i - 2) )
    end do
!----------------------------------------------------------------------------------------------------------
!>                   BC                                BC
!>      _|__.__|__.__||__.__|__.__|__...___|__.__|__.__||__.__|__.__|__.__
!>         -1     0      1     2            nc-1   nc    nc+1   nc+2 
!----------------------------------------------------------------------------------------------------------
    call buildup_ghost_cells_C(fi(:), d1(:), ibc(:), fc(-1:2), fbc(:))
    do i = 1, 2 
      is_bc(i) = (ibc(i) == IBC_INTERIOR   .or. &
                  ibc(i) == IBC_PERIODIC   .or. &
                  ibc(i) == IBC_SYMMETRIC  .or. &
                  ibc(i) == IBC_ASYMMETRIC .or. &
                  ibc(i) == IBC_NEUMANN    .or. &
                  ibc(i) == IBC_DIRICHLET )
    end do
!----------------------------------------------------------------------------------------------------------
    i = 1
    if(is_bc(1)) then
      fo(i) = coeff( 1, 1, ibc(1) ) * ( fi(i + 1) - fc( 0) ) + &
              coeff( 1, 2, ibc(1) ) * ( fi(i + 2) - fc(-1) )

    else
      fo(i) = coeff( 1, 1, IBC_INTRPL) * fi(i    ) + &
              coeff( 1, 2, IBC_INTRPL) * fi(i + 1) + &
              coeff( 1, 3, IBC_INTRPL) * fi(i + 2) 
    end if
!----------------------------------------------------------------------------------------------------------
    i = 2
    if(is_bc(1)) then
      fo(i) = coeff( 2, 1, ibc(1) ) * ( fi(i + 1) - fi(i - 1) ) + &
              coeff( 2, 2, ibc(1) ) * ( fi(i + 2) - fc(0)     )
    else
      fo(i) = coeff( 2, 1, IBC_INTRPL ) * ( fi(i + 1) - fi(i - 1) )
    end if
!----------------------------------------------------------------------------------------------------------
    i = nc  
    if(is_bc(2)) then
      fo(i) = coeff( 5, 1, ibc(2) ) * ( fc(1) - fi(i - 1) ) + &
              coeff( 5, 2, ibc(2) ) * ( fc(2) - fi(i - 2) )
    else
      fo(i) = coeff( 5, 1, IBC_INTRPL) * fi(i    ) + &
              coeff( 5, 2, IBC_INTRPL) * fi(i - 1) + &
              coeff( 5, 3, IBC_INTRPL) * fi(i - 2) 
    end if
!----------------------------------------------------------------------------------------------------------
    i = nc - 1  
    if(is_bc(2)) then
      fo(i) = coeff( 4, 1, ibc(2) ) * ( fi(i + 1) - fi(i - 1) ) + &
              coeff( 4, 2, ibc(2) ) * ( fc(1)     - fi(i - 2) )
    else
      fo(i) = coeff( 4, 1, IBC_INTRPL ) * ( fi(i + 1) - fi(i - 1) )
    end if
!----------------------------------------------------------------------------------------------------------
!   mesh-based scaling
!----------------------------------------------------------------------------------------------------------
    fo(:) = fo(:) * dd

    return
  end subroutine Prepare_TDMA_1deri_C2C_RHS_array
!==========================================================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for 1st derivative.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!----------------------------------------------------------------------------------------------------------
! 1st derivative on collocated grids, C2C/P2P coefficients : Periodic or Symmetric B.C.
! alpha * f'_{i-1} + f'_i + alpha * f'_{i+1} = a/(2h) * ( f_{i+1} - f_{i-1} ) + &
!                                              b/(4h) * ( f_{i+2} - f_{i-2} )
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!> \param[in]   np      the number of unknowns, here is np
!> \param[in]   ibc     the b.c. type at two ends of the unknown array
!> \param[in]   fbc     the b.c. values for the given ibc
!> \param[in]   coeff   the defined TDMA coefficients
!> \param[in]   dd      1/spacing, ie. 1/dx, 1/dy, 1/dz
!> \param[in]   fi      the input variable to build up the RHS array
!> \param[out]  fo      the output RHS array
!==========================================================================================================
  subroutine Prepare_TDMA_1deri_P2P_RHS_array(fi, fo, np, coeff, d1, dd, ibc, fbc)
    use parameters_constant_mod
    implicit none
    real(WP),           intent(in ) :: fi(:)
    integer,            intent(in ) :: np ! unknow numbers
    real(WP),           intent(out) :: fo(np)
    real(WP),           intent(in ) :: coeff(1:NL, 1:NS, NBCS:NBCE)
    real(WP),           intent(in ) :: d1(4)
    real(WP),           intent(in ) :: dd
    integer,            intent(in ) :: ibc(2)
    real(WP), optional, intent(in ) :: fbc(4) ! used for IBC_NEUMANN or interior

    integer :: i
    real(WP) :: fp(-1:2)
    logical  :: is_bc_main(2), is_bc_extd(2)

    fo(:) = ZERO
!----------------------------------------------------------------------------------------------------------
!   i = bulk
!----------------------------------------------------------------------------------------------------------
    do i = 3, np - 2
      fo(i) = coeff( 3, 1, IBC_PERIODIC ) * ( fi(i + 1) - fi(i - 1) ) + &
              coeff( 3, 2, IBC_PERIODIC ) * ( fi(i + 2) - fi(i - 2) )
    end do
!----------------------------------------------------------------------------------------------------------
!>                   BC                                BC
!>      _|__.__|__.__||__.__|__.__|__...___|__.__|__.__||__.__|__.__|__.__
!>      -1     0      1     2     3      np-2   np-1  np    np+1  np+2 (non-periodic)
!>      -1     0      1     2     3      np-1   np    np+1  np+2  np+3 (periodic)
!----------------------------------------------------------------------------------------------------------
    call buildup_ghost_cells_P(fi(:), d1(:), ibc(:), fp(-1:2), fbc(:))
    do i = 1, 2
      is_bc_main(i) = (ibc(i) == IBC_INTERIOR   .or. &
                   ibc(i) == IBC_PERIODIC   .or. &
                   ibc(i) == IBC_SYMMETRIC  .or. &
                   ibc(i) == IBC_ASYMMETRIC .or. &
                   ibc(i) == IBC_DIRICHLET )
      is_bc_extd(i) = (is_bc_main(i) .or. &
                   ibc(i) == IBC_NEUMANN)
    end do
!----------------------------------------------------------------------------------------------------------
    i = 1
    if(is_bc_main(1)) then
      fo(i) = coeff( 1, 1, ibc(1) ) * ( fi(i + 1) - fp( 0) ) + &
              coeff( 1, 2, ibc(1) ) * ( fi(i + 2) - fp(-1) )
    else if(ibc(1) == IBC_NEUMANN) then
      fo(i) = fbc(1)/dd
    else
      fo(i) = coeff( 1, 1, IBC_INTRPL) * fi(i    ) + &
              coeff( 1, 2, IBC_INTRPL) * fi(i + 1) + &
              coeff( 1, 3, IBC_INTRPL) * fi(i + 2) 
    end if

!----------------------------------------------------------------------------------------------------------    
    i = 2
    if(is_bc_extd(1)) then
      fo(i) = coeff( 2, 1, ibc(1) ) * ( fi(i + 1) - fi(i - 1) ) + &
              coeff( 2, 2, ibc(1) ) * ( fi(i + 2) - fp(0)     )
    else
      fo(i) = coeff( 2, 1, IBC_INTRPL ) * ( fi(i + 1) - fi(i - 1) )
    end if
!----------------------------------------------------------------------------------------------------------
    i = np
    if(is_bc_main(2)) then
      fo(i) = coeff( 5, 1, ibc(2) ) * ( fp(1) - fi(i - 1) ) + &
              coeff( 5, 2, ibc(2) ) * ( fp(2) - fi(i - 2) )
    else if(ibc(2) == IBC_NEUMANN) then
      fo(i) = fbc(2)/dd
    else
      fo(i) = coeff( 5, 1, IBC_INTRPL) * fi(i    ) + &
              coeff( 5, 2, IBC_INTRPL) * fi(i - 1) + &
              coeff( 5, 3, IBC_INTRPL) * fi(i - 2) 
    end if
!----------------------------------------------------------------------------------------------------------
    i = np - 1
    if(is_bc_extd(2)) then
      fo(i) = coeff( 4, 1, ibc(2) ) * ( fi(i + 1) - fi(i - 1) ) + &
              coeff( 4, 2, ibc(2) ) * ( fp(1)     - fi(i - 2) )
    else
      fo(i) = coeff( 4, 1, IBC_INTRPL ) * ( fi(i + 1) - fi(i - 1) )
    end if
!----------------------------------------------------------------------------------------------------------
!   mesh-based scaling
!----------------------------------------------------------------------------------------------------------
    fo(:) = fo(:) * dd

    return
  end subroutine Prepare_TDMA_1deri_P2P_RHS_array
!==========================================================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for 1st derivative.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!----------------------------------------------------------------------------------------------------------
! 1st derivative on staggered grids C2P
! C2P ==>
! alpha * f'_{i'-1} + f'_i' + alpha * f'_{i'+1} = a/(h ) * ( f_{i}   - f_{i-1} ) + &
!                                                 b/(3h) * ( f_{i+1} - f_{i-2} )
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!> \param[in]   np      the number of unknowns, here is np
!> \param[in]   ibc     the b.c. type at two ends of the unknown array
!> \param[in]   fbc     the b.c. values for the given ibc
!> \param[in]   coeff   the defined TDMA coefficients
!> \param[in]   dd      1/spacing, ie. 1/dx, 1/dy, 1/dz
!> \param[in]   fi      the input variable to build up the RHS array
!> \param[out]  fo      the output RHS array
!==========================================================================================================
  subroutine Prepare_TDMA_1deri_C2P_RHS_array(fi, fo, np, coeff, d1, dd, ibc, fbc)
    use parameters_constant_mod
    implicit none
    real(WP),           intent(in ) :: fi(:)
    integer,            intent(in ) :: np ! unknow numbers, np
    real(WP),           intent(out) :: fo(np)
    real(WP),           intent(in ) :: coeff(1:NL, 1:NS, NBCS:NBCE)
    real(WP),           intent(in ) :: d1(4)
    real(WP),           intent(in ) :: dd
    integer,            intent(in ) :: ibc(2)
    real(WP), optional, intent(in ) :: fbc(4) ! used for IBC_NEUMANN, and interior

    integer :: i!, m, l
    real(WP) :: fc(-1:2)
    logical  :: is_bc_main(2), is_bc_extd(2)

    fo(:) = ZERO
!----------------------------------------------------------------------------------------------------------
!   i = bulk
!----------------------------------------------------------------------------------------------------------
    do i = 3, np - 2
      fo(i) = coeff( 3, 1, IBC_PERIODIC ) * ( fi(i    ) - fi(i - 1) ) + &
              coeff( 3, 2, IBC_PERIODIC ) * ( fi(i + 1) - fi(i - 2) )
    end do
!----------------------------------------------------------------------------------------------------------
!>                   BC                                BC
!>      _|__.__|__.__||__.__|__.__|__...___|__.__|__.__||__.__|__.__|__.__
!>         -1     0      1     2            nc-1   nc    nc+1   nc+2 
!----------------------------------------------------------------------------------------------------------
    call buildup_ghost_cells_C(fi(:), d1(:), ibc(:), fc(-1:2), fbc(:))
    is_bc_main(1) = (ibc(1) == IBC_INTERIOR   .or. &
                     ibc(1) == IBC_PERIODIC   .or. &
                     ibc(1) == IBC_SYMMETRIC  .or. &
                     ibc(1) == IBC_ASYMMETRIC .or. &
                     ibc(1) == IBC_DIRICHLET )
    is_bc_extd(1) = (is_bc_main(1) .or. &
                     ibc(1) == IBC_NEUMANN)
    is_bc_main(2) = (ibc(2) == IBC_INTERIOR   .or. &
                     ibc(2) == IBC_SYMMETRIC  .or. &
                     ibc(2) == IBC_ASYMMETRIC .or. &
                     ibc(2) == IBC_DIRICHLET )
    is_bc_extd(2) = (is_bc_main(2) .or. &
                     ibc(2) == IBC_NEUMANN)
!----------------------------------------------------------------------------------------------------------
    i = 1
    if(is_bc_main(1)) then
      fo(i) = coeff( 1, 1, ibc(1) ) * ( fi(i    ) - fc( 0) ) + &
              coeff( 1, 2, ibc(1) ) * ( fi(i + 1) - fc(-1) )
    else if(ibc(1) == IBC_NEUMANN) then
      fo(i) = fbc(1)/dd
    else
      fo(i) = coeff( 1, 1, IBC_INTRPL) * fi(i    ) + &
              coeff( 1, 2, IBC_INTRPL) * fi(i + 1) + &
              coeff( 1, 3, IBC_INTRPL) * fi(i + 2) 
    end if
!----------------------------------------------------------------------------------------------------------
    i = 2
    if(is_bc_extd(1)) then
      fo(i) = coeff( 2, 1, ibc(1) ) * ( fi(i    ) - fi(i - 1) ) + &
              coeff( 2, 2, ibc(1) ) * ( fi(i + 1) - fc(0)     )
    else
      fo(i) = coeff( 2, 1, IBC_PERIODIC ) * ( fi(i    ) - fi(i - 1) )
    end if
!----------------------------------------------------------------------------------------------------------
    i = np
    if(is_bc_main(2)) then
      fo(i) = coeff( 5, 1, ibc(2) ) * ( fc(1) - fi(i - 1) ) + &
              coeff( 5, 2, ibc(2) ) * ( fc(2) - fi(i - 2) )
    else if(ibc(2) == IBC_PERIODIC) then
      fo(i) = coeff( 5, 1, IBC_PERIODIC ) * ( fi(i) - fi(i - 1) ) + &
              coeff( 5, 2, IBC_PERIODIC ) * ( fc(1) - fi(i - 2) )
    else if(ibc(2) == IBC_NEUMANN) then
      fo(i) = fbc(2)/dd
    else
      fo(i) = coeff( 5, 1, IBC_INTRPL) * fi(i - 1) + &
              coeff( 5, 2, IBC_INTRPL) * fi(i - 2) + &
              coeff( 5, 3, IBC_INTRPL) * fi(i - 3) 
    end if
!----------------------------------------------------------------------------------------------------------
    i = np - 1
    if(is_bc_extd(2)) then
      fo(i) = coeff( 4, 1, ibc(2) ) * ( fi(i) - fi(i - 1) ) + &
              coeff( 4, 2, ibc(2) ) * ( fc(1) - fi(i - 2) )
    else if(ibc(2) == IBC_PERIODIC) then
      fo(i) = coeff( 4, 1, IBC_PERIODIC ) * ( fi(i)   - fi(i - 1) ) + &
              coeff( 4, 2, IBC_PERIODIC ) * ( fi(i+1) - fi(i - 2) )
    else
      fo(i) = coeff( 4, 1, IBC_PERIODIC ) * ( fi(i) - fi(i - 1) )
    end if
!----------------------------------------------------------------------------------------------------------
!   mesh-based scaling
!----------------------------------------------------------------------------------------------------------
    fo(:) = fo(:) * dd

    return
  end subroutine Prepare_TDMA_1deri_C2P_RHS_array

!==========================================================================================================
!> \brief Preparing the RHS array for the TDMA algorithm for 1st derivative - P2C.
!> This subroutine is called repeatly to update the RHS of the TDMA algorithm
!> \param[in]   nc      the number of unknowns, here is nc
!> \param[in]   ibc     the b.c. type at two ends of the unknown array
!> \param[in]   fbc     the b.c. values for the given ibc
!> \param[in]   coeff   the defined TDMA coefficients
!> \param[in]   dd      1/spacing, ie. 1/dx, 1/dy, 1/dz
!> \param[in]   fi      the input variable to build up the RHS array
!> \param[out]  fo      the output RHS array
!==========================================================================================================
  subroutine Prepare_TDMA_1deri_P2C_RHS_array(fi, fo, nc, coeff, d1, dd, ibc, fbc)
    use parameters_constant_mod
    implicit none
    real(WP), intent(in ) :: fi(:)
    integer,  intent(in ) :: nc ! unknow numbers, nc
    real(WP), intent(out) :: fo(nc)
    real(WP), intent(in ) :: coeff(1:NL, 1:NS, NBCS:NBCE)
    real(WP), intent(in ) :: dd
    real(WP), intent(in ) :: d1(4)
    integer,  intent(in ) :: ibc(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: i!, l, m
    real(WP) :: fp(-1:2)
    logical :: is_bc_main(2), is_bc_extd(2)

    fo(:) = ZERO
!----------------------------------------------------------------------------------------------------------
!   i = bulk
!----------------------------------------------------------------------------------------------------------
    do i = 2, nc - 2
      fo(i) = coeff( 3, 1, IBC_PERIODIC ) * ( fi(i + 1) - fi(i    ) ) + &
              coeff( 3, 2, IBC_PERIODIC ) * ( fi(i + 2) - fi(i - 1) )
    end do
!----------------------------------------------------------------------------------------------------------
!>                   BC                                BC
!>      _|__.__|__.__||__.__|__.__|__...___|__.__|__.__||__.__|__.__|__.__
!>      -1     0      1     2     3      np-2   np-1  np    np+1  np+2 (non-periodic)
!>      -1     0      1     2     3      np-1   np    np+1  np+2  np+3 (periodic)
!----------------------------------------------------------------------------------------------------------
    call buildup_ghost_cells_P(fi(:), d1(:), ibc(:), fp(-1:2), fbc(:))
    do i = 1, 2
      is_bc_main(i) = (ibc(i) == IBC_INTERIOR   .or. &
                       ibc(i) == IBC_SYMMETRIC  .or. &
                       ibc(i) == IBC_ASYMMETRIC .or. &
                       ibc(i) == IBC_NEUMANN    .or. &
                       ibc(i) == IBC_DIRICHLET )
      is_bc_extd(i) = (is_bc_main(i) .or. &
                       ibc(i) == IBC_PERIODIC )
    end do
!----------------------------------------------------------------------------------------------------------
    i = 1
    if(is_bc_extd(1)) then
      fo(i) = coeff( 1, 1, ibc(1) ) * ( fi(i + 1) - fi(i) ) + &
              coeff( 1, 2, ibc(1) ) * ( fi(i + 2) - fp(0) )
    else
      fo(i) = coeff( 1, 1, IBC_INTRPL) * fi(i    ) + &
              coeff( 1, 2, IBC_INTRPL) * fi(i + 1) + &
              coeff( 1, 3, IBC_INTRPL) * fi(i + 2) 
    end if
!----------------------------------------------------------------------------------------------------------
    i = nc
    if(is_bc_main(2)) then
      fo(i) = coeff( 5, 1, ibc(2) ) * ( fi(i + 1) - fi(i    ) ) + &
              coeff( 5, 2, ibc(2) ) * ( fp(1)     - fi(i - 1) )
    else if(ibc(2) == IBC_PERIODIC) then
      fo(i) = coeff( 5, 1, IBC_PERIODIC ) * ( fp(1) - fi(i    ) ) + &
              coeff( 5, 2, IBC_PERIODIC ) * ( fp(2) - fi(i - 1) )
    else
      fo(i) = coeff( 5, 1, IBC_INTRPL) * fi(i + 1) + &
              coeff( 5, 2, IBC_INTRPL) * fi(i    ) + &
              coeff( 5, 3, IBC_INTRPL) * fi(i - 1) 
    end if
!----------------------------------------------------------------------------------------------------------
    i = nc - 1
    if(is_bc_main(2)) then
      fo(i) = coeff( 4, 1, ibc(2) ) * ( fi(i + 1) - fi(i    ) ) + &
              coeff( 4, 2, ibc(2) ) * ( fi(i + 2) - fi(i - 1) )
    else if(ibc(2) == IBC_PERIODIC) then
      fo(i) = coeff( 4, 1, IBC_PERIODIC ) * ( fi(i + 1) - fi(i    ) ) + &
              coeff( 4, 2, IBC_PERIODIC ) * ( fp(1)     - fi(i - 1) )
    else
      fo(i) = coeff( 4, 1, IBC_PERIODIC ) * ( fi(i + 1) - fi(i    ) ) + &
              coeff( 4, 2, IBC_PERIODIC ) * ( fi(i + 2) - fi(i - 1) )
    end if
!----------------------------------------------------------------------------------------------------------
!   mesh-based scaling
!----------------------------------------------------------------------------------------------------------
    fo(:) = fo(:) * dd

    return
  end subroutine Prepare_TDMA_1deri_P2C_RHS_array
!==========================================================================================================
!> \brief To caculate the mid-point interpolation in 1D.
!> This subroutine is called as required to get the mid-point interpolation.
!---------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   pubic
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixsub         x-subdomain index
!> \param[in]     ibc           bc type
!> \param[in]     fbc           bc value
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!_______________________________________________________________________________
  subroutine Get_x_midp_C2P_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: ixsub, nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic

    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @ Get_x_midp_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_DIRICHLET  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_DIRICHLET, degragded to IBC_INTRPL. @ Get_x_midp_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @ Get_x_midp_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    ixsub = dm%idom
    nsz = size(fo)
    fo = ZERO
    d1(:) = dm%h(1)
    call Prepare_TDMA_interp_C2P_RHS_array(fi(:), fo(:), nsz, m1rC2P(:, :, :, iacc), d1(:), ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
            xtdma_lhs(ixsub)%am1x_C2P(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%bm1x_C2P(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%cm1x_C2P(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%dm1x_C2P(:, ibc(1), ibc(2), iacc), &
            nsz)
      end if
    return
  end subroutine Get_x_midp_C2P_1D
!==========================================================================================================
  subroutine Get_x_midp_P2C_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: ixsub, nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_x_midp_P2C_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_x_midp_P2C_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    ixsub = dm%idom
    d1(:) = dm%h(1)

    call Prepare_TDMA_interp_P2C_RHS_array(fi(:), fo(:), nsz, m1rP2C(1:NL, 1:NS, NBCS:NBCE, iacc), d1, ibc(:), fbc)
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
            xtdma_lhs(ixsub)%am1x_P2C(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%bm1x_P2C(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%cm1x_P2C(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%dm1x_P2C(:, ibc(1), ibc(2), iacc), &
            nsz)
    end if

    return
  end subroutine Get_x_midp_P2C_1D
!==========================================================================================================
  subroutine Get_y_midp_C2P_1D(fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_y_midp_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_DIRICHLET  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_DIRICHLET, degragded to IBC_INTRPL. @Get_y_midp_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_y_midp_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    d1(1) = dm%yp(2) - dm%yp(1)
    d1(3) = dm%yp(3) - dm%yp(2)
    d1(2) = dm%yp(dm%np(2)  ) - dm%yp(dm%np(2)-1)
    d1(4) = dm%yp(dm%np(2)-1) - dm%yp(dm%np(2)-2)
    call Prepare_TDMA_interp_C2P_RHS_array(fi(:), fo(:), nsz, m1rC2P(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
            am1y_C2P(:, ibc(1), ibc(2), iacc), &
            bm1y_C2P(:, ibc(1), ibc(2), iacc), &
            cm1y_C2P(:, ibc(1), ibc(2), iacc), &
            dm1y_C2P(:, ibc(1), ibc(2), iacc), &
            nsz)
    end if
! stretching? No stretching conversion
    return
  end subroutine Get_y_midp_C2P_1D
!==========================================================================================================
  subroutine Get_y_midp_P2C_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_y_midp_P2C_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_y_midp_P2C_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    d1(1) = dm%yp(2) - dm%yp(1)
    d1(3) = dm%yp(3) - dm%yp(2)
    d1(2) = dm%yp(dm%np(2))   - dm%yp(dm%np(2)-1)
    d1(4) = dm%yp(dm%np(2)-1) - dm%yp(dm%np(2)-2)
    call Prepare_TDMA_interp_P2C_RHS_array(fi(:), fo(:), nsz, m1rP2C(1:NL, 1:NS, NBCS:NBCE, iacc), d1, ibc(:), fbc(:) )
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
            am1y_P2C(:, ibc(1), ibc(2), iacc), &
            bm1y_P2C(:, ibc(1), ibc(2), iacc), &
            cm1y_P2C(:, ibc(1), ibc(2), iacc), &
            dm1y_P2C(:, ibc(1), ibc(2), iacc), &
            nsz)
    end if

    return
  end subroutine Get_y_midp_P2C_1D
!==========================================================================================================
  subroutine Get_z_midp_C2P_1D(fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_z_midp_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_DIRICHLET  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_DIRICHLET, degragded to IBC_INTRPL. @Get_z_midp_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_z_midp_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    d1(:) =  dm%h(3)
    call Prepare_TDMA_interp_C2P_RHS_array(fi(:), fo(:), nsz, m1rC2P(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
            am1z_C2P(:, ibc(1), ibc(2), iacc), &
            bm1z_C2P(:, ibc(1), ibc(2), iacc), &
            cm1z_C2P(:, ibc(1), ibc(2), iacc), &
            dm1z_C2P(:, ibc(1), ibc(2), iacc), &
            nsz)
    end if

    return
  end subroutine Get_z_midp_C2P_1D
!==========================================================================================================
  subroutine Get_z_midp_P2C_1D(fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_z_midp_P2C_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_z_midp_P2C_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    d1(:) = dm%h(3)
    call Prepare_TDMA_interp_P2C_RHS_array(fi(:), fo(:), nsz, m1rP2C(1:NL, 1:NS, NBCS:NBCE, iacc), d1, ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
            am1z_P2C(:, ibc(1), ibc(2), iacc), &
            bm1z_P2C(:, ibc(1), ibc(2), iacc), &
            cm1z_P2C(:, ibc(1), ibc(2), iacc), &
            dm1z_P2C(:, ibc(1), ibc(2), iacc), &
            nsz)
    end if

    return
  end subroutine Get_z_midp_P2C_1D
!==========================================================================================================
!> \brief To caculate the 1st derivative in 1D.
!> This subroutine is called as required to get the 1st derivative
!---------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   pubic
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixsub         x-subdomain index
!> \param[in]     ibc           bc type
!> \param[in]     fbc           bc value
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!==========================================================================================================
  subroutine Get_x_1st_derivative_C2C_1D(fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: ixsub, nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_x_1st_derivative_C2C_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_DIRICHLET  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_DIRICHLET, degragded to IBC_INTRPL. @Get_x_1st_derivative_C2C_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_x_1st_derivative_C2C_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    ixsub = dm%idom
    nsz = size(fo)
    fo = ZERO
    d1(:) = dm%h(1)

    call Prepare_TDMA_1deri_C2C_RHS_array(fi(:), fo(:), nsz, d1rC2C(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), dm%h1r(1), ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
            xtdma_lhs(ixsub)%ad1x_C2C(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%bd1x_C2C(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%cd1x_C2C(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%dd1x_C2C(:, ibc(1), ibc(2), iacc), &
            nsz)
    end if

    return
  end subroutine Get_x_1st_derivative_C2C_1D
!==========================================================================================================
  subroutine Get_x_1st_derivative_P2P_1D(fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: ixsub, nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_x_1st_derivative_P2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_x_1st_derivative_P2P_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    ixsub = dm%idom
    d1(:) = dm%h(1)
    call Prepare_TDMA_1deri_P2P_RHS_array(fi(:), fo(:), nsz, d1rP2P(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), dm%h1r(1), ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
            xtdma_lhs(ixsub)%ad1x_P2P(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%bd1x_P2P(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%cd1x_P2P(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%dd1x_P2P(:, ibc(1), ibc(2), iacc), &
            nsz)
    end if
    return
  end subroutine Get_x_1st_derivative_P2P_1D
!==========================================================================================================
  subroutine Get_x_1st_derivative_C2P_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: ixsub, nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_x_1st_derivative_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_DIRICHLET  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_DIRICHLET, degragded to IBC_INTRPL. @Get_x_1st_derivative_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_x_1st_derivative_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO

    ixsub = dm%idom
    d1(:) = dm%h(1)
    call Prepare_TDMA_1deri_C2P_RHS_array(fi(:), fo(:), nsz, d1rC2P(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), dm%h1r(1), ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
            xtdma_lhs(ixsub)%ad1x_C2P(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%bd1x_C2P(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%cd1x_C2P(:, ibc(1), ibc(2), iacc), &
            xtdma_lhs(ixsub)%dd1x_C2P(:, ibc(1), ibc(2), iacc), &
            nsz)
    end if
    return
  end subroutine Get_x_1st_derivative_C2P_1D
!==========================================================================================================
  subroutine Get_x_1st_derivative_P2C_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: ixsub, nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_x_1st_derivative_P2C_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_x_1st_derivative_P2C_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO

    ixsub = dm%idom
    d1(:) = dm%h(1)
    call Prepare_TDMA_1deri_P2C_RHS_array(fi(:), fo(:), nsz, d1rP2C(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), dm%h1r(1), ibc(:), fbc(:) )
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
          xtdma_lhs(ixsub)%ad1x_P2C(:, ibc(1), ibc(2), iacc), &
          xtdma_lhs(ixsub)%bd1x_P2C(:, ibc(1), ibc(2), iacc), &
          xtdma_lhs(ixsub)%cd1x_P2C(:, ibc(1), ibc(2), iacc), &
          xtdma_lhs(ixsub)%dd1x_P2C(:, ibc(1), ibc(2), iacc), &
          nsz)
    end if

    return
  end subroutine Get_x_1st_derivative_P2C_1D
!==========================================================================================================
! y - Get_1st_derivative_1D
!==========================================================================================================
  subroutine Get_y_1st_derivative_C2C_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_y_1st_derivative_C2C_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_DIRICHLET  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_DIRICHLET, degragded to IBC_INTRPL. @Get_y_1st_derivative_C2C_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_y_1st_derivative_C2C_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    d1(1) = dm%yp(2) - dm%yp(1)
    d1(3) = dm%yp(3) - dm%yp(2)
    d1(2) = dm%yp(dm%np(2))   - dm%yp(dm%np(2)-1)
    d1(4) = dm%yp(dm%np(2)-1) - dm%yp(dm%np(2)-2)
    call Prepare_TDMA_1deri_C2C_RHS_array(fi(:), fo(:), nsz, d1rC2C(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), dm%h1r(2), ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
          ad1y_C2C(:, ibc(1), ibc(2), iacc), &
          bd1y_C2C(:, ibc(1), ibc(2), iacc), &
          cd1y_C2C(:, ibc(1), ibc(2), iacc), &
          dd1y_C2C(:, ibc(1), ibc(2), iacc), &
          nsz)
    end if

    if(dm%is_stretching(2)) fo(:) = fo(:) * dm%yMappingcc(:, 1)

    return
  end subroutine Get_y_1st_derivative_C2C_1D
!==========================================================================================================
  subroutine Get_y_1st_derivative_P2P_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_y_1st_derivative_P2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_y_1st_derivative_P2P_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    d1(1) = dm%yp(2) - dm%yp(1)
    d1(3) = dm%yp(3) - dm%yp(2)
    d1(2) = dm%yp(dm%np(2))   - dm%yp(dm%np(2)-1)
    d1(4) = dm%yp(dm%np(2)-1) - dm%yp(dm%np(2)-2)
    call Prepare_TDMA_1deri_P2P_RHS_array(fi(:), fo(:), nsz, d1rP2P(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), dm%h1r(2), ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
          ad1y_P2P(:, ibc(1), ibc(2), iacc), &
          bd1y_P2P(:, ibc(1), ibc(2), iacc), &
          cd1y_P2P(:, ibc(1), ibc(2), iacc), &
          dd1y_P2P(:, ibc(1), ibc(2), iacc), &
          nsz)
    end if
    
    if(dm%is_stretching(2)) fo(:) = fo(:) * dm%yMappingpt(:, 1)

    return
  end subroutine Get_y_1st_derivative_P2P_1D
!==========================================================================================================
  subroutine Get_y_1st_derivative_C2P_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_y_1st_derivative_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_DIRICHLET  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_DIRICHLET, degragded to IBC_INTRPL. @Get_y_1st_derivative_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_y_1st_derivative_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    d1(1) = dm%yp(2) - dm%yp(1)
    d1(3) = dm%yp(3) - dm%yp(2)
    d1(2) = dm%yp(dm%np(2))   - dm%yp(dm%np(2)-1)
    d1(4) = dm%yp(dm%np(2)-1) - dm%yp(dm%np(2)-2)
    call Prepare_TDMA_1deri_C2P_RHS_array(fi(:), fo(:), nsz, d1rC2P(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), dm%h1r(2), ibc(:), fbc(:) )
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
            ad1y_C2P(:, ibc(1), ibc(2), iacc), &
            bd1y_C2P(:, ibc(1), ibc(2), iacc), &
            cd1y_C2P(:, ibc(1), ibc(2), iacc), &
            dd1y_C2P(:, ibc(1), ibc(2), iacc), &
            nsz)
    end if
  
    if(dm%is_stretching(2)) fo(:) = fo(:) * dm%yMappingpt(:, 1)

    return
  end subroutine Get_y_1st_derivative_C2P_1D
!==========================================================================================================
  subroutine Get_y_1st_derivative_P2C_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_y_1st_derivative_P2C_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_y_1st_derivative_P2C_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    d1(1) = dm%yp(2) - dm%yp(1)
    d1(3) = dm%yp(3) - dm%yp(2)
    d1(2) = dm%yp(dm%np(2))   - dm%yp(dm%np(2)-1)
    d1(4) = dm%yp(dm%np(2)-1) - dm%yp(dm%np(2)-2)
    call Prepare_TDMA_1deri_P2C_RHS_array(fi(:), fo(:), nsz, d1rP2C(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), dm%h1r(2), ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
          ad1y_P2C(:, ibc(1), ibc(2), iacc), &
          bd1y_P2C(:, ibc(1), ibc(2), iacc), &
          cd1y_P2C(:, ibc(1), ibc(2), iacc), &
          dd1y_P2C(:, ibc(1), ibc(2), iacc), &
          nsz)
    end if

    if(dm%is_stretching(2)) fo(:) = fo(:) * dm%yMappingcc(:, 1)

    return
  end subroutine Get_y_1st_derivative_P2C_1D
!==========================================================================================================
! z - Get_1st_derivative_1D
!==========================================================================================================
  subroutine Get_z_1st_derivative_C2C_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_z_1st_derivative_C2C_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_DIRICHLET  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_DIRICHLET, degragded to IBC_INTRPL. @Get_z_1st_derivative_C2C_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_z_1st_derivative_C2C_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    d1(:) = dm%h(3)
    call Prepare_TDMA_1deri_C2C_RHS_array(fi(:), fo(:), nsz, d1rC2C(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), dm%h1r(3), ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
          ad1z_C2C(:, ibc(1), ibc(2), iacc), &
          bd1z_C2C(:, ibc(1), ibc(2), iacc), &
          cd1z_C2C(:, ibc(1), ibc(2), iacc), &
          dd1z_C2C(:, ibc(1), ibc(2), iacc), &
          nsz)
    end if

    return
  end subroutine Get_z_1st_derivative_C2C_1D
!==========================================================================================================
  subroutine Get_z_1st_derivative_P2P_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_z_1st_derivative_P2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_z_1st_derivative_P2P_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    d1(:) = dm%h(3)
    call Prepare_TDMA_1deri_P2P_RHS_array(fi(:), fo(:), nsz, d1rP2P(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), dm%h1r(3), ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
          ad1z_P2P(:, ibc(1), ibc(2), iacc), &
          bd1z_P2P(:, ibc(1), ibc(2), iacc), &
          cd1z_P2P(:, ibc(1), ibc(2), iacc), &
          dd1z_P2P(:, ibc(1), ibc(2), iacc), &
          nsz)
      end if

    return
  end subroutine Get_z_1st_derivative_P2P_1D
!==========================================================================================================
  subroutine Get_z_1st_derivative_C2P_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_z_1st_derivative_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_DIRICHLET  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_DIRICHLET, degragded to IBC_INTRPL. @Get_z_1st_derivative_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_z_1st_derivative_C2P_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    d1(:) = dm%h(3)
    call Prepare_TDMA_1deri_C2P_RHS_array(fi(:), fo(:), nsz, d1rC2P(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), dm%h1r(3), ibc(:), fbc(:) )
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
          ad1z_C2P(:, ibc(1), ibc(2), iacc), &
          bd1z_C2P(:, ibc(1), ibc(2), iacc), &
          cd1z_C2P(:, ibc(1), ibc(2), iacc), &
          dd1z_C2P(:, ibc(1), ibc(2), iacc), &
          nsz)
    end if

    return
  end subroutine Get_z_1st_derivative_C2P_1D
!==========================================================================================================
  subroutine Get_z_1st_derivative_P2C_1D (fi, fo, dm, iacc, ibc0, fbc)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    
    real(WP),           intent(in ) :: fi(:)
    real(WP),           intent(out) :: fo(:)
    type(t_domain),     intent(in ) :: dm
    integer,            intent(in ) :: iacc
    integer,            intent(in ) :: ibc0(2)
    real(WP), optional, intent(in ) :: fbc(4)

    integer :: nsz
    integer :: i
    integer :: ibc(2)
    real(WP) :: d1(4)
    logical :: is_periodic
    
    ibc = ibc0
    do i = 1, 2
      if (ibc(i) == IBC_INTERIOR  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_INTERIOR, degragded to IBC_INTRPL. @Get_z_1st_derivative_P2C_1D')
        ibc(i) = IBC_INTRPL
      end if
      if (ibc(i) == IBC_NEUMANN  .and. (.not. present(fbc) )) then
        call Print_warning_msg('Lack of fbc info for IBC_NEUMANN, degragded to IBC_INTRPL. @Get_z_1st_derivative_P2C_1D')
        ibc(i) = IBC_INTRPL
      end if
    end do

    nsz = size(fo)
    fo = ZERO
    d1(:) = dm%h(3)
    call Prepare_TDMA_1deri_P2C_RHS_array(fi(:), fo(:), nsz, d1rP2C(1:NL, 1:NS, NBCS:NBCE, iacc), d1(:), dm%h1r(3), ibc(:), fbc(:))
    if (iacc == IACCU_CP4 .or. iacc == IACCU_CP6) then 
      if(ibc(1) == IBC_PERIODIC) then
        is_periodic = .true.
      else
        is_periodic = .false.
      end if
      call Solve_TDMA(is_periodic, fo(:), &
          ad1z_P2C(:, ibc(1), ibc(2), iacc), &
          bd1z_P2C(:, ibc(1), ibc(2), iacc), &
          cd1z_P2C(:, ibc(1), ibc(2), iacc), &
          dd1z_P2C(:, ibc(1), ibc(2), iacc), &
          nsz)
    end if

    return
  end subroutine Get_z_1st_derivative_P2C_1D
!==========================================================================================================
!> \brief To caculate the mid-point interpolation in 3D.
!---------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   pubic
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixsub         x-subdomain index
!> \param[in]     ibc           bc type
!> \param[in]     fbc           bc value
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!==========================================================================================================
  subroutine Get_x_midp_C2P_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use tridiagonal_matrix_algorithm
    use udf_type_mod
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)
    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    real(WP)   :: fbc(4)
    integer :: k, j

!----------------------------------------------------------------------------------------------------------
!  default : x-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(1:4, j, k)
          call Get_x_midp_C2P_1D (fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_x_midp_C2P_1D (fi, fo, dm, iacc, ibc)
        end if
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_midp_C2P_3D
!==========================================================================================================
  subroutine Get_x_midp_P2C_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :) !2 layer each side
    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    integer :: k, j
    real(WP)   :: fbc(4)
!----------------------------------------------------------------------------------------------------------
!  default : x-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(1:4, j, k)
          call Get_x_midp_P2C_1D (fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_x_midp_P2C_1D (fi, fo, dm, iacc, ibc)
        end if
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_midp_P2C_3D
!==========================================================================================================
  subroutine Get_y_midp_C2P_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)
    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    real(WP)   :: fbc(4)
    integer :: k, i

!----------------------------------------------------------------------------------------------------------
!  default : y-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(i, 1:4, k)
          call Get_y_midp_C2P_1D (fi, fo, dm, iacc, ibc, fbc)
        else 
          call Get_y_midp_C2P_1D (fi, fo, dm, iacc, ibc)
        end if
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_midp_C2P_3D
!==========================================================================================================
  subroutine Get_y_midp_P2C_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)
    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    integer :: k, i
    real(WP) :: fbc(4)
!----------------------------------------------------------------------------------------------------------
!  y-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        if(present(fbc2d)) then 
          fbc(1:4) = fbc2d(i, 1:4, k)
          call Get_y_midp_P2C_1D (fi, fo, dm, iacc, ibc, fbc)
        else 
          call Get_y_midp_P2C_1D (fi, fo, dm, iacc, ibc)
        end if
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_midp_P2C_3D
  !==========================================================================================================
  subroutine Get_z_midp_C2P_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)
    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    real(WP)   :: fbc(4)
    integer :: j, i

!----------------------------------------------------------------------------------------------------------
!  default : z-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(i, j, 1:4)
          call Get_z_midp_C2P_1D (fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_z_midp_C2P_1D (fi, fo, dm, iacc, ibc)
        end if
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_midp_C2P_3D
!==========================================================================================================
  subroutine Get_z_midp_P2C_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)
    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    integer :: j, i
    real(WP)   :: fbc(4)
!----------------------------------------------------------------------------------------------------------
!  default : z-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(i, j, 1:4)
          call Get_z_midp_P2C_1D (fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_z_midp_P2C_1D (fi, fo, dm, iacc, ibc)
        end if
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_midp_P2C_3D
!==========================================================================================================
!> \brief To caculate the 1st-deriviate in 3D.
!---------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi            called-freq    xdomain     module
!>       in-given pencil    needed       specified   pubic
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixsub         x-subdomain index
!> \param[in]     ibc           bc type
!> \param[in]     fbc           bc value
!> \param[in]     inbr          the neibouring index of 4 bc nodes
!> \param[in]     fi            the input array of original variable
!> \param[out]    fo            the output array of interpolated variable
!==========================================================================================================
  subroutine Get_x_1st_derivative_C2C_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    real(WP)   :: fbc(4)
    integer :: k, j
!----------------------------------------------------------------------------------------------------------
!  x-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(1:4, j, k)
          call Get_x_1st_derivative_C2C_1D(fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_x_1st_derivative_C2C_1D(fi, fo, dm, iacc, ibc)
        end if
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_1st_derivative_C2C_3D
!==========================================================================================================
  subroutine Get_x_1st_derivative_P2P_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    real(WP)   :: fbc(4)
    integer :: k, j

!----------------------------------------------------------------------------------------------------------
!  x-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        if(present(fbc2d)) then 
          fbc(1:4) = fbc2d(1:4, j, k)
          call Get_x_1st_derivative_P2P_1D(fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_x_1st_derivative_P2P_1D(fi, fo, dm, iacc, ibc)
        end if
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_1st_derivative_P2P_3D
!==========================================================================================================
  subroutine Get_x_1st_derivative_C2P_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)
    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    real(WP)   :: fbc(4)
    integer :: k, j

!----------------------------------------------------------------------------------------------------------
!  x-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        if(present(fbc2d)) then 
          fbc(1:4) = fbc2d(1:4, j, k)
          call Get_x_1st_derivative_C2P_1D(fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_x_1st_derivative_C2P_1D(fi, fo, dm, iacc, ibc)
        end if
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_1st_derivative_C2P_3D
!==========================================================================================================
  subroutine Get_x_1st_derivative_P2C_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)

    real(WP)   :: fi( size(fi3d, 1) )
    real(WP)   :: fo( size(fo3d, 1) )
    real(WP)   :: fbc(4)
    integer :: k, j
!----------------------------------------------------------------------------------------------------------
!  x-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do j = 1, size(fi3d, 2)
        fi(:) = fi3d(:, j, k)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(1:4, j, k)
          call Get_x_1st_derivative_P2C_1D(fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_x_1st_derivative_P2C_1D(fi, fo, dm, iacc, ibc)
        end if
        fo3d(:, j, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_x_1st_derivative_P2C_3D
!==========================================================================================================
  subroutine Get_y_1st_derivative_C2C_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    real(WP)   :: fbc(4)
    integer :: k, i

    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(i, 1:4, k)
          call Get_y_1st_derivative_C2C_1D(fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_y_1st_derivative_C2C_1D(fi, fo, dm, iacc, ibc)
        end if
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_1st_derivative_C2C_3D
!==========================================================================================================
  subroutine Get_y_1st_derivative_P2P_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)
    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    real(WP)   :: fbc(4)
    integer :: k, i

!!----------------------------------------------------------------------------------------------------------
!  y-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        if(present(fbc2d)) then 
          fbc(1:4) = fbc2d(i, 1:4, k)
          call Get_y_1st_derivative_P2P_1D(fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_y_1st_derivative_P2P_1D(fi, fo, dm, iacc, ibc)
        end if
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_1st_derivative_P2P_3D
!==========================================================================================================
  subroutine Get_y_1st_derivative_C2P_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)
    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    real(WP)   :: fbc(4)
    integer :: k, i
!----------------------------------------------------------------------------------------------------------
!  y-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(i, 1:4, k)
          call Get_y_1st_derivative_C2P_1D(fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_y_1st_derivative_C2P_1D(fi, fo, dm, iacc, ibc)
        end if
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_1st_derivative_C2P_3D

!==========================================================================================================
  subroutine Get_y_1st_derivative_P2C_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)

    real(WP)   :: fi( size(fi3d, 2) )
    real(WP)   :: fo( size(fo3d, 2) )
    real(WP)   :: fbc(4)
    integer :: k, i
!----------------------------------------------------------------------------------------------------------
!  y-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do k = 1, size(fi3d, 3)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, :, k)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(i, 1:4, k)
          call Get_y_1st_derivative_P2C_1D(fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_y_1st_derivative_P2C_1D(fi, fo, dm, iacc, ibc)
        end if
        fo3d(i, :, k) = fo(:)
      end do
    end do

    return 
  end subroutine Get_y_1st_derivative_P2C_3D
!==========================================================================================================
  subroutine Get_z_1st_derivative_C2C_3D (fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)
    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    real(WP)   :: fbc(4)
    integer :: j, i
!----------------------------------------------------------------------------------------------------------
!  z-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        if(present(fbc2d)) then 
          fbc(1:4) = fbc2d(i, j, 1:4)
          call Get_z_1st_derivative_C2C_1D(fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_z_1st_derivative_C2C_1D(fi, fo, dm, iacc, ibc)
        end if
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_1st_derivative_C2C_3D

!==========================================================================================================
  subroutine Get_z_1st_derivative_P2P_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)
    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    real(WP)   :: fbc(4)
    integer :: j, i

!----------------------------------------------------------------------------------------------------------
!  z-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(i, j, 1:4)
          call Get_z_1st_derivative_P2P_1D(fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_z_1st_derivative_P2P_1D(fi, fo, dm, iacc, ibc)
        end if
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_1st_derivative_P2P_3D
!==========================================================================================================
  subroutine Get_z_1st_derivative_C2P_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)
    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    real(WP)   :: fbc(4)
    integer :: j, i

!----------------------------------------------------------------------------------------------------------
!  z-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(i, j, 1:4)
          call Get_z_1st_derivative_C2P_1D(fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_z_1st_derivative_C2P_1D(fi, fo, dm, iacc, ibc)
        end if
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_1st_derivative_C2P_3D
  !==========================================================================================================
  subroutine Get_z_1st_derivative_P2C_3D(fi3d, fo3d, dm, iacc, ibc, fbc2d)
    use parameters_constant_mod
    use udf_type_mod
    use tridiagonal_matrix_algorithm
    implicit none
    real(WP),           intent(in) :: fi3d(:, :, :)
    real(WP),           intent(out):: fo3d(:, :, :)
    type(t_domain),     intent(in) :: dm
    integer,            intent(in) :: iacc
    integer,            intent(in) :: ibc(2)
    real(WP), optional, intent(in) :: fbc2d(:, :, :)

    real(WP)   :: fi( size(fi3d, 3) )
    real(WP)   :: fo( size(fo3d, 3) )
    real(WP)   :: fbc(4)
    integer :: j, i
!----------------------------------------------------------------------------------------------------------
!  z-pencil calculation
!----------------------------------------------------------------------------------------------------------
    fo3d(:, :, :) = ZERO
    do j = 1, size(fi3d, 2)
      do i = 1, size(fi3d, 1)
        fi(:) = fi3d(i, j, :)
        if(present(fbc2d)) then
          fbc(1:4) = fbc2d(i, j, 1:4)
          call Get_z_1st_derivative_P2C_1D(fi, fo, dm, iacc, ibc, fbc)
        else
          call Get_z_1st_derivative_P2C_1D(fi, fo, dm, iacc, ibc)
        end if
        fo3d(i, j, :) = fo(:)
      end do
    end do

    return 
  end subroutine Get_z_1st_derivative_P2C_3D

!==========================================================================================================
!==========================================================================================================
!> \brief To test this subroutine for mid-point interpolation.
!>
!> This subroutine is called in \ref Test_algorithms. Define the logicals to choose
!> which test section is required. 
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!_______________________________________________________________________________
  subroutine test_function_setup(ibc, dd, fbc, scale, shift)
    use parameters_constant_mod
    use math_mod
    implicit none
    integer, intent(in) :: ibc(2)
    real(WP), intent(in) :: dd
    real(WP), intent(out) :: fbc(4)
    real(WP), intent(out) :: scale
    real(WP), intent(out) :: shift

    integer :: i
    fbc = MAXP
    do i = 1, 2
      if (ibc(i) == IBC_PERIODIC) then
        ! f = sin(x/scale + shift) : f = sin(x)
        scale = ONE
        shift = ZERO
      else if (ibc(i) == IBC_SYMMETRIC) then
      ! f = sin(x/scale + shift) : f = sin(x + pi/2)
        scale = ONE
        shift = PI * HALF
      else if (ibc(i) == IBC_ASYMMETRIC) then
      ! f = sin(x/scale + shift) : f = sin(x/2)
        scale = TWO
        shift = ZERO
      else if (ibc(i) == IBC_DIRICHLET) then
      ! f = sin(x/scale + shift) : f = sin(x/2)
        scale = TWO
        shift = ZERO
        if(i==1) fbc(1) = ZERO
        if(i==2) fbc(2) = sin_wp(TWOPI / TWO)
      else if (ibc(i) == IBC_NEUMANN) then
      ! f = sin(x/scale + shift) : f = sin(x/2) : f'=1/2cos(x/2)
        scale = TWO
        shift = ZERO
        if(i==1) fbc(1) = HALF * cos_wp(ZERO )
        if(i==2) fbc(2) = HALF * cos_wp(TWOPI * HALF)
      else 
        scale = TWO
        shift = ZERO
        if(i==1) fbc(1) = ZERO
        if(i==2) fbc(2) = sin_wp(TWOPI * HALF)
      end if
    end do
    return
  end subroutine 

  subroutine test_interp_c2p_comparison(nc, np, dd, scale, shift, iacc, ibc, fbc, dm, str)
    use parameters_constant_mod
    use math_mod
    use udf_type_mod
    implicit none
    integer, intent(in) :: nc
    integer, intent(in) :: np
    real(WP), intent(in) :: dd
    real(WP), intent(in) :: scale
    real(WP), intent(in) :: shift
    integer, intent(in) :: iacc
    integer, intent(in) :: ibc(2)
    real(WP), intent(in) :: fbc(4)
    character(1), intent(in) :: str
    type(t_domain), intent(in) :: dm

    integer :: i
    real(WP) :: xc, xp
    real(WP) :: ref
    real(WP) :: err, err_Linf, err_L2
    integer  :: wrt_unit(2)
    real(WP) :: fgxp(np), fxc(nc)

    open (newunit = wrt_unit(1), file = 'test_interp.dat', position="append")
    open (newunit = wrt_unit(2), file = 'test_interp_'//trim(str)//'.dat', position="append")
    
  ! x direction, xc
    do i = 1, nc
      xc = dd * (real(i - 1, WP) + HALF)
      fxc(i) = sin_wp ( xc / scale + shift)
    end do
  ! x: c2p
    if(trim(str)=='x') then
      call Get_x_midp_C2P_1D (fxc, fgxp, dm, iacc, ibc, fbc)
    else if(trim(str)=='y') then
      call Get_y_midp_C2P_1D (fxc, fgxp, dm, iacc, ibc, fbc)
    else if(trim(str)=='z') then
      call Get_z_midp_C2P_1D (fxc, fgxp, dm, iacc, ibc, fbc)
    else
    end if
  ! x: error 
    err_Linf = ZERO
    err_L2   = ZERO
    write(wrt_unit(2), *) '# interp-c2p-'//trim(str), ', iacc=', iacc, ', np=', np
    write(wrt_unit(2), *) '# i, xp, ref, cal, err'
    do i = 1, np
      xp = dd * real(i - 1, WP)
      ref = sin_wp(xp / scale + shift)
      err = abs_wp(fgxp(i) - ref)
      write(wrt_unit(2), *) i, xp, ref, fgxp(i), err !test
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
    end do
    err_L2 = sqrt_wp(err_L2 / np) 
    write(wrt_unit(1), *) '# interp-c2p-'//trim(str), &
                           ', iacc=', iacc, ', np=', np, &
                           ', ibc=', ibc, ', eInf=', err_Linf, ', eL2=', err_L2
    
    close(wrt_unit(1))
    close(wrt_unit(2))
    return
  end subroutine

  subroutine test_interp_p2c_comparison(nc, np, dd, scale, shift, iacc,  ibc, fbc, dm, str)
    use parameters_constant_mod
    use math_mod
    use udf_type_mod
    implicit none
    integer, intent(in) :: nc
    integer, intent(in) :: np
    real(WP), intent(in) :: dd
    real(WP), intent(in) :: scale
    real(WP), intent(in) :: shift
    integer, intent(in) :: iacc
    integer, intent(in) :: ibc(2)
    real(WP), intent(in) :: fbc(4)
    character(1), intent(in) :: str
    type(t_domain), intent(in) :: dm

    integer :: i
    real(WP) :: xc, xp
    real(WP) :: ref
    real(WP) :: err, err_Linf, err_L2
    integer  :: wrt_unit(2)
    real(WP) :: fgxc(nc)
    real(WP) :: fxp(np)

    open (newunit = wrt_unit(1), file = 'test_interp.dat', position="append")
    open (newunit = wrt_unit(2), file = 'test_interp_'//trim(str)//'.dat', position="append")

    ! x direction, xp
    do i = 1, np
      xp = dd * real(i - 1, WP)
      fxp(i) = sin_wp ( xp / scale + shift)
    end do
! x: p2c
    if(trim(str)=='x') then
      call Get_x_midp_P2C_1D (fxp, fgxc, dm, iacc, ibc, fbc)
    else if(trim(str)=='y') then
      call Get_y_midp_P2C_1D (fxp, fgxc, dm, iacc, ibc, fbc)
    else if(trim(str)=='z') then
      call Get_z_midp_P2C_1D (fxp, fgxc, dm, iacc, ibc, fbc)
    else
    end if

    err_Linf = ZERO
    err_L2   = ZERO
    write(wrt_unit(2), *) '# interp-p2c-'//trim(str), ', iacc=', iacc, ', nc=', nc
    write(wrt_unit(2), *) '# i, xc, ref, cal, err'
    do i = 1, nc
      xc = dd * (real(i - 1, WP) + HALF)
      ref = sin_wp(xc / scale + shift)
      err = abs_wp(fgxc(i) - ref)
      write(wrt_unit(2),*) i, xc, ref, fgxc(i), err !test
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
    end do
    err_L2 = sqrt_wp(err_L2 / nc) 
    write(wrt_unit(1), *) '# interp-p2c-'//trim(str), &
                           ', iacc=', iacc, ', nc=', nc, &
                           ', ibc=', ibc, ', eInf=', err_Linf, ', eL2=', err_L2
    
    close(wrt_unit(1))
    close(wrt_unit(2))
    return
  end subroutine

  subroutine Test_interpolation(dm)
    use parameters_constant_mod
    use math_mod
    use EvenOdd_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm

    real(WP) :: scale, shift
    real(WP) :: fbcx(4), fbcy(4), fbcz(4), fbc(4)
    integer  :: ibcx(2), ibcy(2), ibcz(2), ibc(2)
    integer :: iacc, n, i
    character(1) :: str

    dm%h(1) = TWOPI / dm%nc(1)
    dm%h(2) = TWOPI / dm%nc(2)
    dm%h(3) = TWOPI / dm%nc(3)
    dm%h1r(1) = ONE / dm%h(1)
    dm%h1r(2) = ONE / dm%h(2)
    dm%h1r(3) = ONE / dm%h(3)
    ibcx(:) = dm%ibcx_nominal(:, 5)
    ibcy(:) = dm%ibcy_nominal(:, 5)
    ibcz(:) = dm%ibcz_nominal(:, 5)


    do i = 1, 3
      if (i == 1) then
        ibc(:) = ibcx(:)
        fbc(:) = fbcx(:)
        str = 'x'
      else if (i == 2) then
        ibc(:) = ibcy(:)
        fbc(:) = fbcy(:)
        str = 'y'
      else if (i == 3) then
        ibc(:) = ibcz(:)
        fbc(:) = fbcz(:)
        str = 'z'
      else 
      end if
      call test_function_setup(ibc, dm%h(i), fbc, scale, shift)
      do n = 1, NACC
        iacc = n
        call test_interp_p2c_comparison(dm%nc(i), dm%np(i), dm%h(i), scale, shift, iacc, ibc, fbc, dm, str)
      end do
      do n = 1, NACC
        iacc = n
        call test_interp_c2p_comparison(dm%nc(i), dm%np(i), dm%h(i), scale, shift, iacc, ibc, fbc, dm, str)
      end do
    end do 

    return 
  end subroutine

!==========================================================================================================
!==========================================================================================================
!> \brief To test this subroutine for mid-point interpolation.
!>
!> This subroutine is called in \ref Test_algorithms. Define the logicals to choose
!> which test section is required. 
!>
!----------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!_______________________________________________________________________________
  subroutine test_1stder_p2c_comparison(nc, np, dd, scale, shift, iacc,  ibc, fbc, dm, str)
    use parameters_constant_mod
    use math_mod
    use udf_type_mod
    implicit none
    integer, intent(in) :: nc
    integer, intent(in) :: np
    real(WP), intent(in) :: dd
    real(WP), intent(in) :: scale
    real(WP), intent(in) :: shift
    integer, intent(in) :: iacc
    integer, intent(in) :: ibc(2)
    real(WP), intent(in) :: fbc(4)
    character(1), intent(in) :: str
    type(t_domain), intent(in) :: dm

    integer :: i
    real(WP) :: xc, xp
    real(WP) :: ref
    real(WP) :: err, err_Linf, err_L2
    integer  :: wrt_unit(2)
    real(WP) :: fgxc(nc)
    real(WP) :: fxp(np)

    open (newunit = wrt_unit(1), file = 'test_1stder.dat', position="append")
    open (newunit = wrt_unit(2), file = 'test_1stder_'//trim(str)//'.dat', position="append")

    do i = 1, np
      xp = dd * real(i - 1, WP)
      fxp(i) = sin_wp ( xp / scale + shift)
    end do

    if(trim(str)=='x') then
      call Get_x_1st_derivative_P2C_1D(fxp, fgxc, dm, iacc, ibc, fbc)
    else if(trim(str)=='y') then
      call Get_y_1st_derivative_P2C_1D(fxp, fgxc, dm, iacc, ibc, fbc)
    else if(trim(str)=='z') then
      call Get_z_1st_derivative_P2C_1D(fxp, fgxc, dm, iacc, ibc, fbc)
    else
    end if
    
    err_Linf = ZERO
    err_L2   = ZERO
    write(wrt_unit(2), *) '# 1stder-p2c-'//trim(str), ', iacc=', iacc, ', nc=', nc
    write(wrt_unit(2), *) '# i, xc, ref, cal, err'
    do i = 1, nc
      xc = dd * (real(i - 1, WP) + HALF)
      ref = ONE/scale * cos_wp(xc / scale + shift)
      err = abs_wp(fgxc(i) - ref)
      write(wrt_unit(2),*) i, xc, ref, fgxc(i), err
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
    end do
    err_L2 = sqrt_wp(err_L2 / nc) 
    write(wrt_unit(1), *) '# 1stder-p2c-'//trim(str), &
                           ', iacc=', iacc, ', nc=', nc, &
                           ', ibc=', ibc, ', eInf=', err_Linf, ', eL2=', err_L2

    close(wrt_unit(1))
    close(wrt_unit(2))
    
    return
  end subroutine

  subroutine test_1stder_p2p_comparison(nc, np, dd, scale, shift, iacc,  ibc, fbc, dm, str)
    use parameters_constant_mod
    use math_mod
    use udf_type_mod
    implicit none
    integer, intent(in) :: nc
    integer, intent(in) :: np
    real(WP), intent(in) :: dd
    real(WP), intent(in) :: scale
    real(WP), intent(in) :: shift
    integer, intent(in) :: iacc
    integer, intent(in) :: ibc(2)
    real(WP), intent(in) :: fbc(4)
    character(1), intent(in) :: str
    type(t_domain), intent(in) :: dm

    integer :: i
    real(WP) :: xc, xp
    real(WP) :: ref
    real(WP) :: err, err_Linf, err_L2
    integer  :: wrt_unit(2)
    real(WP) :: fgxp(np)
    real(WP) :: fxp(np)

    open (newunit = wrt_unit(1), file = 'test_1stder.dat', position="append")
    open (newunit = wrt_unit(2), file = 'test_1stder_'//trim(str)//'.dat', position="append")


    do i = 1, np
      xp = dd * real(i - 1, WP)
      fxp(i) = sin_wp ( xp / scale + shift)
    end do

    if(trim(str)=='x') then
      call Get_x_1st_derivative_P2P_1D(fxp, fgxp, dm, iacc, ibc, fbc)
    else if(trim(str)=='y') then
      call Get_y_1st_derivative_P2P_1D(fxp, fgxp, dm, iacc, ibc, fbc)
    else if(trim(str)=='z') then
      call Get_z_1st_derivative_P2P_1D(fxp, fgxp, dm, iacc, ibc, fbc)
    else
    end if

    err_Linf = ZERO
    err_L2   = ZERO
    write(wrt_unit(2), *) '# 1stder-p2p-'//trim(str),  ', iacc=', iacc, ', np=', np
    write(wrt_unit(2), *) '# i, xp, ref, cal, err'
    do i = 1, np
      xp = dd * real(i - 1, WP)
      ref = ONE/scale * cos_wp(xp / scale + shift)
      err = abs_wp(fgxp(i) - ref)
      write(wrt_unit(2),*) i, xp, ref, fgxp(i), err
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
    end do
    err_L2 = sqrt_wp(err_L2 / nc) 
    write(wrt_unit(1), *) '# 1stder-p2p-'//trim(str), &
                           ', iacc=', iacc, ', np=', np, &
                           ', ibc=', ibc, ', eInf=', err_Linf, ', eL2=', err_L2
    
    close(wrt_unit(1))
    close(wrt_unit(2))

    return
  end subroutine

  subroutine test_1stder_c2p_comparison(nc, np, dd, scale, shift, iacc,  ibc, fbc, dm, str)
    use parameters_constant_mod
    use math_mod
    use udf_type_mod
    implicit none
    integer, intent(in) :: nc
    integer, intent(in) :: np
    real(WP), intent(in) :: dd
    real(WP), intent(in) :: scale
    real(WP), intent(in) :: shift
    integer, intent(in) :: iacc
    integer, intent(in) :: ibc(2)
    real(WP), intent(in) :: fbc(4)
    character(1), intent(in) :: str
    type(t_domain), intent(in) :: dm

    integer :: i
    real(WP) :: xc, xp
    real(WP) :: ref
    real(WP) :: err, err_Linf, err_L2
    integer  :: wrt_unit(2)
    real(WP) :: fgxp(np)
    real(WP) :: fxc(nc)

    open (newunit = wrt_unit(1), file = 'test_1stder.dat', position="append")
    open (newunit = wrt_unit(2), file = 'test_1stder_'//trim(str)//'.dat', position="append")

    do i = 1, nc
      xc =  dm%h(1) * (real(i - 1, WP) + HALF)
      fxc(i) = sin_wp ( xc / scale + shift)
    end do

    
    if(trim(str)=='x') then
      call Get_x_1st_derivative_C2P_1D(fxc, fgxp, dm, iacc, ibc, fbc)
    else if(trim(str)=='y') then
      call Get_y_1st_derivative_C2P_1D(fxc, fgxp, dm, iacc, ibc, fbc)
    else if(trim(str)=='z') then
      call Get_z_1st_derivative_C2P_1D(fxc, fgxp, dm, iacc, ibc, fbc)
    else
    end if

    err_Linf = ZERO
    err_L2   = ZERO
    write(wrt_unit(2), *) '# 1stder-c2p-'//trim(str), ', iacc=', iacc, ', np=', np
    write(wrt_unit(2), *) '# i, xp, ref, cal, err'
    do i = 1, np
      xp = dd * real(i - 1, WP)
      ref = ONE/scale * cos_wp(xp / scale + shift)
      err = abs_wp(fgxp(i) - ref)
      write(wrt_unit(2),*) i, xp, ref, fgxp(i), err
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
    end do
    err_L2 = sqrt_wp(err_L2 / nc) 
    write(wrt_unit(1), *) '# 1stder-c2p-'//trim(str), &
                           ', iacc=', iacc, ', np=', np, &
                           ', ibc=', ibc, ', eInf=', err_Linf, ', eL2=', err_L2
    
    close(wrt_unit(1))
    close(wrt_unit(2))

    return
  end subroutine

  subroutine test_1stder_c2c_comparison(nc, np, dd, scale, shift, iacc,  ibc, fbc, dm, str)
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    integer, intent(in) :: nc
    integer, intent(in) :: np
    real(WP), intent(in) :: dd
    real(WP), intent(in) :: scale
    real(WP), intent(in) :: shift
    integer, intent(in) :: iacc
    integer, intent(in) :: ibc(2)
    real(WP), intent(in) :: fbc(4)
    character(1), intent(in) :: str
    type(t_domain), intent(in) :: dm

    integer :: i
    real(WP) :: xc, xp
    real(WP) :: ref
    real(WP) :: err, err_Linf, err_L2
    integer  :: wrt_unit(2)
    real(WP) :: fgxc(nc)
    real(WP) :: fxc(nc)

    open (newunit = wrt_unit(1), file = 'test_1stder.dat', position="append")
    open (newunit = wrt_unit(2), file = 'test_1stder_'//trim(str)//'.dat', position="append")


    do i = 1, nc
      xc =  dd * (real(i - 1, WP) + HALF)
      fxc(i) = sin_wp ( xc / scale + shift)
    end do

    
    if(trim(str)=='x') then
      call Get_x_1st_derivative_C2C_1D(fxc, fgxc, dm, iacc, ibc, fbc)
    else if(trim(str)=='y') then
      call Get_y_1st_derivative_C2C_1D(fxc, fgxc, dm, iacc, ibc, fbc)
    else if(trim(str)=='z') then
      call Get_z_1st_derivative_C2C_1D(fxc, fgxc, dm, iacc, ibc, fbc)
    else
    end if


    err_Linf = ZERO
    err_L2   = ZERO
    write(wrt_unit(2), *) '# 1stder-c2c-'//trim(str), ', iacc=', iacc, ', nc=', nc
    write(wrt_unit(2), *) '# i, xc, ref, cal, err'
    do i = 1, nc
      xc =  dd * (real(i - 1, WP) + HALF)
      ref = ONE/scale * cos_wp(xc / scale + shift)
      err = abs_wp(fgxc(i) - ref)
      write(wrt_unit(2),*) i, xc, ref, fgxc(i), err
      if(err > err_Linf) err_Linf = err
      err_L2 = err_L2 + err**2
    end do
    err_L2 = sqrt_wp(err_L2 / nc) 
    write(wrt_unit(1), *) '# 1stder-c2c-'//trim(str), &
                           ', iacc=', iacc, ', nc=', nc, &
                           ', ibc=', ibc, ', eInf=', err_Linf, ', eL2=', err_L2
    
    close(wrt_unit(1))
    close(wrt_unit(2))
    
    return
  end subroutine

  subroutine Test_1st_derivative(dm)
    use parameters_constant_mod
    use EvenOdd_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(inout) :: dm

    real(WP) :: scale, shift
    real(WP) :: fbcx(4), fbcy(4), fbcz(4), fbc(4)
    integer  :: ibcx(2), ibcy(2), ibcz(2), ibc(2)
    integer :: n, iacc, i
    character(1) :: str

    dm%h(1) = TWOPI / dm%nc(1)
    dm%h(2) = TWOPI / dm%nc(2)
    dm%h(3) = TWOPI / dm%nc(3)
    ibcx(:) = dm%ibcx_nominal(:, 5)
    ibcy(:) = dm%ibcy_nominal(:, 5)
    ibcz(:) = dm%ibcz_nominal(:, 5)


    do i = 1, 3
      if (i == 1) then
        ibc(:) = ibcx(:)
        fbc(:) = fbcx(:)
        str = 'x'
      else if (i == 2) then
        ibc(:) = ibcy(:)
        fbc(:) = fbcy(:)
        str = 'y'
      else if (i == 3) then
        ibc(:) = ibcz(:)
        fbc(:) = fbcz(:)
        str = 'z'
      else 
      end if
      call test_function_setup(ibc, dm%h(i), fbc, scale, shift)
      write(*,*) 'input=sin(x/',scale,'+',shift,')'
      do n = 1, NACC
        iacc = n
        call test_1stder_p2p_comparison(dm%nc(i), dm%np(i), dm%h(i), scale, shift, iacc, ibc, fbc, dm, str)
      end do
      do n = 1, NACC
        iacc = n
        call test_1stder_c2c_comparison(dm%nc(i), dm%np(i), dm%h(i), scale, shift, iacc, ibc, fbc, dm, str)
      end do
      do n = 1, NACC
        iacc = n
        call test_1stder_p2c_comparison(dm%nc(i), dm%np(i), dm%h(i), scale, shift, iacc, ibc, fbc, dm, str)
      end do
      do n = 1, NACC
        iacc = n
        call test_1stder_c2p_comparison(dm%nc(i), dm%np(i), dm%h(i), scale, shift, iacc, ibc, fbc, dm, str)
      end do
    end do 
    
    return 
  end subroutine

end module
