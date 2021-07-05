module poisson_mod
  use precision_mod, only : WP
  implicit none
!_______________________________________________________________________________
! Transformation Matrix from \hat{f} to \hat{f''}
!_______________________________________________________________________________
  private 
  real(wp), allocatable, dimension(:)             :: t2x, t2y, t2z
  real(wp), save, allocatable, dimension(:, :, :) :: t2xyz
!_______________________________________________________________________________
! boundary conditions
!_______________________________________________________________________________
  integer(4), save :: bcx, bcy, bcz
  integer(4), save :: nx, ny, nz
!_______________________________________________________________________________
! work arrays, 
! naming convention: cw (complex); rw (real); 
!                    b =     ; c = 
!                    1 = X-pencil; 2 = Y-pencil; 3 = Z-pencil
!_______________________________________________________________________________

  complex(wp), allocatable, dimension(:, :, :) :: cw1
!_______________________________________________________________________________
! FFT library only needs to be initialised once
!_______________________________________________________________________________
  logical, save :: is_fft_initialised = .false.
!_______________________________________________________________________________
! interface for basic poisson solver
!_______________________________________________________________________________
  ABSTRACT INTERFACE
    SUBROUTINE Solve_poisson_xxx(rhs)
      use precision_mod, only : WP
      real(wp), dimension(:, :, :), intent(INOUT) :: rhs
    END SUBROUTINE Solve_poisson_xxx
  END INTERFACE

  PROCEDURE (Solve_poisson_xxx), POINTER ::  Solve_poisson

  public :: Initialize_decomp_poisson, &
            Finalize_decomp_poisson, &
            Solve_poisson

  private :: Transform_2nd_derivative_spectral_1d, &
             Transform_2nd_derivative_spectral_3d

contains
!===============================================================================
!===============================================================================
!> \brief To calcuate all rhs of momentum eq.
!>
!> This subroutine is called everytime when calcuting the rhs of momentum eqs.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  f             flow field
!> \param[inout]  d             domain    
!> \param[in]     isub          the RK iteration to get correct Coefficient 
!_______________________________________________________________________________
  subroutine Initialize_decomp_poisson(d)
    use udf_type_mod,            only : t_domain
    use parameters_constant_mod, only : ZERO
    implicit none
    type(t_domain), intent(in) :: d

    call Print_debug_start_msg("Initializing variables for Poisson Solver ...")
!_______________________________________________________________________________
! set up boundary flags for periodic b.c.
!_______________________________________________________________________________
    if ( d%is_periodic(1) ) then
      bcx = 0
    else
      bcx = 1
    end if

    if ( d%is_periodic(2) ) then
      bcy = 0
    else
      bcy = 1
    end if

    if ( d%is_periodic(3) ) then
      bcz = 0
    else
      bcz = 1
    end if
!_______________________________________________________________________________
! Top level wrapper
! Note: if periodic b.c. exsits, it should be z direction first. 
!   x    y     z
!   0    0     0 
!   1    0     0
!   0    1     0
!   0    0     1 (X, not existing)
!   1    1     0
!   1    0     1 (X, not existing)
!   0    1     1 (X, not existing)
!   1    1     1
!_______________________________________________________________________________
    if      (bcx == 0 .and. bcy == 0 .and. bcz == 0) then
      Solve_poisson => Solve_poisson_000
    else if (bcx == 1 .and. bcy == 0 .and. bcz == 0) then
      !Solve_poisson => Solve_poisson_100
    else if (bcx == 0 .and. bcy == 1 .and. bcz == 0) then
      !Solve_poisson => Solve_poisson_010
    else if (bcx == 1 .and. bcy == 1) then   ! 110 & 111
      !Solve_poisson => Solve_poisson_11x
    else
      stop 'boundary condition not supported'
    end if
!_______________________________________________________________________________
! size of working array
! pressure-grid having 1 fewer point for non-periodic directions
!_______________________________________________________________________________
    nx = d%nc(1)
    ny = d%nc(2)
    nz = d%nc(3)

    if (bcx == 1) nx = nx + 1 ! to check ?
    if (bcy == 1) ny = ny + 1
    if (bcz == 1) nz = nz + 1
!_______________________________________________________________________________
! allocate working arrays
! complex and real part of working array in x-pencil
!_______________________________________________________________________________
    allocate (cw1 (1 : nx, 1 : ny, 1 : nz ) )
!_______________________________________________________________________________
! prepare the transformation \hat{f''}_l = \hat{f}_l * t2x
!_______________________________________________________________________________
    allocate (t2x(nx)) ;  t2x = ZERO
    allocate (t2y(ny)) ;  t2y = ZERO
    allocate (t2z(nz)) ;  t2z = ZERO

    call Transform_2nd_derivative_spectral_1d(d%bc(1, 1), d%nc(1), d%h(1), t2x)
    call Transform_2nd_derivative_spectral_1d(d%bc(1, 2), d%nc(2), d%h(2), t2y)
    call Transform_2nd_derivative_spectral_1d(d%bc(1, 3), d%nc(3), d%h(3), t2z)

    allocate (t2xyz(nx, ny, nz)) ;  t2xyz = ZERO
    call Transform_2nd_derivative_spectral_3d(nx, ny, nz, t2x, t2y, t2z, t2xyz)

    deallocate(t2x)
    deallocate(t2y)
    deallocate(t2z)

    call Print_debug_end_msg
    return
  end subroutine Initialize_decomp_poisson
!===============================================================================
!===============================================================================
  subroutine Finalize_decomp_poisson
    use decomp_2d_fft
    implicit none

    call decomp_2d_fft_finalize
    is_fft_initialised = .false.

    deallocate(t2xyz)
    deallocate(cw1)

    return
  end subroutine Finalize_decomp_poisson
!===============================================================================
!===============================================================================
  subroutine Transform_2nd_derivative_spectral_1d(ibc, nn, dd, t2)
    use operations,              only : d2fC2C, d2rC2C
    use input_general_mod,       only : IBC_PERIODIC
    use parameters_constant_mod, only : FOUR, TWO, ONE, PI
    use math_mod,                only : cos_wp
    implicit none

    integer(4), intent(in)    :: nn
    integer(4), intent(in)    :: ibc
    real(wp),   intent(in)    :: dd
    real(wp),   intent(inout) :: t2(:)

    real(wp) :: a, b, alpha
    real(wp) :: w, cosw
    integer(4) :: i

    if(ibc == IBC_PERIODIC) then
    ! below is for 3 periodic only. due to enrich
    alpha = d2fC2C(3, 1, ibc)
    a     = d2rC2C(3, 1, ibc)
    b     = d2rC2C(3, 2, ibc) * FOUR

    do i = 1, nn
      w = TWO * PI * REAL(i - 1, WP) / REAL(nn, WP)
      cosw = cos_wp(w)
      t2(i) = b * cosw * cosw + TWO * a * cosw - TWO * a - b
      t2(i) = t2(i) / (ONE + TWO * alpha * cosw) / dd / dd
      !write(*, *) i, t2(i)
    end do
    
    end if
    return
  end subroutine Transform_2nd_derivative_spectral_1d
!===============================================================================
!===============================================================================
  subroutine Transform_2nd_derivative_spectral_3d(nnx, nny, nnz, tt2x, tt2y, tt2z, tt2xyz)
    implicit none
    integer(4), intent(in) :: nnx
    integer(4), intent(in) :: nny
    integer(4), intent(in) :: nnz
    real(WP),   intent(in) :: tt2x(:)
    real(WP),   intent(in) :: tt2y(:)
    real(WP),   intent(in) :: tt2z(:)
    real(WP),   intent(out) :: tt2xyz(:, :, :)
    integer(4) :: i, j, k

    do k = 1, nnz
        do j = 1, nny
            do i = 1, nnx
                tt2xyz(i, j, k) = tt2x(i) + tt2y(j) + tt2z(k)
            end do
        end do
    end do
    return 
  end subroutine
!===============================================================================
!===============================================================================
  subroutine Solve_poisson_000(rhs)
    use decomp_2d_fft!, only : decomp_2d_fft_3d
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT) :: rhs
    integer :: nx, ny, nz, i,j,k

    if (.not. is_fft_initialised) then
      call decomp_2d_fft_init()
      is_fft_initialised = .true.
    end if

    ! compute r2c transform 
    call decomp_2d_fft_3d(rhs,cw1)

    ! normalisation
    cw1 = cw1 / real(nx, kind=wp) /real(ny, kind=wp) /real(nz, kind=wp)

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
            cw1(i,j,k) = cw1(i,j,k) / t2xyz(i, j, k)
          end do
       end do
    end do

    ! compute c2r transform
    call decomp_2d_fft_3d(cw1,rhs)

    !   call decomp_2d_fft_finalize

    return
  end subroutine Solve_poisson_000

end module poisson_mod
