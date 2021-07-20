module poisson_mod
    use precision_mod
    use decomp_2d
    use decomp_2d_fft
    use decomp_2d_io
    use MPI
  
    implicit none
    private
  !_______________________________________________________________________________
  ! Variables for debugging and I/O
  !_______________________________________________________________________________
    CHARACTER(LEN=*), PARAMETER :: cmplxfmt = '(ES13.5,SP,ES13.5,"i")'
  !_______________________________________________________________________________
  ! Transformation Matrix from \hat{f} to \hat{f''}
  !_______________________________________________________________________________
    real(wp),       allocatable, dimension(:)       :: t2x, t2y, t2z
    real(wp), save, allocatable, dimension(:, :, :) :: t2xyz
  !_______________________________________________________________________________
  ! boundary conditions and index
  !_______________________________________________________________________________
    integer(4),               save :: bcx, bcy, bcz
    integer(4),               save :: nx,  ny,  nz
    integer(4), dimension(3), save :: fft_start, fft_end, fft_size
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
  
    private :: Transform_2nd_derivative_spectral_1d
  
    public :: Test_poisson_solver
  
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
      use parameters_constant_mod!, only : ZERO, MAXP, TRUNCERR
      implicit none
      type(t_domain), intent(in) :: d
      integer(4) :: i, j, k
  
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
  ! initialise FFT library
  !_______________________________________________________________________________
      call decomp_2d_fft_init
      is_fft_initialised = .true.
      call decomp_2d_fft_get_size( fft_start, fft_end, fft_size)
  !_______________________________________________________________________________
  ! allocate space for wave-space variables
  !_______________________________________________________________________________
      write(*,*) 'Fourier Domain index, l = ', fft_start(1), fft_end(1)
      write(*,*) 'Fourier Domain index, m = ', fft_start(2), fft_end(2)
      write(*,*) 'Fourier Domain index, n = ', fft_start(3), fft_end(3)
  
      allocate ( cw1( fft_start(1) : fft_end(1), &
                      fft_start(2) : fft_end(2), &
                      fft_start(3) : fft_end(3)) )
  !_______________________________________________________________________________
  ! prepare the transformation \hat{f"}_l = \hat{f}_l * t2x
  !_______________________________________________________________________________
      allocate ( t2x(fft_start(1) : fft_end(1)) ) ;  t2x = ZERO
      allocate ( t2y(fft_start(2) : fft_end(2)) ) ;  t2y = ZERO
      allocate ( t2z(fft_start(3) : fft_end(3)) ) ;  t2z = ZERO
      
      call Transform_2nd_derivative_spectral_1d &
            (d%bc(1, 1), fft_start(1), fft_end(1), nx, d%h(1), t2x)
      call Transform_2nd_derivative_spectral_1d &
            (d%bc(1, 2), fft_start(2), fft_end(2), ny, d%h(2), t2y)
      call Transform_2nd_derivative_spectral_1d &
            (d%bc(1, 3), fft_start(3), fft_end(3), nz, d%h(3), t2z)
  
      allocate (t2xyz( fft_start(1) : fft_end(1), &
                       fft_start(2) : fft_end(2), &
                       fft_start(3) : fft_end(3))) ;  t2xyz = ZERO
  
      do k = fft_start(3) , fft_end(3)
        do j = fft_start(2) , fft_end(2)
          do i = fft_start(1) , fft_end(1)
            t2xyz(i, j, k) = t2x(i) + t2y(j) + t2z(k)
            if(dabs(t2xyz(i, j, k)) < TRUNCERR) then
              t2xyz(i, j, k) = MAXP
            end if
          end do
        end do
      end do
  
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
    subroutine Transform_2nd_derivative_spectral_1d(ibc, n0, n1, nn, dd, t2)
      use operations!,              only : d2fC2C, d2rC2C
      use input_general_mod,       only : IBC_PERIODIC
      use parameters_constant_mod!, only : FOUR, TWO, ONE, PI
      use math_mod!,                only : cos_wp
      implicit none
  
      integer(4), intent(in)    :: n0, n1
      integer(4), intent(in)    :: ibc
      integer(4), intent(in)    :: nn
      real(wp),   intent(in)    :: dd
      real(wp),   intent(inout) :: t2(n0 : n1)
  
      real(wp) :: a, b, alpha
      real(wp) :: w, cosw
      integer(4) :: i, i1
      integer(4) :: method = 2
  
      if(ibc == IBC_PERIODIC) then
      ! below is for 3 periodic only. due to enrich
  
      if(method == 2) then
  ! method 1: 2nd deriviate
  
      alpha = d2fC2C(3, 1, ibc)
      a     = d2rC2C(3, 1, ibc)
      b     = d2rC2C(3, 2, ibc) * FOUR
  
      do i = n0, n1
        w = TWO * PI * REAL(i - 1, WP) / REAL(nn, WP) ! check, for 0-n/2, pi or 2pi?
        cosw = cos_wp(w)
        t2(i) = b * cosw * cosw + TWO * a * cosw - TWO * a - b
        t2(i) = t2(i) / (ONE + TWO * alpha * cosw) / dd / dd
        !write(*, *) i, t2(i)
      end do
  
  
    else if (method == 1) then
  ! method 1: 1st deriviate
      alpha = d1fC2C(3, 1, ibc)
      a     = d1rC2C(3, 1, ibc) * TWO
      b     = d1rC2C(3, 2, ibc) * FOUR
  
      do i = n0, n1
        w = TWO * PI * REAL(i  - 1, WP) / REAL(nn, WP) !
        t2(i) = a * sin_wp(w) + b * HALF * sin_wp(TWO * w)
        t2(i) = t2(i) / (ONE + TWO * alpha * cos_wp(w))
        t2(i) = t2(i) / dd
        t2(i) =  - t2(i) * t2(i)
        !write(*, *) i, t2(i)
      end do
  
    else 
    end if
      
      end if
      return
    end subroutine Transform_2nd_derivative_spectral_1d
  !===============================================================================
  !===============================================================================
    subroutine Solve_poisson_000(rhs)
      use decomp_2d_fft!, only : decomp_2d_fft_3d
      implicit none
      real(wp), dimension(:,:,:), intent(INOUT) :: rhs
  !_______________________________________________________________________________
  ! initialize fft if not
  !_______________________________________________________________________________
      if (.not. is_fft_initialised) then
        call decomp_2d_fft_init
        is_fft_initialised = .true.
      end if
  !_______________________________________________________________________________
  ! compute r2c transform, forward FFT
  !_______________________________________________________________________________
      call decomp_2d_fft_3d(rhs,cw1)
  !_______________________________________________________________________________
  ! fft normalisation
  !_______________________________________________________________________________
      cw1 = cw1 / real(nx, kind=wp) / &
                  real(ny, kind=wp) / &
                  real(nz, kind=wp)
  !_______________________________________________________________________________
  ! Fourier domain calculation
  !_______________________________________________________________________________
      cw1(:,:,:) = cw1(:,:,:) / t2xyz(:, :, :) 
  !_______________________________________________________________________________
  ! compute c2r transform, inverse FFT
  !_______________________________________________________________________________
      call decomp_2d_fft_3d(cw1,rhs)
  
      return
    end subroutine Solve_poisson_000
  
  
    subroutine Test_poisson_solver
      use parameters_constant_mod
      use geometry_mod
      implicit none
  
      integer(4) :: k, j, i
      real(WP), allocatable :: rhsphi(:,:,:)
      real(WP) :: solution, x
  
      call Initialize_decomp_poisson(domain)
      allocate(rhsphi(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))); rhsphi = ZERO
      do k = xstart(3),xend(3)
        do j = xstart(2),xend(2)
          do i = xstart(1),xend(1)
            x = domain%h(1)*(real(i, WP)-HALF)
            rhsphi(i, j, k) = dcos(x)*dsin(x)
            ! - dsin( domain%h(3)*(real(k, WP)-HALF) ) &
            !                  - dsin( domain%h(2)*(real(j, WP)-HALF) ) &
            !                  - dsin( domain%h(1)*(real(i, WP)-HALF) )
          end do
        end do
      end do
      
      call Solve_poisson(rhsphi)
  
      do k = xstart(3),xend(3)
        do j = xstart(2),xend(2)
          do i = xstart(1),xend(1)
            x = domain%h(1)*(real(i, WP)-HALF)
            solution = -TWO * dsin( TWO * x )
            !dsin( domain%h(3)*(real(k, WP)-HALF) ) + &
            !           dsin( domain%h(2)*(real(j, WP)-HALF) ) + &
            !           dsin( domain%h(1)*(real(i, WP)-HALF) )
  
            write(*, '(3I5.1, 3ES13.5)') k, j, i, solution, rhsphi(i,j,k), solution / rhsphi(i,j,k)
          end do
        end do
      end do
      deallocate(rhsphi)
  
    end subroutine
  
  end module poisson_mod
  
