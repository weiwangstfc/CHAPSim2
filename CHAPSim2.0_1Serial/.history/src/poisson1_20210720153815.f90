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
! store sine/cosine unit
!_______________________________________________________________________________
  real(mytype), save, allocatable, dimension(:) :: az, bz
  real(mytype), save, allocatable, dimension(:) :: ay, by
  real(mytype), save, allocatable, dimension(:) :: ax, bx
!_______________________________________________________________________________
! FFT library only needs to be initialised once
!_______________________________________________________________________________
  logical, save :: fft_initialised = .false.
!_______________________________________________________________________________
! decomposition object for physical space and spectral space
!_______________________________________________________________________________
  TYPE(DECOMP_INFO), save :: ph
  TYPE(DECOMP_INFO), save :: sp
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
  integer(4), dimension(3), save :: fft_st, fft_en, fft_sz
!_______________________________________________________________________________
! work arrays, 
! naming convention: cw (complex); rw (real); 
!                    b =     ; c = 
!                    1 = X-pencil; 2 = Y-pencil; 3 = Z-pencil
!_______________________________________________________________________________
  real(wp),    allocatable, dimension(:, :, :) :: rhsx
  complex(wp), allocatable, dimension(:, :, :) :: cw1
  complex(wp), allocatable, dimension(:, :, :) :: cwx
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
  private :: Calculate_sine_cosine_unit

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
    logical :: is_half(3)
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
      Solve_poisson => Solve_poisson_100
    else if (bcx == 0 .and. bcy == 0 .and. bcz == 1) then
      Solve_poisson => Solve_poisson_000
    else if (bcx == 0 .and. bcy == 1 .and. bcz == 0) then
      Solve_poisson => Solve_poisson_000
    else if (bcx == 1 .and. bcy == 1) then   ! 110 & 111
      Solve_poisson => Solve_poisson_000
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

    !if (bcx == 1) nx = nx + 1 ! to check ?
    !if (bcy == 1) ny = ny + 1
    !if (bcz == 1) nz = nz + 1
!_______________________________________________________________________________
! initialise 2d-decomp library
! For FFT, by calling decomp_2d_fft_init(PHYSICAL_IN_Z, nx, ny, nz)
!      the physical-space data is stored in Z pencil.
!      the spectral-space data is stored in X pencil after FFT.
! For FFT(real 2 complex)
!      For a 3D real input set of size (nx * ny * nz), the 
!      the complex uotput can be held in an array of size (nx, ny, nz / 2 + 1)
!      for a x-pencil stored spectral data.
!_______________________________________________________________________________
    call decomp_info_init(nx, ny, nz,         ph)
    call decomp_info_init(nx, ny, nz / 2 + 1, sp)

    if (.not. fft_initialised) then
      call decomp_2d_fft_init(PHYSICAL_IN_Z, nx, ny, nz)
      fft_initialised = .true.
    end if
    call decomp_2d_fft_get_size(fft_st, fft_en, fft_sz)

    is_half(1) = .false.
    is_half(2) = .false.
    is_half(3) = .true.

    write(*,*) 'physical domain index, i = ', ph%xst(1), ph%xen(1)
    write(*,*) 'physical domain index, j = ', ph%xst(2), ph%xen(2)
    write(*,*) 'physical domain index, k = ', ph%xst(3), ph%xen(3)

    write(*,*) 'spectral domain index, l = ', sp%xst(1), sp%xen(1)
    write(*,*) 'spectral domain index, m = ', sp%xst(2), sp%xen(2)
    write(*,*) 'spectral domain index, n = ', sp%xst(3), sp%xen(3)

    write(*,*) 'Fourier Domain index,  l = ', fft_st(1), fft_en(1)
    write(*,*) 'Fourier Domain index,  m = ', fft_st(2), fft_en(2)
    write(*,*) 'Fourier Domain index,  n = ', fft_st(3), fft_en(3)

!_______________________________________________________________________________
! preparing sine and cosine factors
!_______________________________________________________________________________
    allocate ( ax( ph%xsz(1) ) ); ax = ZERO
    allocate ( bx( ph%xsz(1) ) ); bx = ZERO

    allocate ( ay( ph%xsz(2) ) ); ay = ZERO
    allocate ( by( ph%xsz(2) ) ); by = ZERO

    allocate ( az( ph%xsz(3) ) ); az = ZERO
    allocate ( bz( ph%xsz(3) ) ); bz = ZERO

    call Calculate_sine_cosine_unit(ax, bx, ph%xsz(1), bcx)
    call Calculate_sine_cosine_unit(ay, by, ph%xsz(2), bcy)
    call Calculate_sine_cosine_unit(az, bz, ph%xsz(3), bcz)

    !_______________________________________________________________________________
! allocate space for wave-space variables
!_______________________________________________________________________________
    allocate ( cw1( sp%xst(1) : sp%xen(1), &
                    sp%xst(2) : sp%xen(2), &
                    sp%xst(3) : sp%xen(3)) )

    allocate ( cwx( sp%xst(1) : sp%xen(1), &
                    sp%xst(2) : sp%xen(2), &
                    sp%xst(3) : sp%xen(3)) )

    if(bcx == 1) then
      allocate ( rhsx(nx, ny, nz) )
    end if
!_______________________________________________________________________________
! prepare the transformation \hat{f"}_l = \hat{f}_l * t2x
!_______________________________________________________________________________
    allocate ( t2x(sp%xst(1) : sp%xen(1)) ) ;  t2x = ZERO
    allocate ( t2y(sp%xst(2) : sp%xen(2)) ) ;  t2y = ZERO
    allocate ( t2z(sp%xst(3) : sp%xen(3)) ) ;  t2z = ZERO
    allocate (t2xyz(  sp%xst(1) : sp%xen(1), &
                      sp%xst(2) : sp%xen(2), &
                      sp%xst(3) : sp%xen(3))) ;  t2xyz = ZERO

    call Transform_2nd_derivative_spectral_1d (bcx, nx, d%h(1), t2x, is_half(1))
    call Transform_2nd_derivative_spectral_1d (bcy, ny, d%h(2), t2y, is_half(2))
    call Transform_2nd_derivative_spectral_1d (bcz, nz, d%h(3), t2z, is_half(3))

    do k = sp%xst(3) , sp%xen(3)
      do j = sp%xst(2) , sp%xen(2)
        do i = sp%xst(1) , sp%xen(1)
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
!> \brief To asign sine and cosine unit
!>
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]    nsz           working array size
!> \param[in]    bc            b.c. flags 
!> \param[out]   afsin         sine unit
!> \param[out]   bfsin         cosine unit
!_______________________________________________________________________________
  subroutine Calculate_sine_cosine_unit(afsin, bfcos, nsz, bc)
    use parameters_constant_mod, only : PI, TWO
    use math_mod
    implicit none
    integer(4), intent(in) :: nsz
    integer(4), intent(in) :: bc
    real(mytype), dimension(:), intent(out) :: afsin
    real(mytype), dimension(:), intent(out) :: bfcos
  
    integer :: i
  
    if (bc == 0) then
  
        do i = 1, nsz
          afsin(i) = sin_wp( real(i - 1, kind = mytype) * PI / &
                             real(nsz,   kind = mytype) )
          bfcos(i) = cos_wp( real(i - 1, kind = mytype) * PI / &
                             real(nsz,   kind = mytype) )
        end do
  
    else if (bc == 1) then
  
        do i = 1, nsz
          afsin(i) = sin_wp( real(i - 1, kind = mytype) * PI / TWO / &
                             real(nsz,   kind = mytype) )
          bfcos(i) = cos_wp( real(i - 1, kind = mytype) * PI / TWO / &
                             real(nsz,   kind = mytype))
        end do
    else
    end if
  
    return
  end subroutine Calculate_sine_cosine_unit
!===============================================================================
!===============================================================================
  subroutine Transform_2nd_derivative_spectral_1d(bc, nn, dd, t2, is_half)
    use operations!,              only : d2fC2C, d2rC2C
    use input_general_mod!,       only : IBC_PERIODIC
    use parameters_constant_mod!, only : FOUR, TWO, ONE, PI
    use math_mod!,                only : cos_wp
    implicit none

    logical :: is_half
    integer(4), intent(in)    :: bc
    integer(4), intent(in)    :: nn
    real(wp),   intent(in)    :: dd
    real(wp),   intent(inout) :: t2(:)

    real(wp) :: a, b, alpha
    real(wp) :: w, cosw, aunit
    integer(4) :: i, iend
    integer(4) :: method = 2

    if(bc == 0) then
      aunit = TWO * PI / REAL(nn, WP)
      iend = nn / 2 + 1
    else
      aunit = PI / REAL(nn, WP)
      iend = nn
    end if

    if(method == 1) then

      alpha = d1fC2C(3, 1, IBC_PERIODIC)
      a     = d1rC2C(3, 1, IBC_PERIODIC) * TWO
      b     = d1rC2C(3, 2, IBC_PERIODIC) * FOUR

      do i = 1, iend
        w = aunit * REAL(i - 1, WP)
        t2(i) = a * sin_wp(w) + b * HALF * sin_wp(TWO * w)
        t2(i) = t2(i) / (ONE + TWO * alpha * cos_wp(w))
        t2(i) = t2(i) / dd
        t2(i) = - t2(i) * t2(i)
        write(*,*) i, t2(i)
      end do

      if((.not. is_half) .and. (bc == 0)) then
        do i = nn/2 + 2, nn
          ! w = aunit * REAL(i - 1, WP)
          ! t2(i) = a * sin_wp(w) + b * HALF * sin_wp(TWO * w)
          ! t2(i) = t2(i) / (ONE + TWO * alpha * cos_wp(w))
          ! t2(i) = t2(i) / dd
          ! t2(i) = - t2(i) * t2(i)

          t2(i) = t2(nn - i + 2)
          write(*,*) i, t2(i), t2(nn - i + 2)
        end do
      end if

    else if(method == 2) then

      alpha = d2fC2C(3, 1, IBC_PERIODIC)
      a     = d2rC2C(3, 1, IBC_PERIODIC)
      b     = d2rC2C(3, 2, IBC_PERIODIC) * FOUR
  
      do i = 1, nn/2 + 1
        w = aunit * REAL(i - 1, WP)! check, for 0-n/2, pi or 2pi?
        cosw = cos_wp(w)
        t2(i) = b * cosw * cosw + TWO * a * cosw - TWO * a - b
        t2(i) = t2(i) / (ONE + TWO * alpha * cosw) / dd / dd
        write(*,*) i, t2(i)
      end do
      
      if((.not. is_half) .and. (bc == 0)) then
        do i = nn/2 + 2, nn
          t2(i) = t2(nn - i + 2)
          write(*,*) i, t2(i)
        end do
      end if

    else
    end if

    if(bc /= 0) then
      t2(1) =  ZERO
    end if

    return
  end subroutine Transform_2nd_derivative_spectral_1d
!===============================================================================
!===============================================================================
  subroutine Solve_poisson_000(rhs)
    use decomp_2d_fft!, only : decomp_2d_fft_3d
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT) :: rhs

    real(wp), dimension(:,:,:), allocatable :: rhs0
    integer :: nn, i, j,k


    allocate(rhs0(nx, ny, nz))
    rhs0 = rhs

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
    !cw1(:,:,:) = cw1(:,:,:) / t2xyz(:, :, :) 
!_______________________________________________________________________________
! compute c2r transform, inverse FFT
!_______________________________________________________________________________
    call decomp_2d_fft_3d(cw1,rhs)

    ! nn = 0
    ! do k = 1, nz
    !   do j = 1, ny
    !     do i = 1, nx
    !       nn = nn + 1
    !       write(*,'(A, 4I4.1, 2ES13.5)') 'check', k, j, i, nn, rhs0(i, j, k), rhs(i, j,k)
    !     end do
    !   end do
    ! end do

    deallocate(rhs0)

    return
  end subroutine Solve_poisson_000

!===============================================================================
  subroutine Solve_poisson_100(rhs)
    use decomp_2d_fft!, only : decomp_2d_fft_3d
    use parameters_constant_mod
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT) :: rhs
    integer(4) :: i, j, k
    real(wp) :: cwRe1, cwRe2, cwIm1, cwIm2
    real(wp) :: aRe1, bRe1, aRe2, bRe2, &
                aIm1, bIm1, aIm2, bIm2

   integer :: nn
   real(wp), dimension(:,:,:), allocatable :: rhs0
   
   allocate(rhs0(nx, ny, nz))
   rhs0 = rhs
!_______________________________________________________________________________
! input data re-organisation to get a periodic sequence in physical domain
!_______________________________________________________________________________
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx/2
          rhsx(i, j, k) = rhs( 2 * i - 1, j,  k)
        end do
        do i = nx/2 + 1, nx
          rhsx(i, j, k) = rhs( 2 * (nx + 1) - 2 * i, j, k)
        end do
      end do
    end do 
!_______________________________________________________________________________
! compute r2c transform, forward FFT
!_______________________________________________________________________________
    call decomp_2d_fft_3d(rhsx,cwx)
!_______________________________________________________________________________
! fft normalisation
!_______________________________________________________________________________
    cwx = cwx / real(nx, kind=wp) / &
                real(ny, kind=wp) / &
                real(nz, kind=wp)
!_______________________________________________________________________________
! Reconstruct FFT(rhs) from the above FFT(rhsx) in spectral domain 
!_______________________________________________________________________________
    do k = sp%xst(3),sp%xen(3)
      do j = sp%xst(2),sp%xen(2)
        cw1(1, j, k) = cwx(1, j, k)
        do i = sp%xst(1) + 1, sp%xen(1)
          cwRe1 = real ( cwx(i, j, k), wp )
          cwIm1 = aimag( cwx(i, j, k) )
          cwRe2 = real ( cwx(sp%xen(1) - i + 2, j, k), wp )
          cwIm2 = aimag( cwx(sp%xen(1) - i + 2, j, k) )

          bRe1 = cwRe1 * bx(i) * HALF
          aRe1 = cwRe1 * ax(i) * HALF

          bIm1 = cwIm1 * bx(i) * HALF
          aIm1 = cwIm1 * ax(i) * HALF

          bRe2 = cwRe2 * bx(i) * HALF
          aRe2 = cwRe2 * ax(i) * HALF

          bIm2 = cwIm2 * bx(i) * HALF
          aIm2 = cwIm2 * ax(i) * HALF

          cw1(i, j, k) = cmplx( bRe1 + aIm1 + bRe2 - aIm2, &
                               -aRe1 + bIm1 + aRe2 + bIm2, WP)

        end do
      end do
    end do
!_______________________________________________________________________________
! Fourier domain calculation
!_______________________________________________________________________________
    cw1(:,:,:) = cw1(:,:,:) / t2xyz(:, :, :) 
!_______________________________________________________________________________
! Reconstruct the func(FFT(rhsx)) from the above func(FFT(rhs))
!_______________________________________________________________________________
    do k = sp%xst(3),sp%xen(3)
      do j = sp%xst(2),sp%xen(2)
        cwx(1, j, k) = cw1(1, j, k)
        do i = sp%xst(1) + 1, sp%xen(1)
          cwRe1 = real ( cw1(i, j, k), wp )
          cwIm1 = aimag( cw1(i, j, k) )
          cwRe2 = real ( cw1(sp%xen(1) - i + 2, j, k), wp )
          cwIm2 = aimag( cw1(sp%xen(1) - i + 2, j, k) )

          bRe1 = cwRe1 * bx(i)
          aRe1 = cwRe1 * ax(i)

          bIm1 = cwIm1 * bx(i)
          aIm1 = cwIm1 * ax(i)

          bRe2 = cwRe2 * bx(i)
          aRe2 = cwRe2 * ax(i)

          bIm2 = cwIm2 * bx(i)
          aIm2 = cwIm2 * ax(i)

          cwx(i, j, k) = cmplx( bRe1 - aIm1 + aRe2 + bIm2, &
                                aRe1 + bIm1 - bRe2 + aIm2, WP)

        end do
      end do
    end do
!_______________________________________________________________________________
! compute c2r transform, inverse FFT
!_______________________________________________________________________________
    call decomp_2d_fft_3d(cwx,rhsx)
!_______________________________________________________________________________
! reorganize the output physical data structure
!_______________________________________________________________________________
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx/2
          rhs( 2 * i - 1, j, k) = rhsx(i, j, k)
        end do
        do i = 1, nx/2
          rhs( 2 * i,     j, k) = rhsx(nx - i + 1, j, k)
        end do
      end do
    end do

    ! nn = 0
    ! do k = 1, nz
    !   do j = 1, ny
    !     do i = 1, nx
    !       nn = nn + 1
    !       write(*,*) 'check', k, j, i, rhs0(i, 1, 1), rhs(i, 1,1), rhs0(i, 1, 1)-rhs(i, 1,1)
    !     end do
    !   end do
    ! end do

    deallocate(rhs0)
    return
  end subroutine Solve_poisson_100
!===============================================================================
  subroutine Test_poisson_solver
    use parameters_constant_mod
    use geometry_mod
    use input_general_mod
    implicit none

    integer(4) :: k, j, i, nn
    real(WP), allocatable :: rhsphi(:,:,:)
    real(WP) :: solution
    real(WP) :: x, y, z
    
    call Initialize_decomp_poisson(domain)
    allocate(rhsphi(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2),ph%xst(3):ph%xen(3))); rhsphi = ZERO
    do k = ph%xst(3),ph%xen(3)
      do j = ph%xst(2),ph%xen(2)
        do i = ph%xst(1),ph%xen(1)
          x = domain%h(1)*(real(i, WP)-HALF)
          y = domain%h(2)*(real(j, WP)-HALF)
          z = domain%h(3)*(real(k, WP)-HALF)
          rhsphi(i, j, k) =  -TWO * dsin( TWO * z )
                           !* dcos( domain%h(2)*(real(j, WP)-HALF) ) !&
                            
          
          !rhsphi(i, j, k) = - dsin( domain%h(3)*(real(k, WP)-HALF) ) &
          !                  - dsin( domain%h(2)*(real(j, WP)-HALF) ) &
          !                  - dsin( domain%h(1)*(real(i, WP)-HALF) )
        end do
      end do
    end do
    
    call Solve_poisson(rhsphi)

    nn = 0
    do k = ph%xst(3),ph%xen(3)
      do j = ph%xst(2),ph%xen(2)
        do i = ph%xst(1),ph%xen(1)
          x = domain%h(1)*(real(i, WP)-HALF)
          z = domain%h(3)*(real(k, WP)-HALF)
          nn = nn + 1
          solution = dcos( z ) * dsin( z ) !- &
                     !dsin( domain%h(1)*(real(i, WP)-HALF) )**2
                     !FOUR * dcos( TWO * domain%h(2)*(real(j, WP)-HALF) ) + &
                     !FOUR * dcos( TWO * domain%h(1)*(real(i, WP)-HALF) )

          write(*, '(4I5.1, 3ES13.5)') i,j,k, nn, solution, rhsphi(i,j,k), solution / rhsphi(i,j,k)
        end do
      end do
    end do
    deallocate(rhsphi)

  end subroutine

end module poisson_mod
