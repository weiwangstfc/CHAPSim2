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
! Transwrtfmt1ion Matrix from \hat{f} to \hat{f''}
!_______________________________________________________________________________
  real(wp),       allocatable, dimension(:)       :: t2x, t2y, t2z
  real(wp), save, allocatable, dimension(:, :, :) :: t2xyz
!_______________________________________________________________________________
! boundary conditions and index
!_______________________________________________________________________________
  logical,                  save :: is_periodic(3)
  integer, dimension(3), save :: fft_st, fft_en, fft_sz
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

  private :: Calculate_compact_coef_in_spectral
  private :: Calculate_sine_cosine_unit

  public :: Test_poisson_solver

contains

!=============================================================================================================================================
!=============================================================================================================================================
!> \brief To calcuate all rhs of momentum eq.
!>
!> This subroutine is called everytime when calcuting the rhs of momentum eqs.
!>
!---------------------------------------------------------------------------------------------------------------------------------------------
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
    integer :: i, j, k

    call Print_debug_start_msg("Initializing variables for Poisson Solver ...")
!_______________________________________________________________________________
! set up boundary flags for periodic b.c.
!_______________________________________________________________________________
    is_periodic(:) = dm%is_periodic(:)
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
    if      (is_periodic(1) .and. &
             is_periodic(2) .and. &
             is_periodic(3)) then

      Solve_poisson => Solve_poisson_000

    else if ((.not. is_periodic(1)) .and. &
                    is_periodic(2)  .and. &
                    is_periodic(3)) then

      Solve_poisson => Solve_poisson_100

    else if (is_periodic(1) .and. &
             is_periodic(2) .and. &
      (.not. is_periodic(3))) then

      Solve_poisson => Solve_poisson_000

    else if (is_periodic(1)  .and. &
      (.not. is_periodic(2)) .and. &
             is_periodic(3)) then

      Solve_poisson => Solve_poisson_000

    else if ((.not. is_periodic(1)) .and. &
             (.not. is_periodic(2))) then   ! 110 & 111

      Solve_poisson => Solve_poisson_000

    else
      stop 'boundary condition not supported'
    end if
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
    call decomp_info_init(dm%nc(1), dm%nc(2), dm%nc(3),         ph)
    call decomp_info_init(dm%nc(1), dm%nc(2), dm%nc(3) / 2 + 1, sp)

    if (.not. fft_initialised) then
      call decomp_2d_fft_init(PHYSICAL_IN_Z, dm%nc(1), dm%nc(2), dm%nc(3))
      fft_initialised = .true.
    end if
    call decomp_2d_fft_get_size(fft_st, fft_en, fft_sz)

#ifdef DEBUG
    write (OUTPUT_UNIT,*) 'physical domain index, i = ', ph%xst(1), ph%xen(1)
    write (OUTPUT_UNIT,*) 'physical domain index, j = ', ph%xst(2), ph%xen(2)
    write (OUTPUT_UNIT,*) 'physical domain index, k = ', ph%xst(3), ph%xen(3)

    write (OUTPUT_UNIT,*) 'spectral domain index, l = ', sp%xst(1), sp%xen(1)
    write (OUTPUT_UNIT,*) 'spectral domain index, m = ', sp%xst(2), sp%xen(2)
    write (OUTPUT_UNIT,*) 'spectral domain index, n = ', sp%xst(3), sp%xen(3)

    write (OUTPUT_UNIT,*) 'Fourier Domain index,  l = ', fft_st(1), fft_en(1)
    write (OUTPUT_UNIT,*) 'Fourier Domain index,  m = ', fft_st(2), fft_en(2)
    write (OUTPUT_UNIT,*) 'Fourier Domain index,  n = ', fft_st(3), fft_en(3)
#endif
!_______________________________________________________________________________
! preparing sine and cosine factors
!_______________________________________________________________________________
    allocate ( ax( ph%xsz(1) ) ); ax = ZERO
    allocate ( bx( ph%xsz(1) ) ); bx = ZERO

    allocate ( ay( ph%xsz(2) ) ); ay = ZERO
    allocate ( by( ph%xsz(2) ) ); by = ZERO

    allocate ( az( ph%xsz(3) ) ); az = ZERO
    allocate ( bz( ph%xsz(3) ) ); bz = ZERO

    call Calculate_sine_cosine_unit(ax, bx, ph%xsz(1), is_periodic(1))
    call Calculate_sine_cosine_unit(ay, by, ph%xsz(2), is_periodic(2))
    call Calculate_sine_cosine_unit(az, bz, ph%xsz(3), is_periodic(3))
!_______________________________________________________________________________
! allocate space for wave-space variables
!_______________________________________________________________________________
    allocate ( cw1( sp%xst(1) : sp%xen(1), &
                    sp%xst(2) : sp%xen(2), &
                    sp%xst(3) : sp%xen(3)) )

    allocate ( cwx( sp%xst(1) : sp%xen(1), &
                    sp%xst(2) : sp%xen(2), &
                    sp%xst(3) : sp%xen(3)) )

    if((.not. is_periodic(1))) then
      allocate ( rhsx(dm%nc(1), dm%nc(2), dm%nc(3)) )
    end if
!_______________________________________________________________________________
! prepare the transwrtfmt1ion \hat{f"}_l = \hat{f}_l * t2x
!_______________________________________________________________________________
    allocate ( t2x(sp%xst(1) : sp%xen(1)) ) ;  t2x = ZERO
    allocate ( t2y(sp%xst(2) : sp%xen(2)) ) ;  t2y = ZERO
    allocate ( t2z(sp%xst(3) : sp%xen(3)) ) ;  t2z = ZERO
    allocate (t2xyz(  sp%xst(1) : sp%xen(1), &
                      sp%xst(2) : sp%xen(2), &
                      sp%xst(3) : sp%xen(3))) ;  t2xyz = ZERO

    call Calculate_compact_coef_in_spectral (is_periodic(1), dm%h(1), ph%xsz(1), sp%xst(1), sp%xen(1), t2x)
    call Calculate_compact_coef_in_spectral (is_periodic(2), dm%h(2), ph%xsz(2), sp%xst(2), sp%xen(2), t2y)
    call Calculate_compact_coef_in_spectral (is_periodic(3), dm%h(3), ph%xsz(3), sp%xst(3), sp%xen(3), t2z)

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
!=============================================================================================================================================
!=============================================================================================================================================
  subroutine Finalize_decomp_poisson
    use decomp_2d_fft
    implicit none

    call decomp_2d_fft_finalize
    is_fft_initialised = .false.

    deallocate(t2xyz)
    deallocate(cw1)

    return
  end subroutine Finalize_decomp_poisson
!=============================================================================================================================================
!=============================================================================================================================================
!> \brief To asign sine and cosine unit
!>
!>
!---------------------------------------------------------------------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]    nsz           working array size
!> \param[in]    bc            b.c. flags 
!> \param[out]   afsin         sine unit
!> \param[out]   bfsin         cosine unit
!_______________________________________________________________________________
  subroutine Calculate_sine_cosine_unit(afsin, bfcos, nsz, is_peri)
    use parameters_constant_mod, only : PI, TWO
    use math_mod
    implicit none
    integer, intent(in) :: nsz
    logical,    intent(in) :: is_peri
    real(mytype), dimension(:), intent(out) :: afsin
    real(mytype), dimension(:), intent(out) :: bfcos
  
    integer :: i
  
    if (is_peri) then
  
      do i = 1, nsz
        afsin(i) = sin_wp( real(i - 1, kind = mytype) * PI / &
                            real(nsz,   kind = mytype) )
        bfcos(i) = cos_wp( real(i - 1, kind = mytype) * PI / &
                            real(nsz,   kind = mytype) )
      end do
  
    else
  
      do i = 1, nsz
        afsin(i) = sin_wp( real(i - 1, kind = mytype) * PI / TWO / &
                            real(nsz,   kind = mytype) )
        bfcos(i) = cos_wp( real(i - 1, kind = mytype) * PI / TWO / &
                            real(nsz,   kind = mytype))
      end do
    end if
  
    return
  end subroutine Calculate_sine_cosine_unit
!=============================================================================================================================================
!=============================================================================================================================================
  subroutine Calculate_compact_coef_in_spectral(is_peri, dd, nn, istt, iend, t2)
    use udf_type_mod,            only : t_domain
    use operations!,              only : d2fC2C, d2rC2C
    use input_general_mod!,       only : IBC_PERIODIC
    use parameters_constant_mod!, only : FOUR, TWO, ONE, PI
    use math_mod!,                only : cos_wp
    implicit none

    logical,  intent(in) :: is_peri
    real(wp), intent(in) :: dd
    integer,  intent(in) :: nn, istt, iend
    real(wp),   intent(inout) :: t2(:)

    real(wp) :: a, b, alpha
    real(wp) :: w, cosw, aunit
    integer :: i

    if(is_peri) then
      aunit = TWO * PI / REAL(nn, WP)
    else
      aunit = PI / REAL(nn, WP)
    end if

    if(ifft2deri == 1) then

      alpha = d1fC2C(3, 1, IBC_PERIODIC)
      a     = d1rC2C(3, 1, IBC_PERIODIC) * TWO
      b     = d1rC2C(3, 2, IBC_PERIODIC) * FOUR

      do i = istt, iend
        w = aunit * REAL(i - 1, WP)
        t2(i) = a * sin_wp(w) + b * HALF * sin_wp(TWO * w)
        t2(i) = t2(i) / (ONE + TWO * alpha * cos_wp(w))
        t2(i) = t2(i) / dd
        t2(i) = - t2(i) * t2(i)
        write (OUTPUT_UNIT,*) i, t2(i)
      end do

    else if(ifft2deri == 2) then

      alpha = d2fC2C(3, 1, IBC_PERIODIC)
      a     = d2rC2C(3, 1, IBC_PERIODIC)
      b     = d2rC2C(3, 2, IBC_PERIODIC) * FOUR
  
      do i = istt, iend
        w = aunit * REAL(i - 1, WP)! check, for 0-n/2, pi or 2pi?
        cosw = cos_wp(w)
        t2(i) = b * cosw * cosw + TWO * a * cosw - TWO * a - b
        t2(i) = t2(i) / (ONE + TWO * alpha * cosw) / dd / dd
        write (OUTPUT_UNIT,*) i, t2(i)
      end do

    else
    end if

    if(.not. is_peri) then
      t2(1) =  ZERO
    end if

    return
  end subroutine Calculate_compact_coef_in_spectral
!=============================================================================================================================================
!=============================================================================================================================================
  subroutine Solve_poisson_000(rhs)
    use decomp_2d_fft!, only : decomp_2d_fft_3d
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT) :: rhs

#ifdef DEBUG
    integer :: nn, i, j,k
    real(wp), allocatable :: rhs0(:, :, :)

    allocate (rhs0 (ph%xsz(1), ph%xsz(2), ph%xsz(3)))
    rhs0 = rhs
#endif    
!_______________________________________________________________________________
! compute r2c transform, forward FFT
!_______________________________________________________________________________
    call decomp_2d_fft_3d(rhs,cw1)
!_______________________________________________________________________________
! fft normalisation
!_______________________________________________________________________________
    cw1 = cw1 / real(ph%xsz(1), kind=wp) / &
                real(ph%xsz(2), kind=wp) / &
                real(ph%xsz(3), kind=wp)

#ifdef DEBUG
    nn = 0
    do k = sp%xst(3), sp%xen(3)
      do j = sp%xst(2), sp%xen(2)
        do i = sp%xst(1), sp%xen(1)
          write (OUTPUT_UNIT,'(A, 3I5.1, 1ES13.5)') 'After F-FFT', k, j, i, cw1(i,j,k)
        end do
      end do
    end do
#endif 
!_______________________________________________________________________________
! Fourier domain calculation
!_______________________________________________________________________________
    cw1(:,:,:) = cw1(:,:,:) / t2xyz(:, :, :) 
!_______________________________________________________________________________
! compute c2r transform, inverse FFT
!_______________________________________________________________________________
    call decomp_2d_fft_3d(cw1,rhs)

#ifdef DEBUG
    do k = ph%xst(3), ph%xen(3)
      do j = ph%xst(2), ph%xen(2)
        do i = ph%xst(1), ph%xen(1)
          write (OUTPUT_UNIT,'(A, 3I4.1, 2ES13.5)') 'After B-FFT', k, j, i, rhs0(i, j, k), rhs(i, j,k)
        end do
      end do
    end do

    deallocate(rhs0)
#endif
    
    return
  end subroutine Solve_poisson_000

!=============================================================================================================================================
  subroutine Solve_poisson_100(rhs)
    use decomp_2d_fft!, only : decomp_2d_fft_3d
    use parameters_constant_mod
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT) :: rhs
    integer :: i, j, k
    real(wp) :: cwRe1, cwRe2, cwIm1, cwIm2
    real(wp) :: aRe1, bRe1, aRe2, bRe2, &
                aIm1, bIm1, aIm2, bIm2

#ifdef DEBUG
    real(wp), allocatable :: rhs0(:, :, :)

    allocate (rhs0 (ph%xsz(1), ph%xsz(2), ph%xsz(3)))
    rhs0 = rhs
#endif  
!_______________________________________________________________________________
! input data re-organisation to get a periodic sequence in physical domain
!_______________________________________________________________________________
    do k = ph%xst(3), ph%xen(3)
      do j = ph%xst(2), ph%xen(2)
        do i = ph%xst(1), ph%xsz(1)/2
          rhsx(i, j, k) = rhs( 2 * i - 1, j,  k)
        end do
        do i = ph%xsz(1)/2 + 1, ph%xen(1)
          rhsx(i, j, k) = rhs( 2 * (ph%xsz(1) + 1) - 2 * i, j, k)
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
    cwx = cwx / real(ph%xsz(1), kind=wp) / &
                real(ph%xsz(2), kind=wp) / &
                real(ph%xsz(3), kind=wp)
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
    do k = ph%xst(3), ph%xen(3)
      do j = ph%xst(2), ph%xen(2)
        do i = ph%xst(1), ph%xsz(1)/2
          rhs( 2 * i - 1, j, k) = rhsx(i, j, k)
        end do
        do i = ph%xst(1), ph%xsz(1)/2
          rhs( 2 * i,     j, k) = rhsx(ph%xsz(1) - i + 1, j, k)
        end do
      end do
    end do

#ifdef DEBUG
    do k = ph%xst(3), ph%xen(3)
      do j = ph%xst(2), ph%xen(2)
        do i = ph%xst(1), ph%xen(1)
          write (OUTPUT_UNIT,'(A, 3I4.1, 2ES13.5)') 'After B-FFT', k, j, i, rhs0(i, j, k), rhs(i, j,k)
        end do
      end do
    end do

    deallocate(rhs0)
#endif
    return
  end subroutine Solve_poisson_100
!=============================================================================================================================================
  subroutine Test_poisson_solver
    use parameters_constant_mod
    use geometry_mod
    use input_general_mod
    implicit none

    integer :: k, j, i, nn
    real(WP), allocatable :: rhsphi(:,:,:)
    real(WP) :: solution
    real(WP) :: x, y, z
    
    !call Initialize_decomp_poisson(domain)
    allocate(rhsphi(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2),ph%xst(3):ph%xen(3))); rhsphi = ZERO
    do i = ph%xst(1),ph%xen(1)
      do j = ph%xst(2),ph%xen(2)
        do k = ph%xst(3),ph%xen(3)
          x = domain%h(1)*(real(i, WP)-HALF)
          y = domain%h(2)*(real(j, WP)-HALF)
          z = domain%h(3)*(real(k, WP)-HALF)
          rhsphi(i, j, k) = -sin(half * x)/FOUR
          !  -TWO * dsin( TWO * z )
                           !* dcos( domain%h(2)*(real(j, WP)-HALF) ) !&
                            
          
          !rhsphi(i, j, k) = - dsin( domain%h(3)*(real(k, WP)-HALF) ) &
          !                  - dsin( domain%h(2)*(real(j, WP)-HALF) ) &
          !                  - dsin( domain%h(1)*(real(i, WP)-HALF) )
        end do
      end do
    end do
    
    call Solve_poisson(rhsphi)

    nn = 0
    do i = ph%xst(1),ph%xen(1)
      do j = ph%xst(2),ph%xen(2)
        do k = ph%xst(3),ph%xen(3)
        !k = 1
          x = domain%h(1)*(real(i, WP)-HALF)
          y = domain%h(2)*(real(j, WP)-HALF)
          z = domain%h(3)*(real(k, WP)-HALF)
          nn = nn + 1
          solution = dsin(HALF * x)
          !dcos( z ) * dsin( z ) !- &
                     !dsin( domain%h(1)*(real(i, WP)-HALF) )**2
                     !FOUR * dcos( TWO * domain%h(2)*(real(j, WP)-HALF) ) + &
                     !FOUR * dcos( TWO * domain%h(1)*(real(i, WP)-HALF) )

          write (OUTPUT_UNIT, '(4I5.1, 3ES13.5)') i,j,k, nn, solution, rhsphi(i,j,k), solution / rhsphi(i,j,k)
        end do
      end do
    end do
    deallocate(rhsphi)
    stop
  end subroutine

end module poisson_mod
