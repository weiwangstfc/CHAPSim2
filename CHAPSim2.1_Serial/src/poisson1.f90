
!===============================================================================
!===============================================================================
!Property Used :
! (1) if the input data in[i] are purely real numbers, in which case the DFT 
! output satisfies the “Hermitian” redundancy: 
! out[i] is the conjugate of out[n-i]
! (2)  DFT of real number :
! the 0th (the “DC”) and n/2-th (the “Nyquist” frequency, when n is even) 
! elements of the cmplx output are purely real.
! (3) For 2D and 3D FFTs, the FFT length along the innermost dimension is used 
! to compute the  value. This is because the FFT along the innermost dimension 
! is computed first and is logically a real-to-hermitian transform. 
! The FFTs along other dimensions are computed next, and they are simply 
! ‘complex-to-complex’ transforms. 
!===============================================================================
!===============================================================================
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
  CHARACTER(LEN=*), PARAMETER :: complexfmt = '(ES13.5,SP,ES13.5,"i")'
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
  cmplx(wp),       allocatable, dimension(:)       :: t2x, t2y, t2z
  cmplx(wp), save, allocatable, dimension(:, :, :) :: t2xyz
!_______________________________________________________________________________
! boundary conditions and index
!_______________________________________________________________________________
  logical,                  save :: is_periodic(3)
  integer(4), dimension(3), save :: fft_st, fft_en, fft_sz
  integer(4), dimension(3), save :: nw
!_______________________________________________________________________________
! work arrays, 
! naming convention: cw (cmplx); rw (real); 
!                    b =     ; c = 
!                    1 = X-pencil; 2 = Y-pencil; 3 = Z-pencil
!_______________________________________________________________________________
  real(wp),    allocatable, dimension(:, :, :) :: rw_xpen, rw_ypen, rw_zpen
  real(wp),    allocatable, dimension(:, :, :) :: rw_recons_xpen, rw_recons_ypen, rw_recons_zpen
  cmplx(wp), allocatable, dimension(:, :, :) :: cw_xpen, cw_ypen, cw_zpen
  cmplx(wp), allocatable, dimension(:, :, :) :: cw_recons_xpen, cw_recons_ypen, cw_recons_zpen
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
  subroutine Calculate_sine_cosine_unit(afsin, bfcos, nsz, is_peri)
    use parameters_constant_mod, only : PI, TWO
    use math_mod
    implicit none
    integer(4), intent(in) :: nsz
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

!===============================================================================
!===============================================================================
  subroutine Calculate_compact_coef_in_spectral(is_peri, dd, nc, t2)
    use udf_type_mod,            only : t_domain
    use operations!,              only : d2fC2C, d2rC2C
    use input_general_mod!,       only : IBC_PERIODIC
    use parameters_constant_mod!, only : FOUR, TWO, ONE, PI
    use math_mod!,                only : cos_wp
    implicit none

    logical,  intent(in) :: is_peri
    real(wp), intent(in) :: dd ! dx
    integer,  intent(in) :: nc ! 
    cmplx(wp),   intent(inout) :: t2(:)

    real(wp) :: a, b, alpha
    real(wp) :: w, cosw, aunit
    real(wp) :: tr, ti
    cmplx(wp) :: tc, ic
    integer(4) :: l

    if(is_peri) then

      aunit = TWO * PI / REAL(nc, WP)
      alpha = d1fC2C(3, 1, IBC_PERIODIC)
      a     = d1rC2C(3, 1, IBC_PERIODIC) * TWO
      b     = d1rC2C(3, 2, IBC_PERIODIC) * FOUR
      !write(*,*) 'alpha, a, b', alpha, a, b
      do l = 1, nc / 2 + 1
        w = aunit * REAL(l - 1, WP)
        tr = a * sin_wp(w) + b * HALF * sin_wp(TWO * w)
        tr = tr / (ONE + TWO * alpha * cos_wp(w))
        tr = tr / dd
        t2(l) = cmplx(- tr * tr, ZERO)
      end do

       do l = nc / 2 + 2, nc
         t2(l) = t2(nc - l + 2)
       end do

    else
!_______________________________________________________________________________
!     non-periodic pressure bc. is Neumann B.C.
!     for pressure stored at cell centre
!     it could be represented by symmetric b.c. ? check!
!_______________________________________________________________________________     
      ic = cmplx(ZERO, ONE)

      aunit = PI / REAL(nc, WP)
!_______________________________________________________________________________
!     bc = 1
!_______________________________________________________________________________
      l = 1
      w = aunit * real(l - 1, WP)

      alpha = d1fC2C(1, 1, IBC_SYMMETRIC)
      a     = d1rC2C(1, 1, IBC_SYMMETRIC) * TWO
      b     = d1rC2C(1, 2, IBC_SYMMETRIC) * FOUR

      tc = (TWO * a - b) / FOUR / dd * exp(ic * w) + &
          b / FOUR / dd * exp(ic * TWO * w) - &
          a / TWO / dd

      tc = tc / (-alpha + ONE + alpha * exp(ic * w))

      t2(l) = tc * tc
!_______________________________________________________________________________
!     bc = 2
!_______________________________________________________________________________
      l = 2
      w = aunit * real(l - 1, WP)

      alpha = d1fC2C(2, 1, IBC_SYMMETRIC)
      a     = d1rC2C(2, 1, IBC_SYMMETRIC) * TWO
      b     = d1rC2C(2, 2, IBC_SYMMETRIC) * FOUR

      tc = b / FOUR / dd * exp(ic * TWO * w) + &
          a / TWO / dd * exp(ic * w) - &
          (TWO * a + b) / FOUR / dd * exp(-ic * w)

      tc = tc / (ONE + TWO * alpha * dcos(w))
      t2(l) = tc * tc
!_______________________________________________________________________________
!     bc = bulk
!_______________________________________________________________________________
      alpha = d1fC2C(3, 1, IBC_SYMMETRIC)
      a     = d1rC2C(3, 1, IBC_SYMMETRIC) * TWO
      b     = d1rC2C(3, 2, IBC_SYMMETRIC) * FOUR

      do l = 3, nc - 2
        w = aunit * REAL(l - 1, WP)
        tr = a * sin_wp(w) + b * HALF * sin_wp(TWO * w)
        tr = tr / (ONE + TWO * alpha * cos_wp(w))
        tr = tr / dd
        t2(l) = cmplx(- tr * tr, ZERO)
      end do
!_______________________________________________________________________________
!     bc = n - 1
!_______________________________________________________________________________
      l = nc - 1
      w = aunit * real(l - 1, WP)
      alpha = d1fC2C(4, 1, IBC_SYMMETRIC)
      a     = d1rC2C(4, 1, IBC_SYMMETRIC) * TWO
      b     = d1rC2C(4, 2, IBC_SYMMETRIC) * FOUR

      tc = -b / FOUR / dd * exp(-ic * TWO * w) - &
          a / TWO / dd * exp(-ic * w) + &
          (TWO * a + b) / FOUR / dd * exp(ic * w)

      tc = tc / (ONE + TWO * alpha * dcos(w))
      t2(l) = tc * tc
!_______________________________________________________________________________
!     bc = n
!_______________________________________________________________________________
      l = nc
      w = aunit * real(l - 1, WP)

      alpha = d1fC2C(5, 1, IBC_SYMMETRIC)
      a     = d1rC2C(5, 1, IBC_SYMMETRIC) * TWO
      b     = d1rC2C(5, 2, IBC_SYMMETRIC) * FOUR

      tc = -(TWO * a - b) / FOUR / dd * exp(-ic * w) - &
          b / FOUR / dd * exp(-ic * TWO * w) + &
          a / TWO / dd

      tc = tc / (-alpha + ONE + alpha * exp(-ic * w))

      t2(l) = tc * tc
    end if

#ifdef DEBUGFFT
    do l = 1, nc
      write(*, *) 'modified wavenumber = ', l, t2(l)
    end do
#endif 

    return
  end subroutine Calculate_compact_coef_in_spectral
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
    use udf_type_mod!,            only : t_domain
    use parameters_constant_mod!, only : ZERO, MAXP, TRUNCERR
    implicit none
    type(t_domain), intent(in) :: d
    integer(4) :: i, j, k
    integer(4) :: it1s, it2s, it3s, it1e, it2e, it3e

    call Print_debug_start_msg("Initializing variables for Poisson Solver ...")
!_______________________________________________________________________________
! set up boundary flags for periodic b.c.
!_______________________________________________________________________________
    is_periodic(:) = d%is_periodic(:)
    nw(:) = d%nc(:) ! only for periodic b.c.
!_______________________________________________________________________________
! boundary conditions, periodic or not
!   x    y     z
!   0    0     0 (x, y, z periodic)
!   1    0     0 (y, z periodic)
!   0    1     0 (x, z periodic)
!   0    0     1 (x, y periodic)
!   1    1     0 (z periodic)
!   1    0     1 (y periodic)
!   0    1     1 (x periodic)
!   1    1     1 (none periodic)
!_______________________________________________________________________________
    if      (is_periodic(1) .and. &
             is_periodic(2) .and. &
             is_periodic(3)) then
      Solve_poisson => Solve_poisson_000
    else if ((.not. is_periodic(1)) .and. &
                    is_periodic(2)  .and. &
                    is_periodic(3)) then
      Solve_poisson => Solve_poisson_100
    else if (       is_periodic(1)  .and. &
             (.not. is_periodic(2)) .and. &
                    is_periodic(3)) then
      Solve_poisson => Solve_poisson_010
    else if (       is_periodic(1)  .and. &
                    is_periodic(2)  .and. &
             (.not. is_periodic(3))) then
      Solve_poisson => Solve_poisson_001
    else if ((.not. is_periodic(1)) .and. &
             (.not. is_periodic(2)) .and. &
                    is_periodic(3))  then
      Solve_poisson => Solve_poisson_110
    else if ((.not. is_periodic(1)) .and. &
                    is_periodic(2)  .and. &
             (.not. is_periodic(3))) then
      Solve_poisson => Solve_poisson_101
    else if (       is_periodic(1)  .and. &
             (.not. is_periodic(2)) .and. &
             (.not. is_periodic(3))) then
      Solve_poisson => Solve_poisson_011
    else if ( (.not. is_periodic(1)) .and. &
              (.not. is_periodic(2)) .and. &
              (.not. is_periodic(3))) then
      Solve_poisson => Solve_poisson_111
    else
      stop 'boundary condition not supported'
    end if
!_______________________________________________________________________________
! preparing sine and cosine factors
!_______________________________________________________________________________
    allocate ( ax( nw(1) ) ); ax = ZERO
    allocate ( bx( nw(1) ) ); bx = ZERO

    allocate ( ay( nw(2) ) ); ay = ZERO
    allocate ( by( nw(2) ) ); by = ZERO

    allocate ( az( nw(3) ) ); az = ZERO
    allocate ( bz( nw(3) ) ); bz = ZERO

    call Calculate_sine_cosine_unit(ax, bx, nw(1), is_periodic(1))
    call Calculate_sine_cosine_unit(ay, by, nw(2), is_periodic(2))
    call Calculate_sine_cosine_unit(az, bz, nw(3), is_periodic(3))
!_______________________________________________________________________________
! initialise 2d-decomp library
! For FFT, by calling decomp_2d_fft_init(PHYSICAL_IN_Z, nx, ny, nz)
!      the physical-space data is stored in Z pencil.
!      the spectral-space data is stored in X pencil after FFT.
! For FFT(real 2 cmplx)
!      For a 3D real input set of size (nx * ny * nz), the 
!      the cmplx uotput can be held in an array of size (nx, ny, nz / 2 + 1)
!      for a x-pencil stored spectral data.
!_______________________________________________________________________________
!_______________________________________________________________________________
! initialise 2d-decomp library using advanced 2d decomposition API
! Derived data type DECOMP_INFO: ph, physical domain of input rhs 
!                                size: nclx, ncly, nclz
! Derived data type DECOMP_INFO: sp, spectral domain of output fft(rhs)
!                                size: nclx, ncly, nclz / 2 + 1
! Note: default FFT API: decomp_2d_fft_init (PHYSICAL_IN_X)
!                        physical-space data stored in x pencil format
!                        spectral-space data stored in z pencil format
!_______________________________________________________________________________
    call decomp_2d_init(nw(1), nw(2), nw(3), 1, 1) ! serial
    call decomp_info_init(nw(1), nw(2), nw(3),         ph)
    call decomp_info_init(nw(1), nw(2), nw(3) / 2 + 1, sp)

#ifdef DEBUGFFT
    !if(myid == 0) write(*, *) 'Physical Domain (Poisson), 2D decompisiton : '
    write(*, *) 'In x-pencil size in x-dir = ', ph%xsz(1), ' x id from ', ph%xst(1), ph%xen(1)
    write(*, *) 'In x-pencil size in y-dir = ', ph%xsz(2), ' y id from ', ph%xst(2), ph%xen(2)
    write(*, *) 'In x-pencil size in z-dir = ', ph%xsz(3), ' z id from ', ph%xst(3), ph%xen(3)
    write(*, *) 'In y-pencil size in x-dir = ', ph%ysz(1), ' x id from ', ph%yst(1), ph%yen(1)
    write(*, *) 'In y-pencil size in y-dir = ', ph%ysz(2), ' y id from ', ph%yst(2), ph%yen(2)
    write(*, *) 'In y-pencil size in z-dir = ', ph%ysz(3), ' z id from ', ph%yst(3), ph%yen(3)
    write(*, *) 'In z-pencil size in x-dir = ', ph%zsz(1), ' x id from ', ph%zst(1), ph%zen(1)
    write(*, *) 'In z-pencil size in y-dir = ', ph%zsz(2), ' y id from ', ph%zst(2), ph%zen(2)
    write(*, *) 'In z-pencil size in z-dir = ', ph%zsz(3), ' z id from ', ph%zst(3), ph%zen(3)
#endif

#ifdef DEBUGFFT
    !if(myid == 0) write(*, *) 'Spectral domain (Poisson), 2D decompisiton : '
    write(*, *) 'In x-pencil size in x-dir = ', sp%xsz(1), ' x id from ', sp%xst(1), sp%xen(1)
    write(*, *) 'In x-pencil size in y-dir = ', sp%xsz(2), ' y id from ', sp%xst(2), sp%xen(2)
    write(*, *) 'In x-pencil size in z-dir = ', sp%xsz(3), ' z id from ', sp%xst(3), sp%xen(3)
    write(*, *) 'In y-pencil size in x-dir = ', sp%ysz(1), ' x id from ', sp%yst(1), sp%yen(1)
    write(*, *) 'In y-pencil size in y-dir = ', sp%ysz(2), ' y id from ', sp%yst(2), sp%yen(2)
    write(*, *) 'In y-pencil size in z-dir = ', sp%ysz(3), ' z id from ', sp%yst(3), sp%yen(3)
    write(*, *) 'In z-pencil size in x-dir = ', sp%zsz(1), ' x id from ', sp%zst(1), sp%zen(1)
    write(*, *) 'In z-pencil size in y-dir = ', sp%zsz(2), ' y id from ', sp%zst(2), sp%zen(2)
    write(*, *) 'In z-pencil size in z-dir = ', sp%zsz(3), ' z id from ', sp%zst(3), sp%zen(3)
#endif

!_______________________________________________________________________________
! allocate space for wave-space variables
!_______________________________________________________________________________
    if      (is_periodic(1) .and. &
             is_periodic(2) .and. &
             is_periodic(3)) then ! 000
    
      allocate (t2xyz( sp%xst(1) : sp%xen(1), &
                       sp%xst(2) : sp%xen(2), &
                       sp%xst(3) : sp%xen(3))) ;  t2xyz = cmplx(ZERO, ZERO)

      allocate ( cw_xpen( sp%xst(1) : sp%xen(1), &
                          sp%xst(2) : sp%xen(2), &
                          sp%xst(3) : sp%xen(3)) )
      
    else if ((.not. is_periodic(1)) .and. &
                    is_periodic(2)  .and. &
                    is_periodic(3)) then ! 100

      allocate (t2xyz( sp%xst(1) : sp%xen(1), &
                       sp%xst(2) : sp%xen(2), &
                       sp%xst(3) : sp%xen(3))) ;  t2xyz = cmplx(ZERO, ZERO)

      allocate ( rw_xpen( ph%xst(1) : ph%xen(1), &
                          ph%xst(2) : ph%xen(2), &
                          ph%xst(3) : ph%xen(3)) )

      allocate ( rw_ypen( ph%yst(1) : ph%yen(1), &
                          ph%yst(2) : ph%yen(2), &
                          ph%yst(3) : ph%yen(3)) )

      allocate ( rw_zpen( ph%zst(1) : ph%zen(1), &
                          ph%zst(2) : ph%zen(2), &
                          ph%zst(3) : ph%zen(3)) )

      allocate ( rw_recons_xpen( ph%xst(1) : ph%xen(1), &
                                 ph%xst(2) : ph%xen(2), &
                                 ph%xst(3) : ph%xen(3)) )

      allocate ( cw_xpen       ( sp%xst(1) : sp%xen(1), &
                                 sp%xst(2) : sp%xen(2), &
                                 sp%xst(3) : sp%xen(3)) )

      allocate ( cw_recons_xpen( sp%xst(1) : sp%xen(1), &
                                 sp%xst(2) : sp%xen(2), &
                                 sp%xst(3) : sp%xen(3)) )
      
    else if (       is_periodic(1)  .and. &
             (.not. is_periodic(2)) .and. &
                    is_periodic(3)) then ! 010

      allocate (t2xyz( sp%yst(1) : sp%yen(1), &
                       sp%yst(2) : sp%yen(2), &
                       sp%yst(3) : sp%yen(3))) ;  t2xyz = cmplx(ZERO, ZERO)

      allocate ( rw_ypen( ph%yst(1) : ph%yen(1), &
                          ph%yst(2) : ph%yen(2), &
                          ph%yst(3) : ph%yen(3)) )

      allocate ( rw_zpen( ph%zst(1) : ph%zen(1), &
                          ph%zst(2) : ph%zen(2), &
                          ph%zst(3) : ph%zen(3)) )

      allocate ( rw_recons_ypen( ph%yst(1) : ph%yen(1), &
                                 ph%yst(2) : ph%yen(2), &
                                 ph%yst(3) : ph%yen(3)) )

      allocate ( cw_xpen       ( sp%xst(1) : sp%xen(1), &
                                 sp%xst(2) : sp%xen(2), &
                                 sp%xst(3) : sp%xen(3)) )

      allocate ( cw_ypen       ( sp%yst(1) : sp%yen(1), &
                                 sp%yst(2) : sp%yen(2), &
                                 sp%yst(3) : sp%yen(3)) )

      allocate ( cw_recons_ypen( sp%yst(1) : sp%yen(1), &
                                 sp%yst(2) : sp%yen(2), &
                                 sp%yst(3) : sp%yen(3)) )
      
    else if (       is_periodic(1)  .and. &
                    is_periodic(2)  .and. &
             (.not. is_periodic(3))) then ! 001
      
    else if ((.not. is_periodic(1)) .and. &
             (.not. is_periodic(2)) .and. &
                    is_periodic(3))  then ! 110
      
    else if ((.not. is_periodic(1)) .and. &
                    is_periodic(2)  .and. &
             (.not. is_periodic(3))) then ! 101
      
    else if (       is_periodic(1)  .and. &
             (.not. is_periodic(2)) .and. &
             (.not. is_periodic(3))) then ! 011
      
    else if ( (.not. is_periodic(1)) .and. &
              (.not. is_periodic(2)) .and. &
              (.not. is_periodic(3))) then ! 111
      
    else
      stop 'boundary condition not supported'
    end if
!_______________________________________________________________________________
! prepare the transformation \hat{f"}_l = \hat{f}_l * t2x
! the operation of spetral data is in x-pencil (from PHYSICAL_IN_Z)
!_______________________________________________________________________________
    allocate ( t2x( nw(1) ) ) ;  t2x = cmplx(ZERO, ZERO)
    allocate ( t2y( nw(2) ) ) ;  t2y = cmplx(ZERO, ZERO)
    allocate ( t2z( nw(3) ) ) ;  t2z = cmplx(ZERO, ZERO)

    call Calculate_compact_coef_in_spectral (is_periodic(1), d%h(1), nw(1), t2x)
    call Calculate_compact_coef_in_spectral (is_periodic(2), d%h(2), nw(2), t2y)
    call Calculate_compact_coef_in_spectral (is_periodic(3), d%h(3), nw(3), t2z)

    if (       is_periodic(1)  .and. &
        (.not. is_periodic(2)) .and. &
               is_periodic(3)) then ! 010
      it1s = sp%yst(1)
      it1e = sp%yen(1)
      it2s = sp%yst(2) 
      it2e = sp%yen(2)
      it3s = sp%yst(3) 
      it3e = sp%yen(3)
    else
      it1s = sp%xst(1)
      it1e = sp%xen(1)
      it2s = sp%xst(2) 
      it2e = sp%xen(2)
      it3s = sp%xst(3) 
      it3e = sp%xen(3)
    end if

    do k = it3s, it3e
      do j = it2s, it2e
        do i = it1s, it1e
          t2xyz(i, j, k) = t2x(i) + t2y(j) + t2z(k) ! (spectral space)
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

    return
  end subroutine Finalize_decomp_poisson
!===============================================================================
!===============================================================================
  subroutine Solve_poisson_000(rhs)
    use decomp_2d_fft!, only : decomp_2d_fft_3d
    use parameters_constant_mod
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT) :: rhs

    integer :: i, j, k
#ifdef DEBUGFFT
    integer :: nn
    real(wp), allocatable :: rhs0(:, :, :)

    allocate (rhs0 (ph%xsz(1), ph%xsz(2), ph%xsz(3)))
    rhs0 = rhs

! cw1 : x pencil format
    do k = sp%xst(3), sp%xen(3)
      do j = sp%xst(2), sp%xen(2)
        do i = sp%xst(1), sp%xen(1)
          write(*, '(A, 3I4.1, 1ES13.5)') 'input', k, j, i, rhs(i,j,k)
        end do
      end do
    end do
#endif 
!_______________________________________________________________________________
! initialize FFT
!_______________________________________________________________________________
    if (.not. fft_initialised) then
      call decomp_2d_fft_init(PHYSICAL_IN_Z)!, nw(1), nw(2), nw(3))
      fft_initialised = .true.
    end if
!_______________________________________________________________________________
! compute r2c transform, forward FFT
! rhs : default z pencil format
! cw1 : x pencil format
!_______________________________________________________________________________
    call decomp_2d_fft_3d(rhs, cw_xpen)
!_______________________________________________________________________________
! fft normalisation
!_______________________________________________________________________________
    cw_xpen = cw_xpen / real(nw(1), kind=wp) / &
                        real(nw(2), kind=wp) / &
                        real(nw(3), kind=wp)

#ifdef DEBUGFFT
! cw1 : x pencil format
    do k = sp%xst(3), sp%xen(3)
      do j = sp%xst(2), sp%xen(2)
        do i = sp%xst(1), sp%xen(1)
          write(*, '(A, 3I5.1, 4ES13.5)') 'After F-FFT', k, j, i, & 
            cw_xpen(i,j,k), t2xyz(i,j,k)
        end do
      end do
    end do
#endif 
!_______________________________________________________________________________
! Fourier domain calculation
!_______________________________________________________________________________
  do k = sp%xst(3), sp%xen(3)
    do j = sp%xst(2), sp%xen(2)
      do i = sp%xst(1), sp%xen(1)
        if ( ( dabs( real(t2xyz(i, j, k), WP)) < TRUNCERR ) .AND. &
             ( dabs(aimag(t2xyz(i, j, k))) < TRUNCERR )) then
          cw_xpen(i, j, k) = cmplx(ZERO, ZERO)
        else 
          cw_xpen(i, j, k) = cw_xpen(i, j, k) / t2xyz(i, j, k)
        end if
      end do
    end do
  end do
!_______________________________________________________________________________
#ifdef DEBUGFFT
! cw1 : x pencil format
    do k = sp%xst(3), sp%xen(3)
      do j = sp%xst(2), sp%xen(2)
        do i = sp%xst(1), sp%xen(1)
          write(*, '(A, 3I4.1, 2ES13.5)') 'D/F in spectral', k, j, i, cw_xpen(i,j,k)
        end do
      end do
    end do
#endif 
!_______________________________________________________________________________
! compute c2r transform, inverse FFT
! rhs : z pencil format
! cw1 : x pencil format
!_______________________________________________________________________________
    call decomp_2d_fft_3d(cw_xpen, rhs)

#ifdef DEBUGFFT
    do k = ph%zst(3), ph%zen(3)
      do j = ph%zst(2), ph%zen(2)
        do i = ph%zst(1), ph%zen(1)
          write(*,'(A, 3I4.1, 2ES13.5)') 'After B-FFT', k, j, i, rhs0(i, j, k), rhs(i, j,k)
        end do
      end do
    end do
    deallocate(rhs0)
#endif
    
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

#ifdef DEBUGFFT
    real(wp), allocatable :: rhs0(:, :, :)

    allocate (rhs0 (ph%xsz(1), ph%xsz(2), ph%xsz(3)))
    rhs0 = rhs
#endif

#ifdef DEBUGFFT
! cw1 : x pencil format
    ! do k = sp%xst(3), sp%xen(3)
    !   do j = sp%xst(2), sp%xen(2)
    !     do i = sp%xst(1), sp%xen(1)
    !       write(*,'(A, 3I5.1, 1ES13.5)') 'input', k, j, i, rhs(i,j,k)
    !     end do
    !   end do
    ! end do
#endif 

!_______________________________________________________________________________
! in z-pencil format
! step 1 : change to x-pencil format
!_______________________________________________________________________________
    call transpose_z_to_y(rhs,     rw_ypen, ph)
    call transpose_y_to_x(rw_ypen, rw_xpen, ph)
!_______________________________________________________________________________
! in x-pencil format
! step 2 : nput data re-organisation to get a periodic sequence in physical domain
!_______________________________________________________________________________
    do k = ph%xst(3), ph%xen(3)
      do j = ph%xst(2), ph%xen(2)
        do i = 1, nw(1)/2
          rw_recons_xpen(i, j, k) = rw_xpen( 2 * i - 1, j,  k)
        end do
        do i = nw(1)/2 + 1, nw(1)
          rw_recons_xpen(i, j, k) = rw_xpen( 2 * nw(1) - 2 * i + 2, j, k)
        end do
      end do
    end do 
!_______________________________________________________________________________
! in x-pencil format
! step 3 : back to z-pencil format
!_______________________________________________________________________________
    call transpose_x_to_y(rw_recons_xpen, rw_ypen, ph)
    call transpose_y_to_z(rw_ypen,        rhs,     ph)
!_______________________________________________________________________________
! initialize FFT
!_______________________________________________________________________________
    if (.not. fft_initialised) then
      call decomp_2d_fft_init(PHYSICAL_IN_Z, nw(1), nw(2), nw(3))
      fft_initialised = .true.
    end if
!_______________________________________________________________________________
! in z-pencil format
! compute r2c transform, forward FFT
! input : rhs     : default z pencil format
! output: cw_xpen : x pencil format
!_______________________________________________________________________________
    call decomp_2d_fft_3d(rhs, cw_xpen)
!_______________________________________________________________________________
! fft normalisation
!_______________________________________________________________________________
    cw_xpen = cw_xpen / real(nw(1), kind=wp) / &
                        real(nw(2), kind=wp) / &
                        real(nw(3), kind=wp)
#ifdef DEBUGFFT
! cw1 : x pencil format
    ! do k = sp%xst(3), sp%xen(3)
    !   do j = sp%xst(2), sp%xen(2)
    !     do i = sp%xst(1), sp%xen(1)
    !       write(*,'(A, 3I5.1, 2ES13.5)') 'After F-FFT', k, j, i, cw_xpen(i,j,k)
    !     end do
    !   end do
    ! end do
#endif 
!_______________________________________________________________________________
! cwx in x pencil format
! Reconstruct FFT(rhs) from the above FFT(rhsx) in spectral domain 
!_______________________________________________________________________________
    do k = sp%xst(3),sp%xen(3)
      do j = sp%xst(2),sp%xen(2)
        cw_recons_xpen(1, j, k) = cw_xpen(1, j, k)
        do i = 2, nw(1)
          cw_recons_xpen(i, j, k) = cw_xpen(i, j, k) + cw_xpen(nw(1) - i + 2, j, k)
        end do
      end do
    end do
!_______________________________________________________________________________
! Fourier domain calculation
!_______________________________________________________________________________
    do k = sp%xst(3), sp%xen(3)
      do j = sp%xst(2), sp%xen(2)
        do i = sp%xst(1), sp%xen(1)
          if ( ( dabs( real(t2xyz(i, j, k))) < TRUNCERR ) .AND. &
               ( dabs(aimag(t2xyz(i, j, k))) < TRUNCERR )) then
            cw_recons_xpen(i, j, k) = cmplx(ZERO, ZERO)
          else 
            cw_recons_xpen(i, j, k) = cw_recons_xpen(i, j, k) / t2xyz(i, j, k)
          end if
        end do
      end do
    end do
!_______________________________________________________________________________
! Reconstruct the func(FFT(rhsx)) from the above func(FFT(rhs))
!_______________________________________________________________________________
    do k = sp%xst(3),sp%xen(3)
      do j = sp%xst(2),sp%xen(2)
        cw_xpen(1, j, k) = cw_recons_xpen(1, j, k)
        do i = 2, nw(1)
          cw_xpen(i, j, k) = cw_recons_xpen(i, j, k) + cw_recons_xpen(nw(1) - i + 2, j, k)
        end do
      end do
    end do
!_______________________________________________________________________________
! compute c2r transform, inverse FFT
! input: cw in x-pencil format
! output: rhs in z pencil
!_______________________________________________________________________________
    call decomp_2d_fft_3d(cw_xpen,rw_zpen)
!_______________________________________________________________________________
! in z-pencil format
! back to x-pencil format
!_______________________________________________________________________________
    call transpose_z_to_y(rw_zpen, rw_ypen,        ph)
    call transpose_y_to_x(rw_ypen, rw_recons_xpen, ph)
!_______________________________________________________________________________
! reorganize the output physical data structure
! rhs in x-pencil format
!_______________________________________________________________________________
    do k = ph%xst(3), ph%xen(3)
      do j = ph%xst(2), ph%xen(2)
        do i = 1, nw(1)/2
          rw_recons_xpen( 2 * i - 1, j, k) = rw_xpen(i, j, k)
        end do
        do i = 1, nw(1)/2
          rw_recons_xpen( 2 * i,     j, k) = rw_xpen(nw(1) - i + 1, j, k)
        end do
      end do
    end do
!_______________________________________________________________________________
! in x-pencil format
! back to z-pencil format
!_______________________________________________________________________________
    call transpose_x_to_y(rw_recons_xpen, rw_ypen, ph)
    call transpose_y_to_z(rw_ypen,        rhs,     ph)

#ifdef DEBUGFFT
    do k = ph%xst(3), ph%xen(3)
      do j = ph%xst(2), ph%xen(2)
        do i = ph%xst(1), ph%xen(1)
          write(*,'(A, 3I4.1, 2ES13.5)') 'After B-FFT', k, j, i, rhs0(i, j, k), rhs(i, j,k)
        end do
      end do
    end do

    deallocate(rhs0)
#endif
    return
  end subroutine Solve_poisson_100

!===============================================================================
  subroutine Solve_poisson_010(rhs)
    use decomp_2d_fft!, only : decomp_2d_fft_3d
    use parameters_constant_mod
    use input_general_mod
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT) :: rhs
    integer(4) :: i, j, k
    real(wp) :: cwRe1, cwRe2, cwIm1, cwIm2
    real(wp) :: aRe1, bRe1, aRe2, bRe2, &
                aIm1, bIm1, aIm2, bIm2

#ifdef DEBUGFFT
    real(wp), allocatable :: rhs0(:, :, :)

    allocate (rhs0 (ph%xsz(1), ph%xsz(2), ph%xsz(3)))
    rhs0 = rhs
#endif

!_______________________________________________________________________________
! in z-pencil format
! change to y-pencil format
!_______________________________________________________________________________
    call transpose_z_to_y(rhs, rw_ypen, ph)
!_______________________________________________________________________________
! in y-pencil format
! nput data re-organisation to get a periodic sequence in physical domain
!_______________________________________________________________________________
    do k = ph%xst(3), ph%xen(3)
      do i = ph%xst(1), ph%xen(1)
        do j = 1, nw(2)/2
          rw_recons_ypen(i, j, k) = rw_ypen(i, 2 * j - 1, k)
        end do
        do j = nw(2)/2 + 1, nw(2)
          rw_recons_ypen(i, j, k) = rw_ypen(i, 2 * nw(2) - 2 * j + 2, k)
        end do
      end do
    end do 
!_______________________________________________________________________________
! in y-pencil format
! change to z-pencil format
!_______________________________________________________________________________
    call transpose_y_to_z(rw_recons_ypen, rhs, ph)
!_______________________________________________________________________________
! initialize FFT
!_______________________________________________________________________________
    if (.not. fft_initialised) then
      call decomp_2d_fft_init(PHYSICAL_IN_Z, nw(1), nw(2), nw(3))
      fft_initialised = .true.
    end if
!_______________________________________________________________________________
! in z-pencil format
! compute r2c transform, forward FFT
! input : rhs     : default z pencil format
! output: cw_xpen : x pencil format
!_______________________________________________________________________________
    call decomp_2d_fft_3d(rhs,cw_xpen)
!_______________________________________________________________________________
! fft normalisation
!_______________________________________________________________________________
    cw_xpen = cw_xpen / real(nw(1), kind=wp) / &
                        real(nw(2), kind=wp) / &
                        real(nw(3), kind=wp)
!_______________________________________________________________________________
! cw in x-pencil format
! change to y-pencil format
!_______________________________________________________________________________
    call transpose_x_to_y(cw_xpen, cw_ypen, sp)         
!_______________________________________________________________________________
! cwx in y pencil format
! Reconstruct FFT(rhs) from the above FFT(rhsx) in spectral domain 
!_______________________________________________________________________________
    do k = sp%xst(3),sp%xen(3)
      do i = sp%xst(1),sp%xen(1)
        cw_recons_ypen(i, 1, k) = cw_ypen(i, 1, k)
        do j = 2, nw(2)
          cw_recons_ypen(i, j, k) = cw_ypen(i, j, k) + cw_ypen(i, nw(2) - j + 2, k)
        end do
      end do
    end do
!_______________________________________________________________________________
! Fourier domain calculation
!_______________________________________________________________________________
    if(istret == ISTRET_NO) then
      do k = sp%yst(3), sp%yen(3)
        do j = sp%yst(2), sp%yen(2)
          do i = sp%yst(1), sp%yen(1)
            if ( ( dabs( real(t2xyz(i, j, k))) < TRUNCERR ) .AND. &
                 ( dabs(aimag(t2xyz(i, j, k))) < TRUNCERR )) then
              cw_recons_ypen(i, j, k) = cmplx(ZERO, ZERO)
            else 
              cw_recons_ypen(i, j, k) = cw_recons_ypen(i, j, k) / t2xyz(i, j, k)
            end if
          end do
        end do
      end do
    else 
      ! to add ...
    end if
!_______________________________________________________________________________
! Reconstruct the func(FFT(rhs)) from the above func(FFT(rhs))
!_______________________________________________________________________________
    do k = sp%xst(3),sp%xen(3)
      do i = sp%xst(1),sp%xen(1)
        cw_ypen(i, 1, k) = cw_recons_ypen(i, 1, k)
        do j = 2, nw(2)
          cw_ypen(i, j, k) = cw_recons_ypen(i, j, k) + cw_recons_ypen(i, nw(2) - j + 2, k)
        end do
      end do
    end do
!_______________________________________________________________________________
! cw in y-pencil format
! change to x-pencil format
!_______________________________________________________________________________
    call transpose_y_to_x(cw_ypen, cw_xpen, sp)
!_______________________________________________________________________________
! compute c2r transform, inverse FFT
! input: cw in y-pencil format
! output: rhs in z pencil
!_______________________________________________________________________________
    call decomp_2d_fft_3d(cw_xpen,rw_zpen)
!_______________________________________________________________________________
! in z-pencil format
! back to y-pencil format
!_______________________________________________________________________________
    call transpose_z_to_y(rw_zpen, rw_recons_ypen,        ph)
!_______________________________________________________________________________
! reorganize the output physical data structure
! rhs in y-pencil format
!_______________________________________________________________________________
    do k = ph%xst(3), ph%xen(3)
      do i = ph%xst(1), ph%xen(1)
        do j = 1, nw(2)/2
          rw_ypen(i, 2 * j - 1, k) = rw_recons_ypen(i, j, k)
        end do
        do j = 1, nw(2)/2
          rw_ypen(i, 2 * j, k) = rw_recons_ypen(i, nw(2) - j + 1, k)
        end do
      end do
    end do

    call transpose_y_to_z(rw_ypen, rhs, ph)

#ifdef DEBUGFFT
    do k = ph%xst(3), ph%xen(3)
      do j = ph%xst(2), ph%xen(2)
        do i = ph%xst(1), ph%xen(1)
          write(*,'(A, 3I4.1, 2ES13.5)') 'After B-FFT', k, j, i, rhs0(i, j, k), rhs(i, j,k)
        end do
      end do
    end do

    deallocate(rhs0)
#endif
    return
  end subroutine Solve_poisson_010

!===============================================================================
  subroutine Solve_poisson_001(rhs)
    use decomp_2d_fft!, only : decomp_2d_fft_3d
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT) :: rhs

    return
  end subroutine Solve_poisson_001

!===============================================================================
  subroutine Solve_poisson_110(rhs)
    use decomp_2d_fft!, only : decomp_2d_fft_3d
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT) :: rhs

    return
  end subroutine Solve_poisson_110
!===============================================================================
  subroutine Solve_poisson_101(rhs)
    use decomp_2d_fft!, only : decomp_2d_fft_3d
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT) :: rhs

    return
  end subroutine Solve_poisson_101
!===============================================================================
  subroutine Solve_poisson_011(rhs)
    use decomp_2d_fft!, only : decomp_2d_fft_3d
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT) :: rhs

    return
  end subroutine Solve_poisson_011
!===============================================================================
  subroutine Solve_poisson_111(rhs)
    use decomp_2d_fft!, only : decomp_2d_fft_3d
    implicit none
    real(wp), dimension(:,:,:), intent(INOUT) :: rhs

    return
  end subroutine Solve_poisson_111
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
    
    write(*, *) ' Test Poisson Solver >>'

    allocate(rhsphi(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3))); rhsphi = ZERO
    ! d^2(phi_i)/dx2_i = rhs
    do i = ph%zst(1),ph%zen(1)
      do j = ph%zst(2),ph%zen(2)
        do k = ph%zst(3),ph%zen(3)
          x = domain%h(1)*(real(i - 1, WP))
          y = domain%h(2)*(real(j - 1, WP))
          z = domain%h(3)*(real(k - 1, WP))
          rhsphi(i, j, k) = -dsin(x)
          !-dsin(TWO*x + TWO*y + TWO*z)*FOUR &
          !                  -dsin(TWO*x + TWO*y + TWO*z)*FOUR &
          !                  -dsin(TWO*x + TWO*y + TWO*z)*FOUR
        end do
      end do
    end do
    
    call Solve_poisson(rhsphi)

    nn = 0
    do k = ph%zst(3),ph%zen(3)
      do j = ph%zst(2),ph%zen(2)
        do i = ph%zst(1),ph%zen(1)
          x = domain%h(1)*(real(i - 1, WP))
          y = domain%h(2)*(real(j - 1, WP))
          z = domain%h(3)*(real(k - 1, WP))

          solution = dsin(x)
          !dsin(TWO*x + TWO*y + TWO*z)

          write(*, *) k, j, i, solution, rhsphi(i,j,k), dabs(rhsphi(i,j,k)-solution)
        end do
      end do
    end do
    deallocate(rhsphi)
    stop
  end subroutine

end module poisson_mod
