module possion_module
  implicit none
  use decomp_2d, only : mytype
  private 
  real(wp), allocatable, dimension(:) :: t2x, t2y, t2z
  real(wp), save, allocatable, dimension(:, :, :) :: t2xyz

!_______________________________________________________________________________
! boundary conditions
!_______________________________________________________________________________
integer, save :: bcx, bcy, bcz
!_______________________________________________________________________________
! decomposition object for physical space (ph) and spectral space (sp)
!_______________________________________________________________________________
  type(DECOMP_INFO), save :: ph
  type(DECOMP_INFO), save :: sp
!_______________________________________________________________________________
! work arrays, 
! naming convention: cw (complex); rw (real); 
!                    b =     ; c = 
!                    1 = X-pencil; 2 = Y-pencil; 3 = Z-pencil
!_______________________________________________________________________________

  complex(mytype), allocatable, dimension(:, :, :) :: cw1
!_______________________________________________________________________________
! FFT library only needs to be initialised once
!_______________________________________________________________________________
  logical, save :: fft_initialised = .false.
!_______________________________________________________________________________
! interface for basic poisson solver
!_______________________________________________________________________________
  ABSTRACT INTERFACE
    SUBROUTINE poisson_xxx(rhs)
      use decomp_2d, only : mytype
      real(mytype), dimension(:, :, :), intent(INOUT) :: rhs
    END SUBROUTINE poisson_xxx
  END INTERFACE

  PROCEDURE (poisson_xxx), POINTER :: poisson

  public :: decomp_2d_poisson_init, &
            decomp_2d_poisson_finalize, &
            poisson

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
subroutine decomp_2d_poisson_init(d)
    use udf_type_mod,            only: t_domain
    use parameters_constant_mod, only: ZERO
    use decomp_2d,               only: nx_global, ny_global, nz_global
    implicit none
    type(t_domain), intent(in) :: d

    integer(4) :: nx, ny, nz, i

!_______________________________________________________________________________
! set up boundary flags for periodic b.c.
!_______________________________________________________________________________
    if ( is_periodic(1) ) then
      bcx = 0
    else
      bcx = 1
    end if

    if ( is_periodic(2) ) then
      bcy = 0
    else
      bcy = 1
    end if

    if ( is_periodic(3) ) then
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
      poisson => poisson_000
    else if (bcx == 1 .and. bcy == 0 .and. bcz == 0) then
      poisson => poisson_100
    else if (bcx == 0 .and. bcy == 1 .and. bcz == 0) then
      poisson => poisson_010
    else if (bcx == 1 .and. bcy == 1) then   ! 110 & 111
      poisson => poisson_11x
    else
      stop 'boundary condition not supported'
    end if

!_______________________________________________________________________________
! size of working array
! pressure-grid having 1 fewer point for non-periodic directions
!_______________________________________________________________________________
    nx = nx_global
    ny = ny_global
    nz = nz_global

    if (bcx == 1) nx = nx - 1 ! to check ?
    if (bcy == 1) ny = ny - 1
    if (bcz == 1) nz = nz - 1

!_______________________________________________________________________________
! decomp_info_init to prepare two working decomp_type
!_______________________________________________________________________________
    call decomp_info_init(nx, ny, nz,       ph)
    call decomp_info_init(nx, ny, nz/2 + 1, sp)
!_______________________________________________________________________________
! allocate working arrays
! complex and real part of working array in x-pencil
!_______________________________________________________________________________
    allocate (cw1 (sp%xst(1) : sp%xen(1), &
                   sp%xst(2) : sp%xen(2), &
                   sp%xst(3) : sp%xen(3) ) )
!_______________________________________________________________________________
! prepare the transformation \hat{f''}_l = \hat{f}_l * t2x
!_______________________________________________________________________________
    allocate (t2x(nx)) ;  t2x = ZERO
    allocate (t2y(ny)) ;  t2y = ZERO
    allocate (t2z(nz)) ;  t2z = ZERO
    call Transform_2nd_derivative_spectral_1d(d%bc(1, 1), d%nc(1), d%h(1), t2x)
    call Transform_2nd_derivative_spectral_1d(d%bc(1, 2), d%nc(2), d%h(2), t2y)
    call Transform_2nd_derivative_spectral_1d(d%bc(1, 3), d%nc(3), d%h(3), t2z)
    allocate (t2xyz(nx, ny, nz)) ;  t2z = ZERO
    call Transform_2nd_derivative_spectral_3d(nx, ny, nz, t2x, t2y, t2z, t2xyz)

    return
  end subroutine decomp_2d_poisson_init


  subroutine Transform_2nd_derivative_spectral_1d(ibc, nx, dx, t2x)
    use operations, only : d2fC2C, d2rC2C
    use input_general_mod, only : IBC_PERIODIC
    implicit none

    integer(wp), intent(in) :: nx
    integer(wp), intent(in) :: ibc
    real(WP),    intent(in) :: dx
    real(WP),    intent(out) :: t2x

    if(ibc == IBC_PERIODIC)
    ! below is for 3 periodic only. due to enrich
    alpha = d2fC2C(3, 1, ibc)
    a     = d2rC2C(3, 1, ibc)
    b     = d2rC2C(3, 2, ibc) * FOUR

    do i = 1, nx
      w = TWO * PI / REAL( nx * (j - 1), WP)
      cosw = cos_wp(w)
      t2x(i) = b * cosw * cosw + TWO * a * cosw - TWO * a - b
      t2x(i) = t2x(i) / (ONE + TWO * alpha * cosw) / dx / dx
    end do
    
    end if
    return
  end subroutine Transform_2nd_derivative_spectral_1d

  subroutine poisson_000(rhs)
    implicit none
    ! right-hand-side of Poisson as input
    ! solution of Poisson as output
    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

    integer, dimension(3) :: fft_start, fft_end, fft_size

    complex(mytype) :: xyzk

    complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
    complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2


    real(mytype) :: tmp1, tmp2,x ,y, z

    integer :: nx,ny,nz, i,j,k

    nx = nx_global
    ny = ny_global
    nz = nz_global

    if (.not. fft_initialised) then
       call decomp_2d_fft_init(PHYSICAL_IN_Z)
       fft_initialised = .true.
    end if

    ! compute r2c transform 
    call decomp_2d_fft_3d(rhs,cw1)

    ! normalisation
    cw1 = cw1 / real(nx, kind=mytype) /real(ny, kind=mytype) &
         / real(nz, kind=mytype)

    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)

            cw1(i,j,k)=cmplx( real(cw1(i,j,k), kind=mytype) / (-tmp1), &
            aimag(cw1(i,j,k))/(-tmp2), kind=mytype)

             !Print result in spectal space after Poisson
             !     if (abs(out(i,j,k)) > 1.0e-4) then
             !        write(*,*) 'AFTER',i,j,k,out(i,j,k),xyzk
             !     end if

             ! post-processing backward

             ! POST PROCESSING IN Z
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)-tmp2*az(k), &
                  -tmp2*bz(k)-tmp1*az(k), kind=mytype)

             ! POST PROCESSING IN Y
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*by(j)+tmp2*ay(j), &
                  tmp2*by(j)-tmp1*ay(j), kind=mytype)
             if (j.gt.(ny/2+1)) cw1(i,j,k)=-cw1(i,j,k)

             ! POST PROCESSING IN X
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bx(i)+tmp2*ax(i), &
                  -tmp2*bx(i)+tmp1*ax(i), kind=mytype)
             if (i.gt.(nx/2+1)) cw1(i,j,k)=-cw1(i,j,k)

          end do
       end do
    end do

    ! compute c2r transform
    call decomp_2d_fft_3d(cw1,rhs)

    !   call decomp_2d_fft_finalize

    return
  end subroutine poisson_000

end module