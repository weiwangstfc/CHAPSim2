module possion_module
  use decomp_2d, only : DECOMP_INFO, mytype, nrank
  implicit none

  private

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16_mytype
#else
  real(mytype), parameter :: epsilon = 1.e-8_mytype
#endif

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
! store sine/cosine factors
!_______________________________________________________________________________
  real(mytype), save, allocatable, dimension(:) :: az, bz
  real(mytype), save, allocatable, dimension(:) :: ay, by
  real(mytype), save, allocatable, dimension(:) :: ax, bx
!_______________________________________________________________________________
! wave numbers: kxyz
! wave numbers for stretching in a pentadiagonal matrice: a, a2, a3
!_______________________________________________________________________________
  complex(mytype), save, allocatable, dimension(:, :, :)   :: kxyz
  complex(mytype), save, allocatable, dimension(:, :, :,:) :: a, a2, a3
!_______________________________________________________________________________
! work arrays, 
! naming convention: cw (complex); rw (real); 
!                    b =     ; c = 
!                    1 = X-pencil; 2 = Y-pencil; 3 = Z-pencil
!_______________________________________________________________________________
  real   (mytype), allocatable, dimension(:, :, :) :: rw1, rw1b, rw2, rw2b, rw3
  complex(mytype), allocatable, dimension(:, :, :) :: cw1, cw1b, cw2, cw2b, cw22, cw2c
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

public :: Initialize_decomp_poisson, &
          decomp_2d_poisson_finalize, &
          poisson

contains
!===============================================================================
!===============================================================================
!> \brief To asign sine and cose factors
!>
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]    nsz           working array size
!> \param[in]    bc            b.c. flags 
!> \param[out]   afsin         sine factors
!> \param[out]   bfsin         cosine factors
!_______________________________________________________________________________
  subroutine Calculate_sine_cosine_factors(afsin, bfcos, nsz, bc)
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
          bfsin(i) = cos_wp( real(i - 1, kind = mytype) * PI / &
                             real(nsz,   kind = mytype) )
        end do

    else if (bc == 1) then

        do i = 1, nsz
          afsin(i) = sin_wp( real(i - 1, kind = mytype) * PI / TWO / &
                             real(nsz,   kind = mytype) )
          bfsin(i) = cos_wp( real(i - 1, kind = mytype) * PI / TWO / &
                             real(nsz,   kind = mytype))
        end do
    else
    end if

    return
  end subroutine Calculate_sine_cosine_factors
!===============================================================================
!===============================================================================
!> \brief To asign sine and cose factors
!>
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]    nsz           working array size
!> \param[in]    bc            b.c. flags 
!> \param[out]   afsin         sine factors
!> \param[out]   bfsin         cosine factors
!_______________________________________________________________________________
  subroutine waves ()
    !
    !***********************************************************

    USE derivX 
    USE derivY 
    USE derivZ 
    USE param
    USE decomp_2d
    USE variables
    use decomp_2d_fft

    implicit none

    integer :: i,j,k
    real(mytype) :: w,wp,w1,w1p 
    complex(mytype) :: xyzk
    complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
    complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2,tmp1,tmp2,tmp3
    complex(mytype) :: tmp4,tmp5,tmp6

    xkx(:)=0. ; xk2(:)=0. ; yky(:)=0. ; yk2(:)=0.
    zkz(:)=0. ; zk2(:)=0.

    !WAVE NUMBER IN X
    if (bcx==0) then
       do i=1,nx/2+1
          w=2.*pi*(i-1)/nx
          wp=acix6*2.*dx*sin(w/2.)+(bcix6*2.*dx)*sin(3./2.*w)
          wp=wp/(1.+2.*alcaix6*cos(w))
          xkx(i)=cmplx( nx*wp/xlx,      nx*wp/xlx,     kind=mytype)
          exs(i)=cmplx( nx*w/xlx,       nx*w/xlx,      kind=mytype)
          xk2(i)=cmplx((nx*wp/xlx)**2, (nx*wp/xlx)**2, kind=mytype)
       enddo
       do i=nx/2+2,nx
          xkx(i)=xkx(nx-i+2)
          exs(i)=exs(nx-i+2)
          xk2(i)=xk2(nx-i+2)
       enddo
    else
       do i=1,nx
          w=2.*pi*0.5*(i-1)/nxm
          wp=acix6*2.*dx*sin(w/2.)+(bcix6*2.*dx)*sin(3./2.*w)
          wp=wp/(1.+2.*alcaix6*cos(w))
          xkx(i)=cmplx( nxm*wp/xlx,     nxm*wp/xlx,     kind=mytype)
          exs(i)=cmplx( nxm*w/xlx,      nxm*w/xlx,      kind=mytype)
          xk2(i)=cmplx((nxm*wp/xlx)**2,(nxm*wp/xlx)**2, kind=mytype)
       enddo
       xkx(1)=0.
       exs(1)=0.
       xk2(1)=0.
    endif

    !WAVE NUMBER IN Y
    if (bcy==0) then
       do j=1,ny/2+1
          w=2.*pi*(j-1)/ny
          wp=aciy6*2.*dy*sin(w/2.)+(bciy6*2.*dy)*sin(3./2.*w)
          wp=wp/(1.+2.*alcaiy6*cos(w))
          if (istret==0) yky(j)=cmplx(ny*wp/yly,ny*wp/yly, kind=mytype)
          if (istret.ne.0) yky(j)=cmplx(ny*wp,ny*wp, kind=mytype)
          eys(j)=cmplx(ny*w/yly,ny*w/yly, kind=mytype)
          yk2(j)=cmplx((ny*wp/yly)**2,(ny*wp/yly)**2, kind=mytype)
       enddo
       do j=ny/2+2,ny
          yky(j)=yky(ny-j+2)
          eys(j)=eys(ny-j+2)
          yk2(j)=yk2(ny-j+2)
       enddo
    else
       do j=1,ny
          w=2.*pi*0.5*(j-1)/nym
          wp=aciy6*2.*dy*sin(w/2.)+(bciy6*2.*dy)*sin(3./2.*w)
          wp=wp/(1.+2.*alcaiy6*cos(w))
          if (istret==0) yky(j)=cmplx(nym*wp/yly,nym*wp/yly, kind=mytype)
          if (istret.ne.0) yky(j)=cmplx(nym*wp,nym*wp, kind=mytype)
          eys(j)=cmplx(nym*w/yly,nym*w/yly, kind=mytype)
          yk2(j)=cmplx((nym*wp/yly)**2,(nym*wp/yly)**2, kind=mytype)
       enddo
       yky(1)=0.
       eys(1)=0.
       yk2(1)=0.
    endif

    !WAVE NUMBER IN Z
    if (bcz==0) then
       do k=1,nz/2+1
          w=2.*pi*(k-1)/nz
          wp=aciz6*2.*dz*sin(w/2.)+(bciz6*2.*dz)*sin(3./2.*w)
          wp=wp/(1.+2.*alcaiz6*cos(w))
          zkz(k)=cmplx(nz*wp/zlz,nz*wp/zlz, kind=mytype)
          ezs(k)=cmplx(nz*w/zlz,nz*w/zlz, kind=mytype)
          zk2(k)=cmplx((nz*wp/zlz)**2,(nz*wp/zlz)**2, kind=mytype)
       enddo
    else
       do k=1,nz/2+1
          w=2.*pi*0.5*(k-1)/nzm
          w1=2.*pi*0.5*(nzm-k+1)/nzm
          wp=aciz6*2.*dz*sin(w/2.)+(bciz6*2.*dz)*sin(3./2.*w)
          wp=wp/(1.+2.*alcaiz6*cos(w))
          w1p=aciz6*2.*dz*sin(w1/2.)+(bciz6*2.*dz)*sin(3./2.*w1)
          w1p=w1p/(1.+2.*alcaiz6*cos(w1))     
          zkz(k)=cmplx(nzm*wp/zlz,-nzm*w1p/zlz, kind=mytype)
          ezs(k)=cmplx(nzm*w/zlz,nzm*w1/zlz, kind=mytype)
          zk2(k)=cmplx((nzm*wp/zlz)**2,(nzm*w1p/zlz)**2, kind=mytype)
       enddo
    endif
    !
    !if (nrank==0) then
    !   do i=1,nx
    !      print *,i,ezs(i)
    !   enddo
    !endif
    !stop

    if ((bcx==0).and.(bcz==0).and.bcy.ne.0) then
       do k = sp%yst(3), sp%yen(3)
          do j = sp%yst(2), sp%yen(2)
             do i = sp%yst(1), sp%yen(1)
                xtt=cmplx( &
                    (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+ &
                     cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+ &
                     dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)), &
                    (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+ &
                     cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+ &
                     dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)), &
                     kind=mytype)
                ytt=cmplx(&
                    (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+ &
                     ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+ &
                     diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)),&
                    (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+ &
                     ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+ &
                     diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)), &
                     kind=mytype)
                ztt=cmplx(&
                    (biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+ &
                     ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)+ &
                     diciz6*2.*cos(real(ezs(k), kind=mytype)*7.*dz/2.)), &
                    (biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+ &
                     ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)+ &
                     diciz6*2.*cos(real(ezs(k), kind=mytype)*7.*dz/2.)), &
                     kind=mytype)
                xtt1=cmplx(&
                     (aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)),&
                     (aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)), &
                     kind=mytype)
                ytt1=cmplx(&
                     (aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)),&
                     (aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)), &
                     kind=mytype)
                ztt1=cmplx(&
                     (aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)),&
                     (aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)), &
                     kind=mytype)
                xt1=cmplx(&
                     (1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)),&
                     (1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)), &
                     kind=mytype)
                yt1=cmplx(&
                     (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)),&
                     (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)), &
                     kind=mytype)
                zt1=cmplx(&
                     (1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)),&
                     (1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)),&
                    kind=mytype)
                xt2= xk2(i) * ( ( ((ytt1+ytt)/yt1) * ((ztt1+ztt)/zt1) )**2 )
                yt2= yk2(j) * ( ( ((xtt1+xtt)/xt1) * ((ztt1+ztt)/zt1) )**2 )
                zt2= zk2(k) * ( ( ((xtt1+xtt)/xt1) * ((ytt1+ytt)/yt1) )**2 )
                xyzk=xt2+yt2+zt2
                kxyz(i,j,k)=xyzk
                !   print *,i,j,k, kxyz(i,j,k)
             enddo
          enddo
       enddo
    else
       if (bcz==0) then
          do k = sp%xst(3),sp%xen(3)
             do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                   xtt=cmplx((bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
                        cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+&
                        dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)),&
                        (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
                        cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+&
                        dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)), kind=mytype)
                   ytt=cmplx((biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
                        ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+&
                        diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)),&
                        (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
                        ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+&
                        diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)), kind=mytype)
                   ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
                        ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)+&
                        diciz6*2.*cos(real(ezs(k), kind=mytype)*7.*dz/2.)),&
                        (biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
                        ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)+&
                        diciz6*2.*cos(real(ezs(k), kind=mytype)*7.*dz/2.)), kind=mytype)
                   xtt1=cmplx((aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)),&
                        (aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)), kind=mytype)
                   ytt1=cmplx((aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)),&
                        (aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)), kind=mytype)
                   ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)),&
                        (aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)), kind=mytype)
                   xt1=cmplx((1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)),&
                        (1.+2.*ailcaix6*cos(real(exs(i))*dx)), kind=mytype)
                   yt1=cmplx((1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)),&
                        (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)), kind=mytype)
                   zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)),&
                        (1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)), kind=mytype)
                   xt2=xk2(i)*((((ytt1+ytt)/yt1)*((ztt1+ztt)/zt1))**2)
                   yt2=yk2(j)*((((xtt1+xtt)/xt1)*((ztt1+ztt)/zt1))**2)
                   zt2=zk2(k)*((((xtt1+xtt)/xt1)*((ytt1+ytt)/yt1))**2)
                   xyzk=xt2+yt2+zt2
                   kxyz(i,j,k)=xyzk
                   !   print *,i,j,k, kxyz(i,j,k)
                enddo
             enddo
          enddo
       else
          do k = sp%xst(3),sp%xen(3)
             do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)  
                   xtt=cmplx((bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
                        cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+&
                        dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)),&
                        (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
                        cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+&
                        dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)), kind=mytype)
                   ytt=cmplx((biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
                        ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+&
                        diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)),&
                        (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
                        ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+&
                        diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)), kind=mytype)
                   !
                   ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
                        ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)+&
                        diciz6*2.*cos(real(ezs(k), kind=mytype)*7.*dz/2.)),&
                        (biciz6*2.*cos(aimag(ezs(k))*3.*dz/2.)+&
                        ciciz6*2.*cos(aimag(ezs(k))*5.*dz/2.)+&
                        diciz6*2.*cos(aimag(ezs(k))*7.*dz/2.)), kind=mytype)
                   !
                   xtt1=cmplx((aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)),&
                        (aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)), kind=mytype)
                   ytt1=cmplx((aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)),&
                        (aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)), kind=mytype)
                   !
                   ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)),&
                        (aiciz6*2.*cos(aimag(ezs(k))*dz/2.)), kind=mytype)
                   !
                   xt1=cmplx((1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)),&
                        (1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)), kind=mytype)
                   yt1=cmplx((1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)),&
                        (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)), kind=mytype)
                   zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)),&
                        (1.+2.*ailcaiz6*cos(aimag(ezs(k))*dz)), kind=mytype)

                   tmp1=cmplx(real(ztt1+ztt, kind=mytype)/real(zt1, kind=mytype),&
                        aimag(ztt1+ztt)/aimag(zt1), kind=mytype)
                   tmp2=cmplx(real(ytt1+ytt, kind=mytype)/real(yt1, kind=mytype),&
                        real(ytt1+ytt, kind=mytype)/real(yt1, kind=mytype), kind=mytype)
                   tmp3=cmplx(real(xtt1+xtt, kind=mytype)/real(xt1, kind=mytype),&
                        real(xtt1+xtt, kind=mytype)/real(xt1, kind=mytype), kind=mytype)

                   tmp4=cmplx((real(tmp1, kind=mytype)*real(tmp2, kind=mytype))**2,(aimag(tmp1)*aimag(tmp2))**2, kind=mytype)
                   tmp5=cmplx((real(tmp1, kind=mytype)*real(tmp3, kind=mytype))**2,(aimag(tmp1)*aimag(tmp3))**2, kind=mytype)
                   tmp6=cmplx((real(tmp3, kind=mytype)*real(tmp2, kind=mytype))**2,(aimag(tmp3)*aimag(tmp2))**2, kind=mytype)

                   tmp1=cmplx(real(tmp4, kind=mytype)*real(xk2(i), kind=mytype),aimag(tmp4)*aimag(xk2(i)), kind=mytype)
                   tmp2=cmplx(real(tmp5, kind=mytype)*real(yk2(j), kind=mytype),aimag(tmp5)*aimag(yk2(j)), kind=mytype)
                   tmp3=cmplx(real(tmp6, kind=mytype)*real(zk2(k), kind=mytype),aimag(tmp6)*aimag(zk2(k)), kind=mytype)

                   xyzk=tmp1+tmp2+tmp3
                   kxyz(i,j,k)=xyzk
                   !         print *,i,j,k,zt1,yt1
                enddo
             enddo
          enddo
       endif
    endif


    !          do k=1,1!nz
    !          do j=1,ny
    !          do i=1,1!!nx
    !             print *,j,a(i,j,k,3),kxyz(i,j,k)
    !          enddo
    !          enddo
    !          enddo

  end subroutine waves
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
  subroutine Initialize_decomp_poisson(is_periodic)
    use udf_type_mod, only : t_domain
    use parameters_constant_mod, only : ZERO
    use decomp_2d,    only : nx_global, ny_global, nz_global
    implicit none
    logical, intent(in) :: is_periodic(3)

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

  #ifdef DEBG 
    if (nrank == 0) print *,'# Initialize_decomp_poisson start'
  #endif
!_______________________________________________________________________________
! preparing sine and cosine factors
!_______________________________________________________________________________
    allocate ( ax(nx), bx(nx) ); ax = ZERO; bx = ZERO
    allocate ( ay(ny), by(ny) ); ay = ZERO; by = ZERO
    allocate ( az(nz), bz(nz) ); az = ZERO; bz = ZERO
    call Calculate_sine_cosine_factors(ax, bx, nx, bcx)
    call Calculate_sine_cosine_factors(ay, by, ny, bcy)
    call Calculate_sine_cosine_factors(az, bz, nz, bcz)

  #ifdef DEBG 
    if (nrank .eq. 0) print *,'# Initialize_decomp_poisson decomp_info_init'
  #endif
!_______________________________________________________________________________
! decomp_info_init to prepare two working decomp_type
!_______________________________________________________________________________
    call decomp_info_init(nx, ny, nz,       ph)
    call decomp_info_init(nx, ny, nz/2 + 1, sp)

  #ifdef DEBG 
    if (nrank .eq. 0) print *,'# Initialize_decomp_poisson decomp_info_init ok'
  #endif

!_______________________________________________________________________________
! allocate working arrays
! complex and real part of working array in x-pencil
!_______________________________________________________________________________
  
  allocate (cw1 (sp%xst(1) : sp%xen(1), &
                 sp%xst(2) : sp%xen(2), &
                 sp%xst(3) : sp%xen(3) ) )

  if( bcx == 1 ) then

    allocate (rw1  ( ph%xst(1) : ph%xen(1), &
                     ph%xst(2) : ph%xen(2), &
                     ph%xst(3) : ph%xen(3) ) )

    allocate (cw1b ( sp%xst(1) : sp%xen(1), &
                     sp%xst(2) : sp%xen(2), &
                     sp%xst(3) : sp%xen(3) ) )

    allocate (rw1b ( ph%xst(1) : ph%xen(1), &
                     ph%xst(2) : ph%xen(2), &
                     ph%xst(3) : ph%xen(3) ) )

  end if
!_______________________________________________________________________________
! allocate working arrays
! complex and real part of working array in y-pencil
!_______________________________________________________________________________
  if( bcy == 1 ) then

    allocate (cw2 (sp%yst(1) : sp%yen(1), &
                   sp%yst(2) : sp%yen(2), &
                   sp%yst(3) : sp%yen(3) ) )

    allocate (rw2b(ph%yst(1) : ph%yen(1), &
                   ph%yst(2) : ph%yen(2), &
                   ph%yst(3) : ph%yen(3) ) )

    allocate (cw2b(sp%yst(1) : sp%yen(1), &
                   sp%yst(2) : sp%yen(2), &
                   sp%yst(3) : sp%yen(3) ) )

    allocate (cw22(sp%yst(1) : sp%yen(1), &
                   sp%yst(2) : sp%yen(2), &
                   sp%yst(3) : sp%yen(3) ) )
     
    allocate (cw2c(sp%yst(1) : sp%yen(1), &
                   sp%yst(2) : sp%yen(2), &
                   sp%yst(3) : sp%yen(3) ) )

  end if

  if (bcx /= 0 .or. bcy /= 0 .or. bcz /= 0) then
    allocate (rw2  ( ph%yst(1) : ph%yen(1), &
                     ph%yst(2) : ph%yen(2), &
                     ph%yst(3) : ph%yen(3) ) )
  end if
!_______________________________________________________________________________
! allocate working arrays
! complex and real part of working array in z-pencil
!_______________________________________________________________________________
  if ( bcz == 1) then  ! 
    allocate (rw3 ( ph%zsz(1), &
                    ph%zsz(2), &
                    ph%zsz(3) ) )
  end if
!_______________________________________________________________________________
! complex type wave number, check size ?
!_______________________________________________________________________________
  if(bcx==0 .and. bcy==1 .and. bcz==0) then
    allocate(kxyz(sp%yst(1) : sp%yen(1), &
                  sp%yst(2) : sp%yen(2), &
                  sp%yst(3) : sp%yen(3) ) )
  else
    allocate(kxyz(sp%xst(1) : sp%xen(1), &
                  sp%xst(2) : sp%xen(2), &
                  sp%xst(3) : sp%xen(3) ) )
  end if
!_______________________________________________________________________________
!complex type wave numbers for stretching in a pentadiagonal matrice
!_______________________________________________________________________________
  if (bcx==1 .and. bcy==1) then
    ysz = nym
  else
    ysz = ny
  end if
  allocate(a   (sp%yst(1) : sp%yen(1), &
                ny/2,                  &
                sp%yst(3) : sp%yen(3), &
                5 ) )
  allocate(a2  (sp%yst(1) : sp%yen(1), &
                ny/2,                  &
                sp%yst(3) : sp%yen(3), &
                5 ) )
  allocate(a3  (sp%yst(1) : sp%yen(1), &
                ysz,                   &
                sp%yst(3) : sp%yen(3), &
                5 ) )
  

  #ifdef DEBG 
    if (nrank .eq. 0) print *,'# Initialize_decomp_poisson before waves'
  #endif

    call waves()

  #ifdef DEBG 
    if (nrank .eq. 0) print *,'# Initialize_decomp_poisson end'
  #endif

    return
  end subroutine Initialize_decomp_poisson


end module possion_module