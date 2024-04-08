module poisson_interface_mod
  use decomp_2d
  use mpi_mod
  use parameters_constant_mod, disabled => WP!, only: zero, half, one, onepfive, two, twopfive, &
                             !        three, pi, threepfive, four, twopi, cx_one_one
  use math_mod, only: cos_prec, abs_prec, sin_prec
  use geometry_mod, only: alpha, beta
  implicit none

  integer :: istret

  integer, parameter :: IFORWARD  = 1
  integer, parameter :: IBACKWARD = -1

!----------------------------------------------------------------------------------------------------------
  real(mytype) :: xlx ! domain length
  real(mytype) :: yly ! physical domain
  real(mytype) :: zlz
!----------------------------------------------------------------------------------------------------------
  logical :: nclx ! logic, whether it is periodic bc
  logical :: ncly
  logical :: nclz
!----------------------------------------------------------------------------------------------------------
  ! below information is from incompact3d.
  ! Boundary conditions : ncl = 2 --> Dirichlet
  ! Boundary conditions : ncl = 1 --> Free-slip
  ! Boundary conditions : ncl = 0 --> Periodic
  ! l: power of 2,3,4,5 and 6
  ! if ncl = 1 or 2, --> n  = 2l+ 1
  !                  --> nm = n - 1
  !                  --> m  = n + 1
  ! If ncl = 0,      --> n  = 2*l
  !                  --> nm = n
  !                  --> m  = n + 2
  integer :: nclx1 ! boundary condition, velocity
  integer :: ncly1 
  integer :: nclz1
!----------------------------------------------------------------------------------------------------------
  integer, save :: nx ! computational node number
  integer, save :: ny
  integer, save :: nz
!----------------------------------------------------------------------------------------------------------
  integer, save :: nxm ! number of spacing 
  integer, save :: nym
  integer, save :: nzm
!----------------------------------------------------------------------------------------------------------
  real(mytype), save :: dx
  real(mytype), save :: dy ! physical grid spacing
  real(mytype), save :: dz
!---------------------------------------------------------------------------------------------------------- 
  !real(mytype) :: alpha
  !real(mytype) :: beta
!----------------------------------------------------------------------------------------------------------
  real(mytype), save :: alcaix6 
  real(mytype), save :: acix6
  real(mytype), save :: bcix6
!----------------------------------------------------------------------------------------------------------
  real(mytype), save :: alcaiy6 
  real(mytype), save :: aciy6
  real(mytype), save :: bciy6
!----------------------------------------------------------------------------------------------------------
  real(mytype), save :: alcaiz6
  real(mytype), save :: aciz6
  real(mytype), save :: bciz6
!----------------------------------------------------------------------------------------------------------
  real(mytype), save :: ailcaix6 
  real(mytype), save :: aicix6
  real(mytype), save :: bicix6
  real(mytype), save :: cicix6
  real(mytype), save :: dicix6
!----------------------------------------------------------------------------------------------------------
  real(mytype), save :: ailcaiy6 
  real(mytype), save :: aiciy6
  real(mytype), save :: biciy6
  real(mytype), save :: ciciy6
  real(mytype), save :: diciy6
!----------------------------------------------------------------------------------------------------------
  real(mytype), save :: ailcaiz6
  real(mytype), save :: aiciz6
  real(mytype), save :: biciz6
  real(mytype), save :: ciciz6
  real(mytype), save :: diciz6
!----------------------------------------------------------------------------------------------------------
  !module waves
  complex(mytype),allocatable,dimension(:), save :: zkz,zk2,ezs
  complex(mytype),allocatable,dimension(:), save :: yky,yk2,eys
  complex(mytype),allocatable,dimension(:), save :: xkx,xk2,exs

  public :: build_up_poisson_interface

contains
!==========================================================================================================
  subroutine build_up_poisson_interface(dm)
    use udf_type_mod
    use parameters_constant_mod
    use operations
    implicit none
    type(t_domain), intent(in) :: dm

    real(WP) :: alcai, aci, bci
    

    if (nrank == 0) call Print_debug_start_msg("Building up the interface for the poisson solver ...")
!----------------------------------------------------------------------------------------------------------
    istret = dm%istret
!----------------------------------------------------------------------------------------------------------
    xlx = dm%lxx
    yly = dm%lyt - dm%lyb ! check computational or physical length?
    zlz = dm%lzz
!----------------------------------------------------------------------------------------------------------
    nclx = dm%is_periodic(1)
    ncly = dm%is_periodic(2)
    nclz = dm%is_periodic(3)
!----------------------------------------------------------------------------------------------------------
!   nclx1, ncly1, nclz1 are not used for poisson solver but only for debugging.
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx(1, 1) == IBC_PERIODIC ) then
      nclx1 = 0
    else if (dm%ibcx(1, 1) == IBC_DIRICHLET ) then
      nclx1 = 2
    else
      nclx1 = 1
    end if

    if(dm%ibcy(1, 1) == IBC_PERIODIC ) then
      ncly1 = 0
    else if (dm%ibcy(1, 1) == IBC_DIRICHLET ) then
      ncly1 = 2
    else
      ncly1 = 1
    end if

    if(dm%ibcz(1, 1) == IBC_PERIODIC ) then
      nclz1 = 0
    else if (dm%ibcz(1, 1) == IBC_DIRICHLET ) then
      nclz1 = 2
    else
      nclz1 = 1
    end if
!----------------------------------------------------------------------------------------------------------
    if (nclx) then
      nx = dm%np_geo(1) - 1
      nxm = nx
    else
      nx = dm%np_geo(1)
      nxm = nx - 1
    end if

    if (ncly) then
      ny = dm%np_geo(2) - 1
      nym = ny
    else
      ny = dm%np_geo(2)
      nym = ny - 1
    end if

    if (nclz) then
      nz = dm%np_geo(3) - 1
      nzm = nz
    else
      nz = dm%np_geo(3)
      nzm = nz - 1
    end if
!----------------------------------------------------------------------------------------------------------
    dx = dm%h(1)
    dy = (dm%lyt - dm%lyb) / real(dm%nc(2), WP) !dm%h(2) ! check, computational or physical grid spacing (yes))?
    dz = dm%h(3)
!----------------------------------------------------------------------------------------------------------
    !alpha, beta from geo 
!----------------------------------------------------------------------------------------------------------
    ! if(dm%iAccuracy == IACCU_CD2) then
    !   alcai = ZERO
    !   aci = ONE 
    !   bci = ZERO
    ! else
    !   alcai = NINE / SIXTYTWO
    !   aci = SIXTYTHREE / SIXTYTWO
    !   bci = SEVENTEEN / SIXTYTWO / THREE
    ! end if
    
    ! alcaix6 = alcai    !d1fC2P(3, 1, IBC_PERIODIC)
    ! acix6   = aci / dx !d1rC2P(3, 1, IBC_PERIODIC) / dx
    ! bcix6   = bci / dx !d1rC2P(3, 2, IBC_PERIODIC) / dx

    ! alcaiy6 = alcai    !d1fC2P(3, 1, IBC_PERIODIC)
    ! aciy6   = aci / dy !d1rC2P(3, 1, IBC_PERIODIC) / dy
    ! bciy6   = bci / dy !d1rC2P(3, 2, IBC_PERIODIC) / dy

    ! alcaiz6 = alcai    !d1fC2P(3, 1, IBC_PERIODIC)
    ! aciz6   = aci / dz !d1rC2P(3, 1, IBC_PERIODIC) / dz
    ! bciz6   = bci / dz !d1rC2P(3, 2, IBC_PERIODIC) / dz

    ! only IBC_PERIODIC is necessary, as all non-period data are converted to periodic data.
    if(dm%ibcx(1, 1) == IBC_PERIODIC ) then
        alcaix6 = d1fC2P(3, 1, IBC_PERIODIC)
        acix6   = d1rC2P(3, 1, IBC_PERIODIC) / dx
        bcix6   = d1rC2P(3, 2, IBC_PERIODIC) / dx
    else if (dm%ibcx(1, 1) == IBC_DIRICHLET ) then
        alcaix6 = d1fC2P(3, 1, IBC_DIRICHLET)
        acix6   = d1rC2P(3, 1, IBC_DIRICHLET) / dx
        bcix6   = d1rC2P(3, 2, IBC_DIRICHLET) / dx
    else 
    ! to add and check
    end if

    

    if(dm%ibcy(1, 2) == IBC_PERIODIC ) then
        alcaiy6 = d1fC2P(3, 1, IBC_PERIODIC)
        aciy6   = d1rC2P(3, 1, IBC_PERIODIC) / dy
        bciy6   = d1rC2P(3, 2, IBC_PERIODIC) / dy
    else if (dm%ibcy(1, 2) == IBC_DIRICHLET ) then
        alcaiy6 = d1fC2P(3, 1, IBC_DIRICHLET)
        aciy6   = d1rC2P(3, 1, IBC_DIRICHLET) / dy
        bciy6   = d1rC2P(3, 2, IBC_DIRICHLET) / dy
    else 
    ! to add and check
    end if

    if(dm%ibcz(1, 3) == IBC_PERIODIC ) then
        alcaiz6 = d1fC2P(3, 1, IBC_PERIODIC)
        aciz6   = d1rC2P(3, 1, IBC_PERIODIC) / dz
        bciz6   = d1rC2P(3, 2, IBC_PERIODIC) / dz
    else if (dm%ibcz(1, 3) == IBC_DIRICHLET ) then
        alcaiz6 = d1fC2P(3, 1, IBC_PERIODIC)
        aciz6   = d1rC2P(3, 1, IBC_PERIODIC) / dz
        bciz6   = d1rC2P(3, 2, IBC_PERIODIC) / dz
    else 
    ! to add and check
    end if

#ifdef DEBUG_STEPS
  write(*,*) 'alcaix6, acix6 * dx, bcix6 * dx = ', alcaix6, acix6 * dx, bcix6 * dx
  write(*,*) 'alcaiy6, aciy6 * dy, bciy6 * dy = ', alcaiy6, aciy6 * dy, bciy6 * dy
  write(*,*) 'alcaiz6, aciz6 * dz, bciz6 * dz = ', alcaiz6, aciz6 * dz, bciz6 * dz
#endif
!----------------------------------------------------------------------------------------------------------
!   only classic interpolation, no optimized schemes added here. check paper S. Lele 1992
!   check pros of optimized schemes, to do (see below info from xcompact3d)
!*``ipinter=1``: conventional sixth-order interpolation coefficients as described in `Lele 1992 <https://www.sciencedirect.com/science/article/pii/002199919290324R>`_\
!*``ipinter=2``: optimal sixth-order interpolation coefficients designed to be as close as possible to spectral interpolators.
!*``ipinter=3``: aggressive sixth-order interpolation coefficients designed to add some numerical dissipation at small scales but they could result in spurious oscillations close to a wall.
    ! if(dm%iAccuracy == IACCU_CD2) then
    !   ailcaix6 = ZERO
    !   aicix6 = HALF 
    !   bicix6 = ZERO
    !   cicix6 = ZERO
    !   dicix6 = ZERO
    ! else
    !   ailcaix6 = THREE * ZPONE
    !   aicix6 = ONEPFIVE * HALF
    !   bicix6 = ONE * ZPONE * HALF
    !   cicix6 = ZERO
    !   dicix6 = ZERO
    ! end if

    if(dm%ibcx(1, 1) == IBC_PERIODIC ) then
        ailcaix6 = m1fC2P(3, 1, IBC_PERIODIC)
        aicix6   = m1rC2P(3, 1, IBC_PERIODIC)
        bicix6   = m1rC2P(3, 2, IBC_PERIODIC) 
        cicix6   = zero
        dicix6   = zero
    else if (dm%ibcx(1, 1) == IBC_DIRICHLET ) then
        ailcaix6 = m1fC2P(3, 1, IBC_DIRICHLET)
        aicix6   = m1rC2P(3, 1, IBC_DIRICHLET)
        bicix6   = m1rC2P(3, 2, IBC_DIRICHLET) 
        cicix6   = zero
        dicix6   = zero
    else 
    ! to add and check
    end if

    if(dm%ibcy(1, 2) == IBC_PERIODIC ) then
        ailcaiy6 = m1fC2P(3, 1, IBC_PERIODIC)
        aiciy6   = m1rC2P(3, 1, IBC_PERIODIC)
        biciy6   = m1rC2P(3, 2, IBC_PERIODIC) 
        ciciy6   = zero
        diciy6   = zero
    else if (dm%ibcy(1, 2) == IBC_DIRICHLET ) then
        ailcaiy6 = m1fC2P(3, 1, IBC_DIRICHLET)
        aiciy6   = m1rC2P(3, 1, IBC_DIRICHLET)
        biciy6   = m1rC2P(3, 2, IBC_DIRICHLET) 
        ciciy6   = zero
        diciy6   = zero
    else 
    ! to add and check
    end if

    if(dm%ibcz(1, 3) == IBC_PERIODIC ) then
        ailcaiz6 = m1fC2P(3, 1, IBC_PERIODIC)
        aiciz6   = m1rC2P(3, 1, IBC_PERIODIC)
        biciz6   = m1rC2P(3, 2, IBC_PERIODIC) 
        ciciz6   = zero
        diciz6   = zero
    else if (dm%ibcz(1, 3) == IBC_DIRICHLET ) then
        ailcaiz6 = m1fC2P(3, 1, IBC_DIRICHLET)
        aiciz6   = m1rC2P(3, 1, IBC_DIRICHLET)
        biciz6   = m1rC2P(3, 2, IBC_DIRICHLET) 
        ciciz6   = zero
        diciz6   = zero
    else 
    ! to add and check
    end if

#ifdef DEBUG_STEPS
  write(*,*) 'ailcaix6, aicix6, bicix6 = ', ailcaix6, aicix6, bicix6
  write(*,*) 'ailcaiy6, aiciy6, biciy6 = ', ailcaiy6, aiciy6, biciy6
  write(*,*) 'ailcaiz6, aiciz6, biciz6 = ', ailcaiz6, aiciz6, biciz6
#endif
!----------------------------------------------------------------------------------------------------------

    !module waves
    allocate(zkz(nz/2+1))
    zkz=zero
    allocate(zk2(nz/2+1))
    zk2=zero
    allocate(ezs(nz/2+1))
    ezs=zero

    allocate(yky(ny))
    yky=zero
    allocate(yk2(ny))
    yk2=zero
    allocate(eys(ny))
    eys=zero

    allocate(xkx(nx))
    xkx=zero
    allocate(xk2(nx))
    xk2=zero
    allocate(exs(nx))
    exs=zero

    if (nrank == 0) call Print_debug_end_msg

    return
  end subroutine build_up_poisson_interface

end module

!==========================================================================================================
! below functions and subroutines are from incompact3d.
! please do not change them except "use xxx"
!==========================================================================================================

!##################################################################
! function rl(complexnumber) from incompact3d
!##################################################################
function rl(complexnumber)

  !use param
  use decomp_2d, only: mytype

  implicit none

  real(mytype) :: rl
  complex(mytype) :: complexnumber

  rl = real(complexnumber, kind=mytype)

end function rl
!##################################################################
! function iy(complexnumber) from incompact3d
!##################################################################
function iy(complexnumber)

  !use param
  use decomp_2d, only: mytype

  implicit none

  real(mytype) :: iy
  complex(mytype) :: complexnumber

  iy = aimag(complexnumber)

end function iy
!##################################################################
! function cx(realpart,imaginarypart) from incompact3d
!##################################################################
function cx(realpart,imaginarypart)

  !use param
  use decomp_2d, only: mytype

  implicit none

  complex(mytype) :: cx
  real(mytype) :: realpart, imaginarypart

  cx = cmplx(realpart, imaginarypart, kind=mytype)

end function cx
!==========================================================================================================
!##################################################################
!##################################################################
subroutine inversion5_v1(aaa_in,eee,spI)

  use decomp_2d
  !use decomp_2d_poisson
  !use variables
  !use param
  !use var
  !use mpi
  !use dbg_schemes, only: abs_prec
  use poisson_interface_mod

  implicit none

  ! decomposition object for spectral space
  TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16_mytype
#else
  real(mytype), parameter :: epsilon = 1.e-8_mytype
#endif

  complex(mytype),dimension(spI%yst(1):spI%yen(1),ny/2,spI%yst(3):spI%yen(3),5) :: aaa, aaa_in
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(2):spI%yen(2),spI%yst(3):spI%yen(3)) :: eee
  integer :: i,j,k,m,mi,jc
  integer,dimension(2) :: ja,jb
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

  real(mytype) :: tmp1,tmp2,tmp3,tmp4

  complex(mytype) :: cx
  real(mytype) :: rl, iy
  external cx, rl, iy

  aaa = aaa_in

  do i = 1, 2
     ja(i) = 4 - i
     jb(i) = 5 - i
  enddo
  do m = 1, ny/2 - 2
     do i = 1, 2
        mi = m + i
        do k = spI%yst(3), spI%yen(3)
           do j = spI%yst(1), spI%yen(1)
              if (rl(aaa(j,m,k,3)) /= zero) tmp1 = rl(aaa(j,mi,k,3-i)) / rl(aaa(j,m,k,3))
              if (iy(aaa(j,m,k,3)) /= zero) tmp2 = iy(aaa(j,mi,k,3-i)) / iy(aaa(j,m,k,3))
              sr(j,k)=cx(tmp1,tmp2)
              eee(j,mi,k)=cx(rl(eee(j,mi,k)) - tmp1 * rl(eee(j,m,k)),&
                             iy(eee(j,mi,k)) - tmp2 * iy(eee(j,m,k)))
           enddo
        enddo
        do jc = ja(i), jb(i)
           do k = spI%yst(3), spI%yen(3)
              do j = spI%yst(1), spI%yen(1)
                 aaa(j,mi,k,jc) = cx(rl(aaa(j,mi,k,jc)) - rl(sr(j,k)) * rl(aaa(j,m,k,jc+i)),&
                                     iy(aaa(j,mi,k,jc)) - iy(sr(j,k)) * iy(aaa(j,m,k,jc+i)))
              enddo
           enddo
        enddo
     enddo
  enddo

  do k = spI%yst(3), spI%yen(3)
     do j = spI%yst(1), spI%yen(1)
        if (abs_prec(rl(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp1 = rl(aaa(j,ny/2,k,2)) / rl(aaa(j,ny/2-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs_prec(iy(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp2 = iy(aaa(j,ny/2,k,2)) / iy(aaa(j,ny/2-1,k,3))
        else
           tmp2 = zero
        endif
        sr(j,k) = cx(tmp1,tmp2)
        b1(j,k) = cx(rl(aaa(j,ny/2,k,3)) - tmp1 * rl(aaa(j,ny/2-1,k,4)),&
                     iy(aaa(j,ny/2,k,3)) - tmp2 * iy(aaa(j,ny/2-1,k,4)))

        if (abs_prec(rl(b1(j,k))) > epsilon) then
           tmp1 = rl(sr(j,k)) / rl(b1(j,k))
           tmp3 = rl(eee(j,ny/2,k)) / rl(b1(j,k)) - tmp1 * rl(eee(j,ny/2-1,k))
        else
           tmp1 = zero
           tmp3 = zero
        endif
        if (abs_prec(iy(b1(j,k))) > epsilon) then
           tmp2 = iy(sr(j,k)) / iy(b1(j,k))
           tmp4 = iy(eee(j,ny/2,k)) / iy(b1(j,k)) - tmp2 * iy(eee(j,ny/2-1,k))
        else
           tmp2 = zero
           tmp4 = zero
        endif
        a1(j,k) = cx(tmp1,tmp2)
        eee(j,ny/2,k) = cx(tmp3,tmp4)

        if (abs_prec(rl(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp1 = one / rl(aaa(j,ny/2-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs_prec(iy(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp2 = one / iy(aaa(j,ny/2-1,k,3))
        else
           tmp2 = zero
        endif
        b1(j,k) = cx(tmp1, tmp2)
        a1(j,k) = cx(rl(aaa(j,ny/2-1,k,4)) * rl(b1(j,k)),&
                     iy(aaa(j,ny/2-1,k,4)) * iy(b1(j,k)))
        eee(j,ny/2-1,k) = cx(rl(eee(j,ny/2-1,k)) * rl(b1(j,k)) - rl(a1(j,k)) * rl(eee(j,ny/2,k)),&
                             iy(eee(j,ny/2-1,k)) * iy(b1(j,k)) - iy(a1(j,k)) * iy(eee(j,ny/2,k)))
     enddo
  enddo

  do i = ny/2 - 2, 1, -1
     do k = spI%yst(3), spI%yen(3)
        do j = spI%yst(1), spI%yen(1)
           if (abs_prec(rl(aaa(j,i,k,3))) > epsilon) then
              tmp1 = one / rl(aaa(j,i,k,3))
           else
              tmp1 = zero
           endif
           if (abs_prec(iy(aaa(j,i,k,3))) > epsilon) then
              tmp2 = one/iy(aaa(j,i,k,3))
           else
              tmp2 = zero
           endif
           sr(j,k) = cx(tmp1,tmp2)
           a1(j,k) = cx(rl(aaa(j,i,k,4)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,4)) * iy(sr(j,k)))
           b1(j,k) = cx(rl(aaa(j,i,k,5)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,5)) * iy(sr(j,k)))
           eee(j,i,k) = cx(rl(eee(j,i,k)) * rl(sr(j,k)) - rl(a1(j,k)) * rl(eee(j,i+1,k)) - rl(b1(j,k)) * rl(eee(j,i+2,k)),&
                           iy(eee(j,i,k)) * iy(sr(j,k)) - iy(a1(j,k)) * iy(eee(j,i+1,k)) - iy(b1(j,k)) * iy(eee(j,i+2,k)))
        enddo
     enddo
  enddo

  return

end subroutine inversion5_v1
!##################################################################
!##################################################################
subroutine inversion5_v2(aaa,eee,spI)

  use decomp_2d
  !use decomp_2d_poisson
  !use variables
  !use param
  !use var
  !use MPI
  !use dbg_schemes, only: abs_prec
  use poisson_interface_mod

  implicit none

  ! decomposition object for spectral space
  TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16_mytype
#else
  real(mytype), parameter :: epsilon = 1.e-8_mytype
#endif

  complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3),5) :: aaa
  complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3)) :: eee
  integer :: i,j,k,m,mi,jc
  integer,dimension(2) :: ja,jb
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

  real(mytype) :: tmp1,tmp2,tmp3,tmp4

  complex(mytype) :: cx
  real(mytype) :: rl, iy
  external cx, rl, iy

  do i = 1, 2
     ja(i) = 4 - i
     jb(i) = 5 - i
  enddo
  do m = 1, nym - 2
     do i = 1, 2
        mi = m + i
        do k = spI%yst(3), spI%yen(3)
           do j = spI%yst(1), spI%yen(1)
              if (rl(aaa(j,m,k,3)) /= zero) tmp1 = rl(aaa(j,mi,k,3-i)) / rl(aaa(j,m,k,3))
              if (iy(aaa(j,m,k,3)) /= zero) tmp2 = iy(aaa(j,mi,k,3-i)) / iy(aaa(j,m,k,3))
              sr(j,k) = cx(tmp1, tmp2)
              eee(j,mi,k) = cx(rl(eee(j,mi,k)) - tmp1 * rl(eee(j,m,k)),&
                               iy(eee(j,mi,k)) - tmp2 * iy(eee(j,m,k)))
           enddo
        enddo
        do jc = ja(i), jb(i)
           do k = spI%yst(3), spI%yen(3)
              do j = spI%yst(1), spI%yen(1)
                 aaa(j,mi,k,jc) = cx(rl(aaa(j,mi,k,jc)) - rl(sr(j,k)) * rl(aaa(j,m,k,jc+i)),&
                                     iy(aaa(j,mi,k,jc)) - iy(sr(j,k)) * iy(aaa(j,m,k,jc+i)))
              enddo
           enddo
        enddo
     enddo
  enddo
  do k = spI%yst(3), spI%yen(3)
     do j = spI%yst(1), spI%yen(1)
        if (abs_prec(rl(aaa(j,nym-1,k,3))) > epsilon) then
           tmp1 = rl(aaa(j,nym,k,2)) / rl(aaa(j,nym-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs_prec(iy(aaa(j,nym-1,k,3))) > epsilon) then
           tmp2 = iy(aaa(j,nym,k,2)) / iy(aaa(j,nym-1,k,3))
        else
           tmp2 = zero
        endif
        sr(j,k) = cx(tmp1,tmp2)
        b1(j,k) = cx(rl(aaa(j,nym,k,3)) - tmp1 * rl(aaa(j,nym-1,k,4)),&
                     iy(aaa(j,nym,k,3)) - tmp2 * iy(aaa(j,nym-1,k,4)))
        if (abs_prec(rl(b1(j,k))) > epsilon) then
           tmp1 = rl(sr(j,k)) / rl(b1(j,k))
           tmp3 = rl(eee(j,nym,k)) / rl(b1(j,k)) - tmp1 * rl(eee(j,nym-1,k))
        else
           tmp1 = zero
           tmp3 = zero
        endif
        if (abs_prec(iy(b1(j,k))) > epsilon) then
           tmp2 = iy(sr(j,k)) / iy(b1(j,k))
           tmp4 = iy(eee(j,nym,k)) / iy(b1(j,k)) - tmp2 * iy(eee(j,nym-1,k))
        else
           tmp2 = zero
           tmp4 = zero
        endif
        a1(j,k) = cx(tmp1, tmp2)
        eee(j,nym,k) = cx(tmp3, tmp4)

        if (abs_prec(rl(aaa(j,nym-1,k,3))) > epsilon) then
           tmp1 = one / rl(aaa(j,nym-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs_prec(iy(aaa(j,nym-1,k,3))) > epsilon) then
           tmp2 = one / iy(aaa(j,nym-1,k,3))
        else
           tmp2 = zero
        endif
        b1(j,k) = cx(tmp1,tmp2)
        a1(j,k) = cx(rl(aaa(j,nym-1,k,4)) * rl(b1(j,k)),&
                     iy(aaa(j,nym-1,k,4)) * iy(b1(j,k)))
        eee(j,nym-1,k) = cx(rl(eee(j,nym-1,k)) * rl(b1(j,k)) - rl(a1(j,k)) * rl(eee(j,nym,k)),&
                            iy(eee(j,nym-1,k)) * iy(b1(j,k)) - iy(a1(j,k)) * iy(eee(j,nym,k)))
     enddo
  enddo

  do i = nym - 2, 1, -1
     do k = spI%yst(3), spI%yen(3)
        do j = spI%yst(1), spI%yen(1)
           if (abs_prec(rl(aaa(j,i,k,3))) > epsilon) then
              tmp1 = one / rl(aaa(j,i,k,3))
           else
              tmp1 = zero
           endif
           if (abs_prec(iy(aaa(j,i,k,3))) > epsilon) then
              tmp2 = one / iy(aaa(j,i,k,3))
           else
              tmp2 = zero
           endif
           sr(j,k) = cx(tmp1,tmp2)
           a1(j,k) = cx(rl(aaa(j,i,k,4)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,4)) * iy(sr(j,k)))
           b1(j,k) = cx(rl(aaa(j,i,k,5)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,5)) * iy(sr(j,k)))
           eee(j,i,k) = cx(rl(eee(j,i,k)) * rl(sr(j,k)) - rl(a1(j,k)) * rl(eee(j,i+1,k)) -rl(b1(j,k)) * rl(eee(j,i+2,k)),&
                           iy(eee(j,i,k)) * iy(sr(j,k)) - iy(a1(j,k)) * iy(eee(j,i+1,k)) -iy(b1(j,k)) * iy(eee(j,i+2,k)))
        enddo
     enddo
  enddo

  return

end subroutine inversion5_v2



module decomp_extended_mod
  use parameters_constant_mod
  implicit none


  public :: ypencil_index_lgl2ggl
  public :: zpencil_index_llg2ggg
  public :: zpencil_index_ggg2llg

  contains
!==========================================================================================================
  subroutine ypencil_index_lgl2ggl(vin, vou, dtmp)
    use decomp_2d
    implicit none

    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), dimension(dtmp%ysz(1),               dtmp%ysz(2), dtmp%ysz(3)), intent(in)  :: vin
    real(WP), dimension(dtmp%yst(1) : dtmp%yen(2), dtmp%ysz(2), dtmp%zsz(3)), intent(out) :: vou

    integer :: i, j, k, ii
    vou = ZERO
    do k = 1, dtmp%ysz(3)
      do j = 1, dtmp%ysz(2)
        do i = 1, dtmp%ysz(1)
          ii = dtmp%yst(1) + i - 1
          vou(ii, j, k) = vin(i, j, k)
        end do
      end do
    end do
    return
  end subroutine 
!==========================================================================================================
  subroutine zpencil_index_llg2ggg(vin, vou, dtmp)
    use decomp_2d
    implicit none

    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), dimension(dtmp%zsz(1),               dtmp%zsz(2),               dtmp%zsz(3)), intent(in)  :: vin
    real(WP), dimension(dtmp%zst(1) : dtmp%zen(1), dtmp%zst(2) : dtmp%zen(2), dtmp%zsz(3)), intent(out) :: vou

    integer :: i, j, k, jj, ii

    vou = ZERO
    do k = 1, dtmp%zsz(3)
      do j = 1, dtmp%zsz(2)
        jj = dtmp%zst(2) + j - 1
        do i = 1, dtmp%zsz(1)
          ii = dtmp%zst(1) + i - 1
          vou(ii, jj, k) = vin(i, j, k)
        end do
      end do
    end do
    return
  end subroutine
!==========================================================================================================  
  subroutine zpencil_index_ggg2llg(vin, vou, dtmp)
    use decomp_2d
    implicit none

    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), dimension(dtmp%zst(1) : dtmp%zen(1), dtmp%zst(2) : dtmp%zen(2), dtmp%zsz(3)), intent(in)   :: vin
    real(WP), dimension(dtmp%zsz(1),               dtmp%zsz(2),               dtmp%zsz(3)), intent(out)  :: vou
    

    integer :: i, j, k, jj, ii
!write(*,*) 'vin', nrank, size(vin, 1), size(vin, 2),size(vin, 3)
!write(*,*) 'vou', nrank, size(vou, 1), size(vou, 2),size(vou, 3)

    vou = ZERO
    do k = 1, dtmp%zsz(3)
      do j = 1, dtmp%zsz(2)
        jj = dtmp%zst(2) + j - 1
        do i = 1, dtmp%zsz(1)
          ii = dtmp%zst(1) + i - 1
          vou(i, j, k) = vin(ii, jj, k)
        end do
      end do
    end do
    return
  end subroutine

  

end module 