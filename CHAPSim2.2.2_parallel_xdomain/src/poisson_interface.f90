module poisson_interface_mod
  use decomp_2d
  use mpi_mod
  use parameters_constant_mod, only: zero, , half, one, onepfive, two, twopfive, &
                                     three, pi, threepfive, four, twopi, cx_one_one
  use math_mod, only: cos_prec, abs_prec
  implicit none

  integer :: istret
!-------------------------------------------------------------------------------
  logical :: nclx
  logical :: ncly
  logical :: nclz
!-------------------------------------------------------------------------------  
  integer :: nx
  integer :: ny
  integer :: nz
!-------------------------------------------------------------------------------
  integer :: nxm
  integer :: nym
  integer :: nzm
!-------------------------------------------------------------------------------
  integer :: nclx1 
  integer :: ncly1 
  integer :: nclz1
!-------------------------------------------------------------------------------
  real(mytype) :: xlx
  real(mytype) :: yly
  real(mytype) :: zlz

!-------------------------------------------------------------------------------
  real(mytype) :: alpha
  real(mytype) :: beta
!-------------------------------------------------------------------------------
  real(mytype) :: acix6
  real(mytype) :: aciy6
  real(mytype) :: aciz6

  real(mytype) :: alcaix6 
  real(mytype) :: alcaiy6 
  real(mytype) :: alcaiz6

  real(mytype) :: aicix6
  real(mytype) :: aiciy6
  real(mytype) :: aiciz6

  real(mytype) :: ailcaix6 
  real(mytype) :: ailcaiy6 
  real(mytype) :: ailcaiz6
!-------------------------------------------------------------------------------
  real(mytype) :: bcix6
  real(mytype) :: bciy6
  real(mytype) :: bciz6

  real(mytype) :: bicix6
  real(mytype) :: biciy6
  real(mytype) :: biciz6
!-------------------------------------------------------------------------------
  real(mytype) :: cicix6
  real(mytype) :: ciciy6
  real(mytype) :: ciciz6
!-------------------------------------------------------------------------------
  real(mytype) :: dicix6
  real(mytype) :: diciy6
  real(mytype) :: diciz6
!-------------------------------------------------------------------------------
  real(mytype) :: dx
  real(mytype) :: dy
  real(mytype) :: dz


!-------------------------------------------------------------------------------
  !module waves
  complex(mytype),allocatable,dimension(:) :: zkz,zk2,ezs
  complex(mytype),allocatable,dimension(:) :: yky,yk2,eys
  complex(mytype),allocatable,dimension(:) :: xkx,xk2,exs

  public :: build_up_poisson_interface
contains


  subroutine build_up_poisson_interface(dm)
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    type(t_domain), intent(in) :: dm
    

    nclx = dm%is_periodic(1)
    ncly = dm%is_periodic(1)
    nclz = dm%is_periodic(1)

    nx = dm%np_geo(1)
    ny = dm%np_geo(2)
    nz = dm%np_geo(3)

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

    return
  end subroutine build_up_poisson_interface


end module