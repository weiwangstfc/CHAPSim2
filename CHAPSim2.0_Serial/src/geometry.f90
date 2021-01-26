module geometry_mod
  use precision_mod
  implicit none

  type domain_t
    integer :: bcx(2)
    integer :: bcy(2)
    integer :: bcz(2)
    logical :: is_periodic(3)
    integer :: np(3)
    integer :: nc(3)
    real(wp) :: dx
    real(wp) :: dz
    real(wp) :: dx2
    real(wp) :: dz2
    real(wp) :: dxi
    real(wp) :: dzi
    real(wp), allocatable :: yp(:)
    real(wp), allocatable :: yc(:)
  end type domain_t

  type(domain_t), save :: domain

  ! node location, mapping 
  real(wp), allocatable :: mpd1c1(:) ! first coefficient in first deriviation. 1/h'
  real(wp), allocatable :: mpd2c1(:) ! first coefficient in second deriviation 1/h'^2
  real(wp), allocatable :: mpd2c2(:) ! second coefficient in second deriviation -h"/h'^3
  ! cell center location, mapping
  real(wp), allocatable :: mcd1c1(:) ! first coefficient in first deriviation. 1/h'
  real(wp), allocatable :: mcd2c1(:) ! first coefficient in second deriviation 1/h'^2
  real(wp), allocatable :: mcd2c2(:) ! second coefficient in second deriviation -h"/h'^3

  !private
  private :: Buildup_grid_mapping_1D
  public :: Display_vtk_mesh
  public :: Initialize_geometry_variables
  
contains

  subroutine Display_vtk_mesh(d, str)
    type(domain_t), intent( in ) :: d
    character( len = *), intent( in ) :: str

    integer :: output_unit
    integer :: i, j, k
    integer :: n1, n2, n3
    real(WP) :: x, y, z

    n1 = d%np(1)
    n2 = d%np(2)
    n3 = d%np(3)

    if ( trim( str ) == 'xy' ) then
      n3 = 1
    end if
    if ( trim( str ) == 'yz' ) then
      n1 = 1
    end if
    if ( trim( str ) == 'zx' ) then
      n2 = 1
    end if

    open(newunit = output_unit, file = 'mesh_'//str//'.vtk', action = "write", status = "replace")
    write(output_unit, '(A)') '# vtk DataFile Version 2.0'
    write(output_unit, '(A)') str//'_mesh'
    write(output_unit, '(A)') 'ASCII'
    write(output_unit, '(A)') 'DATASET STRUCTURED_GRID'
    write(output_unit, '(A, 3I10.1)') 'DIMENSIONS', n1, n2, n3
    write(output_unit, '(A, I10.1, X, A)') 'POINTS', n1 * n2 * n3, 'float'
    do k = 1, n3
      z = d%dz * real(k - 1, WP)
      do j = 1, n2
        y = d%yp(j)
        do i = 1, n1
          x = d%dx * real(i - 1, WP)
          write(output_unit, *) x, y, z
        end do
      end do
    end do
    close(output_unit)

  end subroutine Display_vtk_mesh

  subroutine Initialize_geometry_variables ()
    use mpi_mod
    use input_general_mod
    use math_mod
    use parameters_constant_mod, only : ONE, HALF, ZERO, MAXP, MINP, TRUNCERR

    ! Build up domain info
    domain%bcx(:) = ifbcx(:)
    domain%bcy(:) = ifbcy(:)
    domain%bcz(:) = ifbcz(:)

    domain%is_periodic(:) = is_periodic(:)
    
    domain%nc(1) = ncx
    domain%nc(2) = ncy
    domain%nc(3) = ncz

    ! to note: geometric node number is always cell numbe + 1,
    ! different from the flow node number
    domain%np(1) = ncx + 1
    domain%np(2) = ncy + 1
    domain%np(3) = ncz + 1

    domain%dx = lxx / real(domain%nc(1), WP)
    domain%dz = lzz / real(domain%nc(3), WP)

    domain%dx2 = domain%dx * domain%dx
    domain%dz2 = domain%dz * domain%dz

    domain%dxi = ONE / domain%dx
    domain%dzi = ONE / domain%dz

    ! allocate  variables for mapping physical domain to computational domain
    allocate ( domain%yp( domain%np(2) ) ); domain%yp(:) = ZERO
    allocate ( domain%yc( domain%nc(2) ) ); domain%yc(:) = ZERO

    allocate ( mpd1c1( domain%np(2) ) ); mpd1c1(:) =  ONE
    allocate ( mpd2c1( domain%np(2) ) ); mpd2c1(:) =  ONE
    allocate ( mpd2c2( domain%np(2) ) ); mpd2c1(:) =  ONE
    
    allocate ( mcd1c1( domain%nc(2) ) ); mcd1c1(:) =  ONE
    allocate ( mcd2c1( domain%nc(2) ) ); mcd2c1(:) =  ONE
    allocate ( mcd2c2( domain%nc(2) ) ); mcd2c1(:) =  ONE
    
    call Buildup_grid_mapping_1D ('nd', domain%np(2), domain%yp(:), mpd1c1(:), mpd2c1(:), mpd2c1(:))
    call Buildup_grid_mapping_1D ('cl', domain%nc(2), domain%yc(:), mcd1c1(:), mcd2c1(:), mcd2c1(:))

    ! to print geometry/mesh domain
    if(myid == 0) then ! test
      call Display_vtk_mesh ( domain, 'xy' )
      call Display_vtk_mesh ( domain, 'yz' )
      call Display_vtk_mesh ( domain, 'zx' )
      call Display_vtk_mesh ( domain, '3D' )
    end if

  end subroutine  Initialize_geometry_variables


  subroutine Buildup_grid_mapping_1D (str, n, y, mp1, mp2, mp3)
    use math_mod
    use input_general_mod
    use parameters_constant_mod

    character(len = *), intent(in) :: str
    integer(4), intent(in) :: n
    real(WP), intent( out ) :: y(n), mp1(n), mp2(n), mp3(n)

    integer(4) :: j
    real(WP) :: eta_shift
    real(WP) :: eta_delta
    real(WP) :: alpha, beta, gamma, delta, cc, dd, ee, st1, st2, mm
    real(WP), dimension(n) :: eta

    eta_shift = ZERO
    eta_delta = ONE
    if ( trim( str ) == 'nd' ) then
      eta_shift = ZERO
      eta_delta = ONE / real( n - 1, WP )
    else if ( trim( str ) == 'cl' ) then
      eta_shift = ONE / ( real(n, WP) ) * HALF
      eta_delta = ONE / real( n, WP )
    else 
      call Print_error_msg(101, 'Grid stretching location not defined.')
    end if

    ! to build up the computational domain \eta \in [0, 1] uniform mesh
    eta(1) = ZERO + eta_shift

    do j = 2, n
      eta(j) = eta(1) + real(j - 1, WP) * eta_delta
    end do

    ! to build up the physical domain y stretching grids based on Eq(53) of Leizet2009JCP
    ! and to build up the derivates based on Eq(53) and (47) in Leizet2009JCP
    gamma = ONE
    delta = ZERO
    if (istret == ISTRET_NO) then
      y(:) = eta(:)
      y(:) = y(:) * (lyt - lyb) + lyb
      mp1(:) = ONE
      mp2(:) = ONE
      mp3(:) = ONE
      return
    else if (istret == ISTRET_CENTRE) then
      gamma = ONE
      delta = ZERO
    else if (istret == ISTRET_2SIDES) then
      gamma = ONE
      delta = HALF
    else if (istret == ISTRET_BOTTOM) then
      gamma = HALF
      delta = HALF
    else if (istret == ISTRET_TOP) then
      gamma = HALF
      delta = ZERO
    else
      call Print_error_msg('Grid stretching flag is not valid.')
    end if

    beta = rstret
    alpha =  ( -ONE + sqrt_wp( ONE + FOUR * PI * PI * beta * beta ) ) / beta * HALF

    cc = sqrt_wp( alpha * beta + ONE ) / sqrt_wp( beta )
    dd = cc / sqrt_wp( alpha )
    ee = cc * sqrt_wp( alpha )

    st1 = (ONE   - TWO * delta) / gamma * HALF
    st2 = (THREE - TWO * delta) / gamma * HALF

    do j = 1, n
      mm = PI * (gamma * eta(j) + delta)

      ! y \in [0, 1]
      y(j) = atan_wp ( dd * tan_wp( mm ) ) - &
             atan_wp ( dd * tan_wp( PI * delta) ) + &
             PI * ( heaviside_step( eta(j) - st1 ) + heaviside_step( eta(j) - st2 ) )
      y(j) = ONE / (gamma * ee) * y(j)
      ! y \in [lyb, lyt]
      y(j) = y(j) * (lyt - lyb) + lyb

      ! 1/h'
      mp1(j) = (alpha / PI + sin_wp(mm) * sin_wp(mm) / PI / beta)  / (lyt - lyb)

      ! (1/h')^2
      mp2(j) = mp1(j) * mp1(j)

      ! -h"/(h'^3) = 1/h' * [ d(1/h') / d\eta]
      mp3(j) = gamma / (lyt - lyb) / beta * sin_wp(TWO * mm) * mp1(j)

    end do

    return
  end subroutine Buildup_grid_mapping_1D


end module geometry_mod

