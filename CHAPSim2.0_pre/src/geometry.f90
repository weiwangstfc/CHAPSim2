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
  end type domain_t

  type cell_t
    real(WP) :: x
    real(WP) :: y
    real(WP) :: z
  end type cell_t

  type node_t
    real(WP) :: x
    real(WP) :: y
    real(WP) :: z
  end type node_t

  type(domain_t), save :: domain
  type(node_t), save, allocatable, dimension(:, :, :) :: node
  type(cell_t), save, allocatable, dimension(:, :, :) :: cell

  ! to define continuous variable for efficiency
  real(WP), save, allocatable, dimension(:) :: yp
  real(WP), save, allocatable, dimension(:) :: yc
  real(WP), save, allocatable, dimension(:) :: dycc
  real(WP), save, allocatable, dimension(:) :: dypp
  real(WP), save, allocatable, dimension(:) :: rpi
  real(WP), save, allocatable, dimension(:) :: rci

  real(WP), save, allocatable, dimension(:) :: fyc2p0
  real(WP), save, allocatable, dimension(:) :: fyc2p1
  real(WP) :: fxc2p, fzc2p
  real(WP) :: fxp2c, fzp2c, fyp2c


  private
  public :: Display_vtk_mesh
  public :: Initialize_geometry_variables
  
contains

  subroutine Display_vtk_mesh(nd, str)
    type(node_t), intent( in ) :: nd(:, :, :)
    character( len = *), intent( in ) :: str

    integer :: output_unit
    integer :: i, j, k
    integer :: n1, n2, n3

    n1 = size( nd(:, 1, 1) )
    n2 = size( nd(1, :, 1) )
    n3 = size( nd(1, 1, :) )

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
      do j = 1, n2
        do i = 1, n1
          write(output_unit, *) nd(i, j, k)%x, nd(i, j, k)%y, nd(i, j, k)%z 
        end do
      end do
    end do
    close(output_unit)

  end subroutine Display_vtk_mesh

  subroutine Initialize_geometry_variables ()
    use mpi_mod
    use input_general_mod
    use math_mod
    use chapsim_abort_mod
    use parameters_constant_mod, only : ONE, HALF, ZERO, MAXP, MINP, TRUNCERR
    integer :: i, j, k
    real(WP) :: s, yy, c1, c2, c3, c4

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

    ! to build up node and cell info
    allocate ( node (domain%np(1), domain%np(2), domain%np(3)) )
    allocate ( cell (domain%nc(1), domain%nc(2), domain%nc(3)) )
    ! to initialize 
    do k = 1, domain%np(3)
      do j = 1, domain%np(2)
        do i = 1, domain%np(1)
          node(i, j, k)%x = ZERO
          node(i, j, k)%y = ZERO
          node(i, j, k)%z = ZERO
        end do
      end do
    end do
    
    do k = 1, domain%nc(3)
      do j = 1, domain%nc(2)
        do i = 1, domain%nc(1)
          cell(i, j, k)%x = ZERO
          cell(i, j, k)%y = ZERO
          cell(i, j, k)%z = ZERO
        end do
      end do
    end do
    
    ! x 
    do i = 1, domain%np(1)
      node(i, :, :)%x = real( (i - 1), WP ) * domain%dx
    end do

    do i = 1, domain%nc(1)
      cell(i, :, :)%x = real( (i - 1), WP ) * domain%dx + domain%dx * HALF
    end do

    ! z
    do k = 1, domain%np(3)
      node(:, :, k)%z = real( (k - 1), WP ) * domain%dz
    end do

    do k = 1, domain%nc(3)
      cell(:, :, k)%z = real( (k - 1), WP ) * domain%dz + domain%dz * HALF
    end do

    ! y
    if (istret == ISTRET_SIDES) then
      c1 = rstret * HALF
      c2 = tanh_wp (c1)
      c3 = ONE
      c4 = HALF
    else if (istret == ISTRET_BOTTOM) then
      c1 = rstret * ONE
      c2 = tanh_wp (c1)
      c3 = ONE
      c4 = ONE
    else if (istret == ISTRET_TOP) then
      c1 = rstret * ZERO
      c2 = tanh_wp (rstret)
      c3 = ZERO
      c4 = ONE
    else
      c1 = rstret * HALF
      c2 = tanh_wp (c1)
      c3 = ONE
      c4 = HALF
    end if

    do j = 1, domain%np(2)
      yy = real ((j - 1), WP) / real ( (npy - 1), WP)
      if(istret == ISTRET_NO) then
        s = yy
      else 
        s = (tanh_wp( (rstret * yy) - c1 ) / c2 + c3) * c4
      end if
      node(:, j, :)%y = s * (lyt - lyb) + lyb
    end do

    do j = 1, domain%nc(2)
      cell(:, j, :)%y = ( node(1, j, 1)%y + node(1, j + 1, 1)%y ) * HALF
    end do

    ! to print geometry/mesh domain
    if(myid == 0) then ! test
      call Display_vtk_mesh ( node(:, :, :), 'xy' )
      call Display_vtk_mesh ( node(:, :, :), 'yz' )
      call Display_vtk_mesh ( node(:, :, :), 'zx' )
      call Display_vtk_mesh ( node(:, :, :), '3D' )
    end if

    ! to build up the variables to be used in governing equations
    
    ! yp
    allocate ( yp( domain%np(2) ) ); yp = ZERO
    do j = 1, domain%np(2)
      yp(j) = node(1, j, 1)%y
    end do

    ! yc
    allocate ( yc( domain%nc(2) ) ); yc = ZERO
    do j = 1, domain%nc(2)
      yc(j) = cell(1, j, 1)%y
    end do

    ! rci
    allocate ( rpi( domain%np(2) ) ); rpi = ONE
    allocate ( rci( domain%nc(2) ) ); rci = ONE
    rci(:) = ONE
    rpi(:) = ONE
    if(icoordinate == ICYLINDRICAL) then
      do j = 1, domain%np(2)
        if ( abs_wp( yp(j) ) < MINP ) then
          rpi(j) = MAXP
        else
          rpi(j) = ONE / yp(j)
        end if
      end do 
      do j = 1, domain%nc(2)
        rci(j) = ONE / yc(j)
      end do
    end if

    ! dy
    allocate ( dypp( domain%nc(2) ) ); dypp = ZERO
    allocate ( dycc( domain%np(2) ) ); dycc = ZERO
    ! dy from point to point
    do j = 1, domain%nc(2)
      dypp(j) = yp(j + 1) - yp(j)
    end do
    ! dy from cell centre to cell cetre
    do j = 2, domain%np(2) - 1
      dycc(j) = ( yp(j + 1) - yp(j - 1) ) * HALF
    end do

    dycc(1) = ( yp(2) - yp(1) ) * HALF
    dycc(domain%np(2)) = ( yp( domain%np(2) ) - yc( domain%np(2) - 1) ) * HALF
    
    if ( domain%is_periodic(2) ) then
      ! for example, TGV
      dycc(1) = abs_wp( yc( 1            ) - yp( 1          ) ) + &
                abs_wp( yp( domain%np(2) ) - yc( domain%np(2) - 1 ) )
      dycc(domain%np(2)) = dycc(1)
    end if

    if (domain%bcy(1) == IBC_INTERIOR) then
      ! for example, pipe
      dycc(1) = yp(2) - yp(1)
    end if

    ! to calculate interpolation coefficients
    fxp2c = HALF
    fyp2c = HALF
    fzp2c = HALF

    fxc2p = HALF
    fzc2p = HALF

    allocate ( fyc2p0( domain%np(2) ) ); fyc2p0 = ZERO
    allocate ( fyc2p1( domain%np(2) ) ); fyc2p1 = ZERO
    do j = 2, domain%np(2) - 1
      fyc2p0( j ) = ( yc( j ) - yp( j    ) ) / ( yc( j ) - yc( j - 1 ) )
      fyc2p1( j ) = ( yp( j ) - yc( j - 1) ) / ( yc( j ) - yc( j - 1 ) )
    end do

    fyc2p0( 1            ) = ZERO
    fyc2p1( 1            ) = ONE
    fyc2p0( domain%np(2) ) = ONE
    fyc2p1( domain%np(2) ) = ZERO

    if ( domain%is_periodic(2) ) then
      fyc2p0( 1 ) = ( yc( 1            ) - yp( 1          ) ) / &
                    ( yc( 1 ) - yp( 1 ) + yp( domain%np(2) ) - yc( domain%nc(2) ) )
      fyc2p1( 1 ) = ( yp( domain%np(2) ) - yc( domain%nc(2) ) ) /  &
                    ( yc( 1 ) - yp( 1 ) + yp( domain%np(2) ) - yc( domain%nc(2) ) )
    end if

    if (domain%bcy(1) == IBC_INTERIOR) then
      fyc2p0( 1 ) = HALF
      fyc2p0( 2 ) = HALF
    end if

    ! to validate the y distance
    yy = 0
    do j = 1, domain%nc(2)
      yy = yy + dypp(j)
    end do
    if ( abs_wp(yy - (lyt - lyb)) > TRUNCERR ) &
    call Print_error_msg(1, "Error in meshing.")

  end subroutine  Initialize_geometry_variables

end module geometry_mod

