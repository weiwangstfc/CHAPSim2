subroutine Print_error_msg_mpi(errorcode, msg)
  use mpi_mod
  implicit none
  integer, intent(IN) :: errorcode
  character(len=*), intent(IN) :: msg
  
  if (myid==0) then
    write(*,*) 'CHAPSim ERROR - errorcode: ', errorcode
    write(*,*) 'ERROR: ' // msg
  end if
  call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)

  return
end subroutine Print_error_msg_mpi

subroutine Print_error_msg(msg)
  use mpi_mod
  implicit none
  character(len=*), intent(IN) :: msg
  
  write(*,*) 'ERROR: ' // msg

  STOP

  return
end subroutine Print_error_msg


subroutine Print_warning_msg(msg)
  use mpi_mod
  implicit none
  character(len=*), intent(IN) :: msg
  
  write(*,*) 'WARNNING: ' // msg

  return
end subroutine Print_warning_msg

subroutine Display_vtk_slice(d, str, varnm, var1, var2, var3, var4)
  use udf_type_mod
  type(domain_t), intent( in ) :: d
  character( len = *), intent( in ) :: str
  character( len = *), intent( in ) :: varnm
  real(WP), dimension(:, :, :), intent(in):: var1, var2, var3, var4

  integer :: output_unit
  integer :: i, j, k
  integer :: nd1, nd2, nd3
  integer :: nc1, nc2, nc3
  real(WP) :: x, y, z
  character( len = 128) :: filename

  nd1 = d%np(1)
  nd2 = d%np(2)
  nd3 = d%np(3)

  nc1 = d%nc(1)
  nc2 = d%nc(2)
  nc3 = d%nc(3)

  if ( trim( str ) == 'xy' ) then
    nd3 = 1
    nc3 = 1
  end if
  if ( trim( str ) == 'yz' ) then
    nd1 = 1
    nc1 = 1
  end if
  if ( trim( str ) == 'zx' ) then
    nd2 = 1
    nc2 = 1
  end if

  filename = 'display_'//varnm//'_'//str//'.vtk'
  open(newunit = output_unit, file = trim(filename), action = "write", status = "replace")
  write(output_unit, '(A)') '# vtk DataFile Version 2.0'
  write(output_unit, '(A)') str//'_mesh'
  write(output_unit, '(A)') 'ASCII'
  write(output_unit, '(A)') 'DATASET STRUCTURED_GRID'
  write(output_unit, '(A, 3I10.1)') 'DIMENSIONS', nd1, nd2, nd3
  write(output_unit, '(A, I10.1, X, A)') 'POINTS', nd1 * nd2 * nd3, 'float'
  do k = 1, nd3
    z = d%dz * real(k - 1, WP)
    do j = 1, nd2
      y = d%yp(j)
      do i = 1, nd1
        x = d%dx * real(i - 1, WP)
        write(output_unit, *) x, y, z
      end do
    end do
  end do

  ! data convert to points data...




  ! write out points data...
  write(output_unit, '(A, 1I10.1)') 'CELL_DATA', nc1 * nc2 * nc3
  write(output_unit, '(A, 1I10.1)') 'SCALARS '//varnm(1:1)//' float',  nc1 * nc2 * nc3
  write(output_unit, '(A)') 'LOOKUP_TABLE default'

  do k = 1, nc3
    do j = 1, nc2
      do i = 1, nc1
        write(output_unit, *) var1(i, j, k)
      end do
    end do
  end do

  write(output_unit, '(A, 1I10.1)') 'SCALARS '//varnm(2:2)//' float',  nc1 * nc2 * nc3
  write(output_unit, '(A)') 'LOOKUP_TABLE default'

  do k = 1, nc3
    do j = 1, nc2
      do i = 1, nc1
        write(output_unit, *) var2(i, j, k)
      end do
    end do
  end do

  write(output_unit, '(A, 1I10.1)') 'SCALARS '//varnm(3:3)//' float',  nc1 * nc2 * nc3
  write(output_unit, '(A)') 'LOOKUP_TABLE default'

  do k = 1, nc3
    do j = 1, nc2
      do i = 1, nc1
        write(output_unit, *) var3(i, j, k)
      end do
    end do
  end do


    write(output_unit, '(A, 1I10.1)') 'SCALARS '//varnm(4:4)//' float',  nc1 * nc2 * nc3
    write(output_unit, '(A)') 'LOOKUP_TABLE default'

    do k = 1, nc3
      do j = 1, nc2
        do i = 1, nc1
          write(output_unit, *) var4(i, j, k)
        end do
      end do
    end do


  close(output_unit)

  return

end subroutine Display_vtk_slice



module random_number_generation_mod
  use precision_mod
  implicit none
  private
  public :: Initialize_random_number
  public :: Generate_rvec_random
  public :: Generate_r_random

contains
  subroutine Initialize_random_number ( seed )
    !*******************************************************************************
    !
    !! random_initialize initializes the FORTRAN 90 random number seed.
    !
    !
    !  Discussion:
    !
    !    If you don't initialize the random number generator, its behavior
    !    is not specified.  If you initialize it simply by:
    !
    !      CALL random_seed
    !
    !    its behavior is not specified.  On the DEC ALPHA, If that's all you
    !    do, the same random number sequence is returned.  In order to actually
    !    try to scramble up the random number generator a bit, this routine
    !    goes through the tedious process of getting the size of the random
    !    number seed, making up values based on the current time, and setting
    !    the random number seed.
    !
    !  Modified:
    !
    !    19 December 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  parameters:
    !
    !    Input/output, INTEGER(4) seed.
    !    IF seed is zero on input, THEN you're asking this routine to come up
    !    with a seed value, whICh is RETURNed as output.
    !    IF seed is nonzero on input, THEN you're asking this routine to
    !    USE the input value of seed to initialize the random number generator,
    !    and seed is not changed on output.
    !
    implicit none
    !
    integer(4) :: count
    integer(4) :: count_max
    integer(4) :: count_rate
    logical, parameter :: debug = .false.
    integer(4) :: i
    integer(4) :: seed
    integer(4), allocatable :: seed_vector(:)
    integer(4) :: seed_size
    real(wp) :: t
    !
    !  Initialize the random number seed.
    !
    call random_seed
    !
    !  determine the size of the random number seed.
    !
    call random_seed ( size = seed_size )
    !
    !  allocate a seed of the right size.
    !
    allocate ( seed_vector(seed_size) ); seed_vector = 0

    if ( seed /= 0 ) then

        if ( debug ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'random_initialize'
            write ( *, '(a, i20)' ) '  initialize random_number, user seed = ', seed
        end if

    else

        call system_clock ( count, count_rate, count_max )

        seed = count

        if ( debug ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'random_initialize'
            write ( *, '(a, i20)' ) '  initialize random_number, arbitrary seed = ', &
            seed
        end if

    end if
    !
    !  now set the seed.
    !
    seed_vector(1:seed_size) = seed

    call random_seed ( put = seed_vector(1:seed_size) )
    !
    !  free up the seed space.
    !
    deallocate ( seed_vector )
    !
    !  call the random number routine a bunch of times.
    !random_initialize
    do i = 1, 100
        call random_number ( harvest = t )
    end do

    return
  end subroutine Initialize_random_number

  !**********************************************************************************************************************************
  subroutine Generate_rvec_random ( alo, ahi, n, a )
    !
    !*******************************************************************************
    !
    !! RVEC_random RETURNs a random REAL(WP) vector in a given range.
    !
    !
    !  ModIFied:
    !
    !    04 FebruARy 2001
    !
    !  Author:
    !
    !    John BurkARdt
    !
    !  parameters:
    !
    !    Input, REAL(WP) ALO, AHI, the range allowed for the entries.
    !
    !    Input, INTEGER(4) N, the number of entries in the vector.
    !
    !    Output, REAL(WP) A(N), the vector of randomly chosen values.
    !
    implicit none
    !
    integer(4) n
    !
    real(wp) a(n)
    real(wp) ahi
    real(wp) alo
    integer(4) i
    !
    do i = 1, n
        call Generate_r_random ( alo, ahi, a(i) )
    end do

    return
  end subroutine Generate_rvec_random

!**********************************************************************************************************************************
  subroutine Generate_r_random ( rlo, rhi, r )
    !
    !*******************************************************************************
    !
    !! R_random RETURNs a random REAL(WP) in a given range.
    !
    !
    !  ModIFied:
    !
    !    06 April 2001
    !
    !  Author:
    !
    !    John BurkARdt
    !
    !  parameters:
    !
    !    Input, REAL(WP) RLO, RHI, the minimum and maximum values.
    !
    !    Output, REAL(WP) R, the randomly chosen value.
    !
    implicit none
    !
    real(wp) :: r
    real(wp) :: rhi
    real(wp) :: rlo
    real(wp) :: t
    !
    !  pick t, a random number in (0, 1).
    !
    call random_number ( harvest = t )
    !
    !  set r in ( rlo, rhi ).
    !
    r = ( 1.0e+00 - t ) * rlo + t * rhi

    return
  end subroutine Generate_r_random
end module random_number_generation_mod

