!===============================================================================
  subroutine Print_error_msg(msg)
    implicit none
    character(len=*), intent(IN) :: msg
    !integer(4), intent(in) :: myid
    
    write(*,*) 'ERROR: ' // msg

    write(*,*) 'Code is terminated in processor = '!, myid
    STOP

    return
  end subroutine Print_error_msg
!===============================================================================
  subroutine Print_warning_msg(msg)
    implicit none
    character(len=*), intent(IN) :: msg
    !integer(4), intent(in) :: myid
    
    write(*,*) 'WARNNING: ' // msg // ' in processor = '!, myid

    return
  end subroutine Print_warning_msg
  !===============================================================================
  subroutine Print_debug_start_msg(msg)
    !use mpi_mod
    implicit none
    character(len=*), intent(IN) :: msg

    !if(myid /= 0) return
    write(*,*) "==============================================================================="
    write(*,*) msg

    return
  end subroutine Print_debug_start_msg
!===============================================================================
  subroutine Print_debug_mid_msg(msg)
    !use mpi_mod
    implicit none
    character(len=*), intent(IN) :: msg
    !if(myid /= 0) return
    write(*,*) msg
    return
  end subroutine Print_debug_mid_msg
!===============================================================================
  subroutine Print_debug_end_msg
    !use mpi_mod
    implicit none
    !if(myid /= 0) return
    write(*,*) "... done."
    return
  end subroutine Print_debug_end_msg

!===============================================================================
  subroutine Print_3d_array(var, nx, ny, nz, str)
    !use mpi_mod
    use precision_mod
    implicit none
    integer(4), intent(in) :: nx, ny, nz
    real(wp), intent(in) :: var(nx, ny, nz)
    character(len=*),  intent(in) :: str

    integer(4) :: i, j, k

    write(*, *) str
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write(*, *) k, j, i, var(i, j, k)
        end do
      end do
    end do

    return
  end subroutine Print_3d_array

!===============================================================================
!===============================================================================
module code_performance_mod
  use precision_mod
  implicit none
  
  integer(4), parameter :: CPU_TIME_CODE_START = 1, &
                           CPU_TIME_ITER_START = 2, &
                           CPU_TIME_ITER_END   = 3, &
                           CPU_TIME_CODE_END   = 4

  real(wp), save :: t_code_start
  real(wp), save :: t_iter_start
  real(wp), save :: t_iter_end
  real(wp), save :: t_code_end

  private :: Convert_sec_to_hms
  public :: Call_cpu_time

  contains

  subroutine Convert_sec_to_hms (s, hrs, mins, secs)
    use precision_mod
    use parameters_constant_mod
    implicit none
    real(wp), intent(in) :: s
    integer, intent(out) :: hrs
    integer, intent(out) :: mins
    real(wp), intent(out) :: secs

    secs = s

    hrs = floor(secs / SIXTY / SIXTY)
    
    secs = secs - real(hrs, WP) * SIXTY * SIXTY
    mins = floor(secs / SIXTY)

    secs = secs - real(mins, WP) * SIXTY
    return
  end subroutine 

  subroutine Call_cpu_time(itype, nrsttckpt, niter, iter)
    use parameters_constant_mod
    use typeconvert_mod
    implicit none
    integer(4), intent(in) :: itype
    integer(4), intent(in) :: nrsttckpt, niter
    integer(4), intent(in), optional :: iter
    integer(4) :: hrs, mins
    real(wp) :: secs
    real(WP) :: t_total, t_elaspsed,t_remaining, t_aveiter, t_this_iter

    if(itype == CPU_TIME_CODE_START) then

      t_code_start = ZERO
      t_iter_start = ZERO
      t_iter_end   = ZERO
      t_code_end   = ZERO
      call cpu_time(t_code_start)

    else if (itype == CPU_TIME_ITER_START) then

      call cpu_time(t_iter_start)

      call Print_debug_start_msg ("Time Step = "//trim(int2str(iter))// &
          '/'//trim(int2str(niter-nrsttckpt)))

    else if (itype == CPU_TIME_ITER_END) then
      if(.not.present(iter)) call Print_error_msg("Error in calculating CPU Time.")
      call cpu_time(t_iter_end)

      t_this_iter = t_iter_end - t_iter_start
      call Print_debug_mid_msg ("  Code Performance Info :")
      call Print_debug_mid_msg ("    Time for this time step : " // &
          trim(real2str(t_this_iter))//' s')

      t_elaspsed  = t_iter_end - t_code_start
      call Convert_sec_to_hms (t_elaspsed, hrs, mins, secs)
      call Print_debug_mid_msg ("    Elaspsed Wallclock Time : "// &
           trim(int2str(hrs)) // ' h ' // &
           trim(int2str(mins)) // ' m ' // &
           trim(real2str(secs)) // ' s ')

      t_remaining= t_elaspsed / real(iter - nrsttckpt, wp) * real(niter - iter, wp)
      call Convert_sec_to_hms (t_remaining, hrs, mins, secs)
      call Print_debug_mid_msg ("    Remaning Wallclock Time : "// &
           trim(int2str(hrs)) // ' h ' // &
           trim(int2str(mins)) // ' m ' // &
           trim(real2str(secs)) // ' s ')

    else if (itype == CPU_TIME_CODE_END) then

      call cpu_time(t_code_end)
      t_total = t_code_end - t_code_start
      t_aveiter = t_total / real(niter - nrsttckpt, WP)
      call Print_debug_start_msg("CHAPSim Simulation is finished successfully.")
      call Print_debug_mid_msg ("Average wallclock time per step  : "// &
           trim(real2str(t_aveiter))//' s')
      
      call Convert_sec_to_hms (t_total, hrs, mins, secs)
      call Print_debug_mid_msg ("Total wallclock time of this run : "// &
           trim(int2str(hrs)) // ' h ' // &
           trim(int2str(mins)) // ' m ' // &
           trim(real2str(secs)) // ' s ')
      call Print_debug_start_msg(' ')
    else
    end if

    return
  end subroutine

end module



!===============================================================================
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

