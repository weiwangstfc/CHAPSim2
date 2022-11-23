!==========================================================================================================
  subroutine Print_error_msg(msg)
    !use iso_fortran_env
    implicit none
    character(len=*), intent(IN) :: msg
    
    write (*, *) 'ERROR: ' // msg

    write (*, *) 'Code is terminated in processor = '
    STOP

    return
  end subroutine Print_error_msg
!==========================================================================================================
  subroutine Print_warning_msg(msg)
    !use iso_fortran_env
    implicit none
    character(len=*), intent(IN) :: msg
    
    write (*, *) 'WARNNING: ' // msg

    return
  end subroutine Print_warning_msg
  !==========================================================================================================
  subroutine Print_debug_start_msg(msg)
    !use iso_fortran_env
    implicit none
    character(len=*), intent(IN) :: msg

    write (*, *) "=========================================================================================================="
    write (*, *) msg

    return
  end subroutine Print_debug_start_msg
!==========================================================================================================
  subroutine Print_debug_mid_msg(msg)
    !use iso_fortran_env
    implicit none
    character(len=*), intent(IN) :: msg

    write (*, *) msg
    return
  end subroutine Print_debug_mid_msg
!==========================================================================================================
  subroutine Print_debug_end_msg
    !use iso_fortran_env
    implicit none

    write (*, *) "... done."
    return
  end subroutine Print_debug_end_msg

!==========================================================================================================
  subroutine Print_3d_array(var, nx, ny, nz, str)
    use precision_mod
    !use iso_fortran_env
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(wp), intent(in) :: var(nx, ny, nz)
    character(len=*),  intent(in) :: str

    integer :: i, j, k

    write (*, *) str
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write (*, *) k, j, i, var(i, j, k)
        end do
      end do
    end do

    return
  end subroutine Print_3d_array

!==========================================================================================================
!==========================================================================================================
module code_performance_mod
  use precision_mod
  implicit none
  
  integer, parameter :: CPU_TIME_CODE_START = 1, &
                        CPU_TIME_STEP_START = 2, &
                        CPU_TIME_ITER_START = 3, &
                        CPU_TIME_ITER_END   = 4, &
                        CPU_TIME_STEP_END   = 5, &
                        CPU_TIME_CODE_END   = 6

  real(wp), save :: t_code_start
  real(wp), save :: t_step_start
  real(wp), save :: t_iter_start
  real(wp), save :: t_iter_end
  real(wp), save :: t_step_end
  real(wp), save :: t_code_end

  private :: Convert_sec_to_hms
  public :: call_cpu_time

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

  subroutine call_cpu_time(itype, iterfrom, niter, iter)
    use parameters_constant_mod
    use typeconvert_mod
    use mpi_mod
    use decomp_2d
    implicit none
    integer, intent(in) :: itype
    integer, intent(in) :: iterfrom, niter
    integer, intent(in), optional :: iter
    integer :: hrs, mins
    real(wp) :: secs
    real(WP) :: t_total, t_elaspsed, t_remaining, t_aveiter, t_this_iter, t_preparation, t_postprocessing
    real(WP) :: t_total0, t_elaspsed0,t_remaining0, t_aveiter0, t_this_iter0, t_preparation0, t_postprocessing0
!----------------------------------------------------------------------------------------------------------
    if(itype == CPU_TIME_CODE_START) then
      t_code_start = ZERO
      call cpu_time(t_code_start)
!----------------------------------------------------------------------------------------------------------
    else if (itype == CPU_TIME_STEP_START) then
      t_step_start = ZERO
      call cpu_time(t_step_start)
      t_preparation = t_step_start - t_code_start
      call mpi_barrier(MPI_COMM_WORLD, ierror)
      call mpi_allreduce(t_preparation, t_preparation0, 1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
      if(nrank == 0) call Print_debug_start_msg ("  Code Performance Info :")
      if(nrank == 0) call Print_debug_mid_msg ("    Time for code preparation : " // &
          trim(real2str(t_preparation0))//' s')
!----------------------------------------------------------------------------------------------------------
    else if (itype == CPU_TIME_ITER_START) then
      t_iter_start = ZERO
      call cpu_time(t_iter_start)
      if(nrank == 0) call Print_debug_start_msg ("Time Step = "//trim(int2str(iter))// &
          '/'//trim(int2str(niter-iterfrom)))
!----------------------------------------------------------------------------------------------------------
    else if (itype == CPU_TIME_ITER_END) then
      if(.not.present(iter)) call Print_error_msg("Error in calculating CPU Time.")
      call cpu_time(t_iter_end)

      t_this_iter = t_iter_end - t_iter_start
      call mpi_barrier(MPI_COMM_WORLD, ierror)
      call mpi_allreduce(t_this_iter, t_this_iter0, 1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
      if(nrank == 0) call Print_debug_mid_msg ("  Code Performance Info :")
      if(nrank == 0) call Print_debug_mid_msg ("    Time for this time step : " // &
          trim(real2str(t_this_iter0))//' s')

      t_elaspsed  = t_iter_end - t_step_start
      call mpi_barrier(MPI_COMM_WORLD, ierror)
      call mpi_allreduce(t_elaspsed, t_elaspsed0, 1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
      call Convert_sec_to_hms (t_elaspsed0, hrs, mins, secs)
      if(nrank == 0) call Print_debug_mid_msg ("    Elaspsed Wallclock Time : "// &
           trim(int2str(hrs)) // ' h ' // &
           trim(int2str(mins)) // ' m ' // &
           trim(real2str(secs)) // ' s ')

      t_aveiter   = t_elaspsed / real(iter - iterfrom, WP)
      t_remaining = t_aveiter * real(niter - iter, wp)
      call mpi_barrier(MPI_COMM_WORLD, ierror)
      call mpi_allreduce(t_aveiter,   t_aveiter0,   1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
      call mpi_allreduce(t_remaining, t_remaining0, 1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
      
      call Convert_sec_to_hms (t_remaining0, hrs, mins, secs)

      if(nrank == 0) then
        call Print_debug_mid_msg ("    Remaning Wallclock Time : "// &
           trim(int2str(hrs)) // ' h ' // &
           trim(int2str(mins)) // ' m ' // &
           trim(real2str(secs)) // ' s ')
        call Print_debug_mid_msg ("    Moving averaged time per iteration  : "// &
           trim(real2str(t_aveiter0))//' s')
        
      end if
!----------------------------------------------------------------------------------------------------------
    else if (itype == CPU_TIME_STEP_END) then

      call cpu_time(t_step_end)
      t_total = t_step_end - t_step_start
      t_aveiter= t_total / real(niter - iterfrom, WP)
      call mpi_barrier(MPI_COMM_WORLD, ierror)
      call mpi_allreduce(t_total,   t_total0,   1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
      call mpi_allreduce(t_aveiter, t_aveiter0, 1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
      
      call Convert_sec_to_hms (t_total0, hrs, mins, secs)
      if(nrank == 0) then
        call Print_debug_start_msg ("  Code Performance Info :")
        call Print_debug_mid_msg   ("    Averaged time per iteration  : "// &
           trim(real2str(t_aveiter0))//' s')
        call Print_debug_mid_msg ("    Wallclock time of all iterations : "// &
           trim(int2str(hrs)) // ' h ' // &
           trim(int2str(mins)) // ' m ' // &
           trim(real2str(secs)) // ' s ')
      end if
!----------------------------------------------------------------------------------------------------------
    else if (itype == CPU_TIME_CODE_END) then

      call cpu_time(t_code_end)
      t_total  = t_code_end - t_code_start
      t_postprocessing = t_code_end - t_step_end
      call mpi_barrier(MPI_COMM_WORLD, ierror)
      call mpi_allreduce(t_total,          t_total0,          1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
      call mpi_allreduce(t_postprocessing, t_postprocessing0, 1, MPI_REAL_WP, MPI_MAX, MPI_COMM_WORLD, ierror)
      
      call Convert_sec_to_hms (t_total0, hrs, mins, secs)
      if(nrank == 0) then
        call Print_debug_start_msg ("  Code Performance Info :")
        call Print_debug_mid_msg    ("    Wallclock time for postprocessing : "// &
           trim(real2str(t_postprocessing0))//' s')
        call Print_debug_mid_msg ("    Total wallclock time of this run : "// &
           trim(int2str(hrs)) // ' h ' // &
           trim(int2str(mins)) // ' m ' // &
           trim(real2str(secs)) // ' s ')
        call Print_debug_mid_msg("CHAPSim Simulation is finished successfully.")
      end if
    else
    end if

    return
  end subroutine

end module



!==========================================================================================================
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
    !    Input/output, integer seed.
    !    IF seed is zero on input, THEN you're asking this routine to come up
    !    with a seed value, whICh is RETURNed as output.
    !    IF seed is nonzero on input, THEN you're asking this routine to
    !    USE the input value of seed to initialize the random number generator,
    !    and seed is not changed on output.
    !
    implicit none
    !
    integer :: count
    integer :: count_max
    integer :: count_rate
    logical, parameter :: debug = .false.
    integer :: i
    integer :: seed
    integer, allocatable :: seed_vector(:)
    integer :: seed_size
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
    !    Input, integer N, the number of entries in the vector.
    !
    !    Output, REAL(WP) A(N), the vector of randomly chosen values.
    !
    implicit none
    !
    integer n
    !
    real(wp) a(n)
    real(wp) ahi
    real(wp) alo
    integer i
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

