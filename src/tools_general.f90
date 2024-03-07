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

    write (*, *) "    "//msg
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
module decomp_operation_mod
  implicit none 
contains
    
  function is_decomp_same ( a, b ) result(f)
    use decomp_2d
    type(DECOMP_INFO), intent(in) :: a, b
    logical :: f
    integer :: i

    f = .true.
    do i = 1, 3
      if(a%xst(i) /= b%xst(i)) f = .false.
      if(a%xen(i) /= b%xen(i)) f = .false.
      if(a%yst(i) /= b%yst(i)) f = .false.
      if(a%yen(i) /= b%yen(i)) f = .false.
      if(a%zst(i) /= b%zst(i)) f = .false.
      if(a%zen(i) /= b%zen(i)) f = .false.
    end do
  end function
end module

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
  integer :: cpu_nfre 

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


module cubic_spline_interpolation


  public :: cubic_spline
  public :: spline_interpolation


contains

  !**********************************************************************************************************************************
  subroutine cubic_spline (n, x, y, b, c, d)
  !---------------------------------------------------------------------
  !     this subroutine calculates the coefficients b, c, d of a cubic
  !     spline to best approximate a discreet fonction given by n points
  !
  !     inputs:
  !     n       number of given points
  !     x, y    vectors of dimension n, storing the coordinates
  !             of function f(x)
  !
  !     outputs:
  !     b,c, d   vectors of dimension n, storing the coefficients
  !             of the cubic spline
  !     function:
  !     y =  x
  !     reference:
  !     forsythe, g.e. (1977) computer methods for mathematical
  !     computations. prentice - hall, inc.
  !---------------------------------------------------------------------
      use precision_mod
      implicit none
      integer(4), intent(in) :: n
      real(wp), intent(in) :: x(n), y(n)
      real(wp), intent(out) :: b(n), c(n), d(n)

      integer(4) :: nm1, i, l
      real(wp) :: t


      if (n < 2) return
      if (n < 3) then
          b(1) = (y(2) - y(1)) / (x(2) - x(1))
          c(1) = 0.0_wp
          d(1) = 0.0_wp
          b(2) = b(1)
          c(2) = 0.0_wp
          d(2) = 0.0_wp
          return
      end if

      ! step 1: preparation
      !        build the tridiagonal system
      !        b (diagonal), d (upperdiagonal), c (second member)
      nm1 = n - 1

      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1)) / d(1)
      do i = 2, nm1
          d(i) = x(i + 1) - x(i)
          b(i) = 2.0_wp * (d(i - 1) + d(i))
          c(i + 1) = (y(i + 1) - y(i)) / d(i)
          c(i) = c(i + 1) - c(i)
      end do

      ! step 2: end conditions
      !     conditions at limits
      !     third derivatives obtained by divided differences
      b(1) = - d(1)
      b(n) = - d(n - 1)
      c(1) = 0.0_wp
      c(n) = 0.0_wp

      if(n /= 3) then
          c(1) = c(3) / (x(4) - x(2)) - c(2) / (x(3) - x(1))
          c(n) = c(n - 1) / (x(n) - x(n - 2)) - c(n - 2) / (x(n - 1) - x(n - 3))
          c(1) = c(1) * d(1) * d(1) / (x(4) - x(1))
          c(n) = - c(n) * d(n - 1)**2 / (x(n) - x(n - 3))
      end if

      ! step 3:     forward elimination
      do i = 2, n
          t = d(i - 1) / b(i - 1)
          b(i) = b(i) - t * d(i - 1)
          c(i) = c(i) - t * c(i - 1)
      end do

      !step 4:     back substitution
      c(n) = c(n) / b(n)
      do  l = 1, nm1
          i = n - l
          c(i) = (c(i) - d(i) * c(i + 1)) / b(i)
      end do

      !step 5: coefficients of 3rd degree polynomial
      b(n) = (y(n) - y(nm1)) / d(nm1) + d(nm1) * (c(nm1) + 2.0_wp * c(n))
      do  i = 1, nm1
          b(i) = (y(i + 1) - y(i)) / d(i) - d(i) * (c(i + 1) + 2.0_wp * c(i))
          d(i) = (c(i + 1) -c(i)) / d(i)
          c(i) = 3.0_wp * c(i)
      end do
      c(n) = 3.0_wp * c(n)
      d(n) = d(nm1)

      return
  end subroutine

  !**********************************************************************************************************************************
  function spline_interpolation(n, yprofile, b, c, d, y) result(eval)
    use precision_mod
    implicit none
    integer, intent(in) :: n
    real(WP), intent(in) :: yprofile(n)
    real(wp), intent(in) :: b(n), c(n), d(n)
    real(WP), intent(in) :: y
    real(WP) :: eval
    
    integer :: i, j, k
    real(WP) :: dy


    !*
    !  binary search for for i, such that x(i) <= u <= x(i + 1)
    !*
    i = 1
    j = n + 1
    do while (j > i + 1)
        k = (i + j) / 2
        if(y < yprofile(k)) then
            j = k
        else
            i = k
        end if
    end do
    !*
    !  evaluate spline interpolation
    !*
    dy = y - yprofile(i)
    eval = yprofile(i) + dy * (b(i) + dy * (c(i) + dy * d(i)))

  end function

end module


subroutine map_1d_profile_to_case(nin, yin, uin, nout, ycase, ucase)
  use cubic_spline_interpolation
  use precision_mod
  implicit none
  integer, intent(in) :: nin
  real(WP), dimension(nin), intent(in) :: yin
  real(WP), dimension(nin), intent(in) :: uin
  integer, intent(in) :: nout
  real(WP), dimension(nout), intent(in)  :: ycase
  real(WP), dimension(nout), intent(out) :: ucase

  integer :: i
  real(WP), allocatable :: cs_b(:), cs_c(:), cs_d(:)

  allocate(cs_b(nin))
  allocate(cs_c(nin))
  allocate(cs_d(nin))

  call cubic_spline (nin, yin, uin, cs_b, cs_c, cs_d)


  do i = 1, nout
    ucase(i) = spline_interpolation(nin, uin, cs_b, cs_c, cs_d, ycase(i))
  end do 

  deallocate(cs_b)
  deallocate(cs_c)
  deallocate(cs_d)

  return
end subroutine


!==========================================================================================================

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



subroutine wrt_3d_pt_debug(var, dtmp, iter, irk, str, loc)
  use precision_mod
  use udf_type_mod
  implicit none 
  type(DECOMP_INFO), intent(in) :: dtmp
  real(wp), intent(in)     :: var(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3))
  character(*), intent(in) :: str
  character(*), intent(in) :: loc
  
  integer, intent(in) :: iter, irk

  integer, parameter :: npt = 4
  integer, parameter :: nfil = 20
  integer  :: nid(4, 3), a(12)

  character(1) :: pntim
  character(128) :: flnm
  logical :: file_exists
  integer :: n, i, j, k, jj, kk

 ! based on x pencil

  a = (/8, 16, 32, 40, 8, 16, 32, 40, 8, 16, 32, 40/)
  nid = reshape(a, (/4, 3/))
  do n = 1, npt
      write(pntim,'(i1.1)') n
      flnm = 'chapsim2_p'//pntim//'_'//trim(str)//'.dat'   
      do k =1, dtmp%xsz(3)
          kk = dtmp%xst(3) + k - 1
          if(kk == nid(n, 3)) then
              do j = 1, dtmp%xsz(2)
                  jj = dtmp%xst(2) + j - 1
                  if(jj == nid(n, 2)) then
                    do i = 1, dtmp%xsz(1)
                        if(i == nid(n, 1)) then
                          inquire(file=trim(adjustl(flnm)), exist=file_exists) 
                          if(file_exists) then
                            open(nfil+n,file=trim(adjustl(flnm)), position='append')
                            !write(nfil+n,*) '# iter = ', iter
                          else
                            open(nfil+n,file=trim(adjustl(flnm)) )
                            !write(nfil+n,*) '# iter = ', iter
                          end if
                          write(nfil+n, *) trim(str), trim(loc), iter, irk, i, jj, kk, var(i, j, k)
                          close(nfil+n)
                        end if
                      end do
                  end if
              end do
          end if
      end do
    end do

return

end subroutine


subroutine wrt_3d_all_debug(var, dtmp, iter, irk, str, loc)
  use precision_mod
  use udf_type_mod
  implicit none 
  type(DECOMP_INFO), intent(in) :: dtmp
  real(wp), intent(in)     :: var(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3))
  character(*), intent(in) :: str
  character(*), intent(in) :: loc
  
  integer, intent(in) :: iter, irk
  integer, parameter :: nfil = 20

  character(128) :: flnm
  logical :: file_exists
  integer :: n, i, j, k, jj, kk


  flnm = 'chapsim2_'//trim(str)//'.dat'   
  inquire(file=trim(adjustl(flnm)), exist=file_exists) 
  if(file_exists) then
    open(nfil,file=trim(adjustl(flnm)), position='append')
    !write(nfil+n,*) '# iter = ', iter
  else
    open(nfil,file=trim(adjustl(flnm)) )
    !write(nfil+n,*) '# iter = ', iter
  end if

  do j = 1, dtmp%xsz(2)
    jj = dtmp%xst(2) + j - 1
    do k =1, dtmp%xsz(3)
      kk = dtmp%xst(3) + k - 1
      do i = 1, dtmp%xsz(1)
        write(nfil, *) jj, kk, i, var(i, j, k)  
      end do
    end do
  end do
  close(nfil)

  return
end subroutine

