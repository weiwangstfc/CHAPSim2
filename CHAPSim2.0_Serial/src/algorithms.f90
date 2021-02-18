module tridiagonal_matrix_algorithm
  
  implicit none
contains

  subroutine Prepare_TDMA_coeffs(a, b, c, d, n)
    use math_mod
    use parameters_constant_mod, only: ONE
    use precision_mod
    implicit none
    integer(4), intent(in) :: n
    real(WP), intent(in)    :: a(n), b(n)
    real(WP), intent(inout) :: c(n)
    real(WP), intent(out)   :: d(n)
    integer(4) :: i

    ! prepare coefficients
    c(1) = c(1) / b(1)
    
    do i = 2, n
      d(i) = ONE / ( b(i) - a(i) * c(i - 1) )
      if (i < n) c(i) = c(i) * d(i)
    end do

    return
  end subroutine Prepare_TDMA_coeffs



  subroutine Solve_TDMA_basic(x, a, b, c, d, n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solution of a tridiagnal system of n equations of the form
! 
!  a(i) * x(i-1) + b(i) * x(i) + c(i) * x(i+1) = R(i), i = 1, ..., n
!  a(1) and c(n) are not used. 
!  The solution x(i) is restored in R(i).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use precision_mod
    implicit none

    integer(4), intent(in) :: n
    real(WP), intent(inout) :: x(n) ! R in, X out.
    real(WP), intent(in) :: a(n), b(n)
    real(WP), intent(in) :: c(n), d(n)
    integer(4) :: i
    
    x(1) = x(1) / b(1)
    
    ! forward substitution
    do i = 2, n
      x(i) = ( x(i) - a(i) * x(i-1) ) * d(i)
    end do

    ! backward substitution
    do i = n-1, 1, -1
      x(i) = x(i) - c(i) * x(i+1)  
    end do

    return
  end subroutine Solve_TDMA_basic



  subroutine Solve_TDMA_cyclic(x, a, b, c, d, n)

    use precision_mod
    implicit none
    
    integer(4), intent(in) :: n
    real(WP), intent(inout) :: x(n) ! R in, X out.
    real(WP), intent(in) :: a(n), b(n)
    real(WP), intent(in) :: c(n), d(n)
    real(WP) :: x1(n)

    call Solve_TDMA_basic(x(1:n-1), a(1:n-1), b(1:n-1), c(1:n-1), d(1:n-1), n-1)
    
    x1(:) = 0.0
    x1(1) = - a(1)
    x1(n-1) =  - c(n-1)
    call Solve_TDMA_basic(x1(1:n-1), a(1:n-1), b(1:n-1), c(1:n-1), d(1:n-1), n-1)

    x(n) = (x(n) - c(n) * x(1) - a(n) * x(n-1)) / &
           (b(n) + c(n) * x1(1) + a(n) * x1(n-1))
    x(1:n-1) = x(1:n-1) + x1(1:n-1) * x(n)
    return
  end subroutine Solve_TDMA_cyclic



  subroutine test_TDMA_noncyclic
    use precision_mod
    implicit none
    integer(4), parameter :: n = 10
    real(WP) :: a(n), b(n), c(n), d(n), r(n)
    real(WP) :: ref(n)
    integer(4) :: i

    a(1: n) = [3.0, 1.0, 1.0, 7.0, 6.0, 3.0, 8.0, 6.0, 5.0, 4.0]
    b(1: n) = [2.0, 3.0, 3.0, 2.0, 2.0, 4.0, 1.0, 2.0, 4.0, 5.0]
    c(1: n) = [1.0, 2.0, 1.0, 6.0, 1.0, 3.0, 5.0, 7.0, 3.0, 5.0]
    r(1: n) = [1.0, 2.0, 6.0, 34.0, 10.0, 1.0, 4.0, 22.0, 25.0, 3.0]

    d(:) = 0.0

    call Prepare_TDMA_coeffs(a(:), b(:), c(:), d(:), n)

    call Solve_TDMA_basic(r(:), a(:), b(:), c(:), d(:), n)

    
    ! data output
    ref=[1.0, -1.0, 2.0, 1.0, 3.0, -2.0, 0.0, 4.0, 2.0, -1.0]
    write(*, '(A)') "TDMA Basic: I, X, Xref, A, B, C, D" 
    do i = 1, n
      write(*, '(I3, 6F8.4)') i, r(i), ref(i), a(i), b(i), c(i), d(i)
    end do
    
    return
  end subroutine test_TDMA_noncyclic


  subroutine test_TDMA_cyclic
    use precision_mod
    implicit none
    integer(4), parameter :: n = 10
    real(WP) :: a(n), b(n), c(n), d(n), r(n)
    real(WP) :: ref(n)
    integer(4) :: i

    a(1: n) = [3.0, 1.0, 1.0, 7.0, 6.0, 3.0, 8.0, 6.0, 5.0, 4.0]
    b(1: n) = [2.0, 3.0, 3.0, 2.0, 2.0, 4.0, 1.0, 2.0, 4.0, 5.0]
    c(1: n) = [1.0, 2.0, 1.0, 6.0, 1.0, 3.0, 5.0, 7.0, 3.0, 5.0]
    r(1: n) = [1.0, 2.0, 6.0, 34.0, 10.0, 1.0, 4.0, 22.0, 25.0, 3.0]

    d(:) = 0.0

    call Prepare_TDMA_coeffs(a(1:n-1), b(1:n-1), c(1:n-1), d(1:n-1), n-1)
    call Solve_TDMA_cyclic(r(:), a(:), b(:), c(:), d(:), n)

    ! data output
    ref=[518663./174746., -299297./174746., 182180./87373., &
         5419./3718., 480243./174746., -370592./87373., 566251./174746., &
         1212441./174746., -76./47., -187761./174746.]
    write(*, '(A)') "TDMA-cyclic:  I, X, Xref, A, B, C, D" 
    do i = 1, n
      write(*, '(I3, 6F8.4)') i, r(i), ref(i), a(i), b(i), c(i), d(i)
    end do

    return
  end subroutine test_TDMA_cyclic


end module tridiagonal_matrix_algorithm


subroutine Test_algorithms
  use tridiagonal_matrix_algorithm

  call test_TDMA_noncyclic
  call test_TDMA_cyclic

  return
end subroutine