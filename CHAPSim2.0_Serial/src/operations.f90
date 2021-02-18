module operations
  use geometry_mod
  implicit none

  integer(4), parameter :: I1Deri2CD = 1 ! 2nd order central difference
  integer(4), parameter :: I1Deri4CD = 2 ! 4th order central difference
  integer(4), parameter :: I1Deri4CP = 3 ! 4th order compact scheme, Pade 
  integer(4), parameter :: I1Deri6CP = 4 ! 6th order compact scheme

  integer(4), parameter :: I2Deri2CD = 1
  integer(4), parameter :: I2Deri4CD = 2
  integer(4), parameter :: I2Deri4CP = 3
  integer(4), parameter :: I2Deri6CP = 4

  integer(4) :: compact1FD(3, 4) ! collocated, 1st derivative
  integer(4) :: alpha1FDs, a1FDs, b1FDs ! staggered, 1st derivative
  integer(4) :: alphaIP,   aIP,   bIP   ! mid-point interpolation


  private
  public :: Interpolate_stencil

  abstract interface
    subroutine Get_derivative_x()
    end subroutine Get_derivative_x

    subroutine Get_derivative_y()
    end subroutine Get_derivative_y

    subroutine Get_derivative_z()
    end subroutine Get_derivative_z

    subroutine Get_derivative_yy()
    end subroutine Get_derivative_yy
  end interface

  procedure (Get_derivative_x), pointer :: derx, derxx
  procedure (Get_derivative_x) derx_00, derxx_00

  procedure (Get_derivative_y), pointer :: dery
  procedure (Get_derivative_y) dery_00, dery_11

  procedure (Get_derivative_z), pointer :: derz
  procedure (Get_derivative_z) derz_00, derzz_00

  procedure (Get_derivative_yy), pointer :: deryy
  procedure (Get_derivative_yy) deryy_00, deryy_11

contains

  subroutine set_scheme_coefficient(fip, f1, f2)

    integer(4), intent(in) :: fip, f1, f2

    ! 
    ! Finite Difference collocated:
    ! alpha * f'(i-1) + f'(i) + f'(i+1) = [B/4h] * ( f(i+2) - f(i-2) ) + [A/2h] * ( f(i+1) - f(i-1) )
    
    ! Finite Difference staggered:
    ! alpha * f'(i-1) + f'(i) + f'(i+1) = [B/3h] * ( f(i'+2) - f(i'-1) ) + [A/h] * ( f(i'+1) - f(i') )
    ! alpha * f'(i'-1) + f'(i') + f'(i'+1) = [B/3h] * ( f(i+2) - f(i-1) ) + [A/h] * ( f(i+1) - f(i) )
    
    ! Mid-point interpolation
    ! alpha * f(i'-1) + f(i') + f(i'+1) = [B/2] * ( f(i+1) + f(i-2) ) + [A/2] * ( f(i) + f(i-1) )
    ! alpha * f(i-1) + f(i) + f(i+1) = [B/2] * ( f(i'+2) + f(i'-1) ) + [A/2] * ( f(i') + f(i+1) )
    if(i1DerScheme == I1Deri2CD) then ! 2nd order central difference
      alphaiFDc = ZERO
      aiFDc = HALF
      biFDc = ZERO

      alphaiFDs = ZERO
      aiFDs = HALF
      biFDs = ZERO

      alphaiIP = ZERO
      aiIP = HALF
      biIP = ZERO

    else if (i1DerScheme == I1Deri4CD) then
      alphaiFDc = ZERO
      aiFDc = TWO / THREE
      biFDc = - ONE / TWELVE

      alphaiFDs = ZERO
      aiFDs = NINE / EIGHT
      biFDs = - ONE / TWENTYFOUR

      alphaiIP = ZEROS
      aiIP = NINE / SIXTEEN
      biIP = - ONE / SIXTEEN

    else if (i1DerScheme == I1Deri4CP) then
      alpha1FDc = ONE / FOUR 
      a1FDc = THREE / FOUR
      b1FDc = ZERO

      alpha1FDs = ONE / TWENTYTWO
      a1FDs = TWELVE / ELEVEN
      b1FDs = ZERO

      alphaIP = ONE / SIX
      aIP = TWO / THREE
      bIP = ZERO

    else if (i1DerScheme == I1Deri6CP) then
      alpha1FDc = ONE / THREE
      a1FDc = SEVEN / NINE
      b1FDc = ONE / THIRTYSIX

      alpha1FDs = NINE / SIXTYTWO
      a1FDs = SIXTYTHREE / SIXTYTWO
      b1FDs = SEVENTEEN / (SIXTYTWO * THREE)

      alphaIP = THREE / TEN
      aIP = THREE / FOUR
      bIP = ONE / TWENTY
      
    else 
      call Print_error_msg("No such scheme defined.")
    end if

    ! add boundary

  end subroutine set_compact_coefficients

  subroutine derx_00




  end subroutine derx_00



end module