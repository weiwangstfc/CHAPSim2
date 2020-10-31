module time_stepping_mod
  use precision_mod
  implicit none

  integer :: ntInner
  real(WP) :: tGamma(0:3)
  real(WP) :: tZeta (0:3)
  real(WP) :: tAlpha(0:3)

  public :: Set_timestepping_coefficients

contains

  subroutine Set_timestepping_coefficients()
    use parameters_constant_mod
    use input_mod, only : itimesteping

    if(itimesteping == ITIME_RK3) then
      
      ntInner = 3
      tGamma(0) = ONE
      tGamma(1) = EIGHT / FIFTEEN
      tGamma(2) = FIVE / TWELVE
      tGamma(3) = THREE / FOUR

      tZeta (0) = ZERO
      tZeta (1) = ZERO
      tZeta (2) = -SEVENTEEN / SIXTY
      tZeta (3) = -FIVE / TWELVE

    else if (itimesteping == ITIME_AB1) then

      ntInner = 1
      tGamma(0) = ONE
      tGamma(1) = ONEPFIVE
      tGamma(2) = ZERO
      tGamma(3) = ZERO

      tZeta (0) = ZERO
      tZeta (1) = -ZPFIVE
      tZeta (2) = ZERO
      tZeta (3) = ZERO

    else 

      ntInner = 0
      tGamma(:) = ZERO
      tZeta (:) = ZERO
 
    end if 
    
    tAlpha(:) = tGamma(:) + tZeta(:)

  end subroutine Set_timestepping_coefficients

end module time_stepping_mod