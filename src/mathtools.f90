module math_mod
  use precision_mod

  interface sqrt_wp
    module procedure sqrt_sp
    module procedure sqrt_dp
  end interface sqrt_wp

  interface tanh_wp
    module procedure tanh_sp
    module procedure tanh_dp
  end interface tanh_wp

contains

  function sqrt_sp ( r ) result(d)
    real(kind = SP), intent(in) :: r
    real(kind = SP) :: d
    d = sqrt ( r )
  end function

  function sqrt_dp ( r ) result (d)
    real(kind = DP), intent(in) :: r
    real(kind = DP) :: d
    d = dsqrt ( r ) 
  end function


  function tanh_sp ( r ) result(d)
    real(kind = SP), intent(in) :: r
    real(kind = SP) :: d
    d = tanh ( r )
  end function

  function tanh_dp ( r ) result (d)
    real(kind = DP), intent(in) :: r
    real(kind = DP) :: d
    d = dtanh ( r ) 
  end function

end module math_mod


