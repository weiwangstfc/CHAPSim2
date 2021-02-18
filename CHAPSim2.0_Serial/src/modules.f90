!##############################################################################
module precision_mod
  
  integer, parameter :: I4 = selected_int_kind( 4 )
  integer, parameter :: I8 = selected_int_kind( 8 )
  integer, parameter :: I15 = selected_int_kind( 15 )
  integer, parameter :: SP = selected_real_kind( p = 6, r = 37 )
  integer, parameter :: DP = selected_real_kind( p = 15, r = 307 )
  integer, parameter :: QP = selected_real_kind( p = 33, r = 4931 )

  integer, parameter :: WP = DP

end module precision_mod

!##############################################################################
module parameters_constant_mod
  use precision_mod

  real(WP), parameter :: ZPONE     = 0.1_WP
  real(WP), parameter :: ZPTWO     = 0.2_WP
  real(WP), parameter :: ZPTHREE   = 0.3_WP
  real(WP), parameter :: ZPFOUR    = 0.4_WP
  real(WP), parameter :: ZPFIVE    = 0.5_WP
  real(WP), parameter :: ZPSIX     = 0.6_WP
  real(WP), parameter :: ZPSEVEN   = 0.7_WP
  real(WP), parameter :: ZPEIGHT   = 0.8_WP
  real(WP), parameter :: ZPNINE    = 0.9_WP

  real(WP), parameter :: HALF      = 0.5_WP
  real(WP), parameter :: ZERO      = 0.0_WP
  real(WP), parameter :: ONE       = 1.0_WP
  real(WP), parameter :: ONEPFIVE  = 1.5_WP
  real(WP), parameter :: TWO       = 2.0_WP
  real(WP), parameter :: THREE     = 3.0_WP
  real(WP), parameter :: FOUR      = 4.0_WP
  real(WP), parameter :: FIVE      = 5.0_WP
  real(WP), parameter :: SIX       = 6.0_WP
  real(WP), parameter :: SEVEN     = 7.0_WP
  real(WP), parameter :: EIGHT     = 8.0_WP
  real(WP), parameter :: NINE      = 9.0_WP

  real(WP), parameter :: TEN       = 10.0_WP
  real(WP), parameter :: ELEVEN    = 11.0_WP
  real(WP), parameter :: TWELVE    = 12.0_WP 
  real(WP), parameter :: SIXTEEN   = 16.0_WP
  real(WP), parameter :: FIFTEEN   = 15.0_WP
  real(WP), parameter :: SEVENTEEN = 17.0_WP

  real(WP), parameter :: TWENTYTWO = 22.0_WP
  real(WP), parameter :: TWENTYFOUR= 24.0_WP

  real(WP), parameter :: THIRTYSIX = 36.0_WP

  real(WP), parameter :: SIXTY     = 60.0_WP
  real(WP), parameter :: SIXTYTWO  = 62.0_WP
  real(WP), parameter :: SIXTYTHREE= 63.0_WP

  real(WP), parameter :: MINP      = 1.0E-20_WP
  real(WP), parameter :: MAXP      = 1.0E20_WP

  real(WP), parameter :: MINN      = -1.0E20_WP
  real(WP), parameter :: MAXN      = -1.0E-20_WP


  real(WP), parameter :: TRUNCERR = 1.0E-15_WP

  real(WP), parameter :: PI = dacos( -ONE )
  real(WP), parameter :: TWOPI = TWO * dacos( -ONE )

  real(WP), parameter, dimension(3, 3) :: KRONECKER_DELTA = &
                                            reshape( (/ &
                                            ONE, ZERO, ZERO, &
                                            ZERO, ONE, ZERO, &
                                            ZERO, ZERO, ONE  /), &
                                            (/3, 3/) )
end module parameters_constant_mod


module udf_type_mod
  use precision_mod
  type domain_t
    integer :: case
    integer :: bcx(2)
    integer :: bcy(2)
    integer :: bcz(2)
    logical :: is_periodic(3)
    integer :: np(3) ! geometric points
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

end module

!##############################################################################
module math_mod
  use precision_mod
  use parameters_constant_mod, only: ONE, ZERO, MINP

  interface sqrt_wp
    module procedure sqrt_sp
    module procedure sqrt_dp
  end interface sqrt_wp

  interface tanh_wp
    module procedure tanh_sp
    module procedure tanh_dp
  end interface tanh_wp

  interface abs_wp
    module procedure abs_sp
    module procedure abs_dp
  end interface abs_wp

  interface sin_wp
    module procedure sin_sp
    module procedure sin_dp
  end interface sin_wp

  interface cos_wp
    module procedure cos_sp
    module procedure cos_dp
  end interface cos_wp

  interface tan_wp
    module procedure tan_sp
    module procedure tan_dp
  end interface tan_wp

  interface atan_wp
    module procedure atan_sp
    module procedure atan_dp
  end interface atan_wp
  
contains

  ! abs
  pure function abs_sp ( r ) result(d)
  real(kind = SP), intent(in) :: r
  real(kind = SP) :: d
    d = abs ( r )
  end function

  pure function abs_dp ( r ) result (d)
  real(kind = DP), intent(in) :: r
  real(kind = DP) :: d
    d = dabs ( r ) 
  end function

  ! sqrt
  pure function sqrt_sp ( r ) result(d)
    real(kind = SP), intent(in) :: r
    real(kind = SP) :: d
    d = sqrt ( r )
  end function

  pure function sqrt_dp ( r ) result (d)
    real(kind = DP), intent(in) :: r
    real(kind = DP) :: d
    d = dsqrt ( r ) 
  end function

  ! sin
  pure function sin_sp ( r ) result(d)
    real(kind = SP), intent(in) :: r
    real(kind = SP) :: d
    d = sin ( r )
  end function

  pure function sin_dp ( r ) result (d)
    real(kind = DP), intent(in) :: r
    real(kind = DP) :: d
    d = dsin ( r ) 
  end function

  ! cos
  pure function cos_sp ( r ) result(d)
    real(kind = SP), intent(in) :: r
    real(kind = SP) :: d
    d = cos ( r )
  end function

  pure function cos_dp ( r ) result (d)
    real(kind = DP), intent(in) :: r
    real(kind = DP) :: d
    d = dcos ( r ) 
  end function

  ! tanh
  pure function tanh_sp ( r ) result(d)
    real(kind = SP), intent(in) :: r
    real(kind = SP) :: d
    d = tanh ( r )
  end function

  pure function tanh_dp ( r ) result (d)
    real(kind = DP), intent(in) :: r
    real(kind = DP) :: d
    d = dtanh ( r ) 
  end function

  ! tan
  pure function tan_sp ( r ) result(d)
    real(kind = SP), intent(in) :: r
    real(kind = SP) :: d
    d = tan ( r )
  end function

  pure function tan_dp ( r ) result (d)
    real(kind = DP), intent(in) :: r
    real(kind = DP) :: d
    d = tan ( r ) 
  end function

  ! atan
  pure function atan_sp ( r ) result(d)
    real(kind = SP), intent(in) :: r
    real(kind = SP) :: d
    d = atan ( r )
  end function

  pure function atan_dp ( r ) result (d)
    real(kind = DP), intent(in) :: r
    real(kind = DP) :: d
    d = atan ( r ) 
  end function

  ! heaviside_step
  pure function heaviside_step ( r ) result (d)
    real(kind = WP), intent(in) :: r
    real(kind = WP) :: d
    d = ZERO
    if (r > MINP) d = ONE
  end function

end module math_mod


