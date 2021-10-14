!##############################################################################
module precision_mod
  use decomp_2d, only : mytype
  implicit none

  integer, parameter :: I4 = selected_int_kind( 4 )
  integer, parameter :: I8 = selected_int_kind( 8 )
  integer, parameter :: I15 = selected_int_kind( 15 )
  integer, parameter :: S6P = selected_real_kind( p = 6, r = 37 )
  integer, parameter :: D15P = selected_real_kind( p = 15, r = 307 )
  integer, parameter :: Q33P = selected_real_kind( p = 33, r = 4931 )

  integer, parameter :: WP = D15P
  !integer, parameter :: WP = mytype ! inherit from decomp_2d, flag of -DDOUBLE_PREC is required.

end module precision_mod

!##############################################################################
module parameters_constant_mod
  use precision_mod
  implicit none

  real(WP), parameter :: ZPONE       = 0.1_WP
  real(WP), parameter :: ZPTWO       = 0.2_WP
  real(WP), parameter :: ZPTHREE     = 0.3_WP
  real(WP), parameter :: ZPFOUR      = 0.4_WP
  real(WP), parameter :: HALF        = 0.5_WP
  real(WP), parameter :: ZPSIX       = 0.6_WP
  real(WP), parameter :: ZPSEVEN     = 0.7_WP
  real(WP), parameter :: ZPEIGHT     = 0.8_WP
  real(WP), parameter :: ZPNINE      = 0.9_WP

  real(WP), parameter :: ZERO        = 0.0_WP
  real(WP), parameter :: ONE         = 1.0_WP
  real(WP), parameter :: ONEPFIVE    = 1.5_WP
  real(WP), parameter :: TWO         = 2.0_WP
  real(WP), parameter :: THREE       = 3.0_WP
  real(WP), parameter :: FOUR        = 4.0_WP
  real(WP), parameter :: FIVE        = 5.0_WP
  real(WP), parameter :: SIX         = 6.0_WP
  real(WP), parameter :: SEVEN       = 7.0_WP
  real(WP), parameter :: EIGHT       = 8.0_WP
  real(WP), parameter :: NINE        = 9.0_WP
  real(WP), parameter :: ONE_THIRD   = ONE / THREE
  real(WP), parameter :: TWO_THIRD   = TWO / THREE

  real(WP), parameter :: TEN         = 10.0_WP
  real(WP), parameter :: ELEVEN      = 11.0_WP
  real(WP), parameter :: TWELVE      = 12.0_WP
  real(WP), parameter :: THIRTEEN    = 13.0_WP
  real(WP), parameter :: FOURTEEN    = 14.0_WP
  real(WP), parameter :: FIFTEEN     = 15.0_WP
  real(WP), parameter :: SIXTEEN     = 16.0_WP
  real(WP), parameter :: SEVENTEEN   = 17.0_WP

  real(WP), parameter :: TWENTYTWO   = 22.0_WP
  real(WP), parameter :: TWENTYTHREE = 23.0_WP
  real(WP), parameter :: TWENTYFOUR  = 24.0_WP
  real(WP), parameter :: TWENTYFIVE  = 25.0_WP
  real(WP), parameter :: TWENTYSIX   = 26.0_WP
  real(WP), parameter :: TWENTYSEVEN = 27.0_WP

  real(WP), parameter :: THIRTYFIVE  = 35.0_WP
  real(WP), parameter :: THIRTYSIX   = 36.0_WP
  real(WP), parameter :: THIRTYSEVEN = 37.0_WP

  real(WP), parameter :: SIXTY       = 60.0_WP
  real(WP), parameter :: SIXTYTWO    = 62.0_WP
  real(WP), parameter :: SIXTYTHREE  = 63.0_WP

  real(WP), parameter :: EIGHTYSEVEN = 87.0_WP

  real(WP), parameter :: MINP        = 1.0E-20_WP
  real(WP), parameter :: MAXP        = 1.0E20_WP

  real(WP), parameter :: MINN        = -1.0E20_WP
  real(WP), parameter :: MAXN        = -1.0E-20_WP


  real(WP), parameter :: TRUNCERR    = 1.0E-16_WP

  

  real(WP), parameter :: PI          = dacos( -ONE )
  real(WP), parameter :: TWOPI       = TWO * dacos( -ONE )

  real(WP), parameter, dimension(3, 3) :: KRONECKER_DELTA = &
                                            reshape( (/ &
                                            ONE, ZERO, ZERO, &
                                            ZERO, ONE, ZERO, &
                                            ZERO, ZERO, ONE  /), &
                                            (/3, 3/) )

  real(WP), parameter :: GRAVITY     = 9.80665_WP

end module parameters_constant_mod
module udf_type_mod
  use precision_mod
  use mpi_mod
  implicit none
  type t_domain
    logical :: is_periodic(3)
    logical :: is_stretching(3)
    integer :: case
    integer :: np_geo(3) ! geometric points
    integer :: np(3) ! calculated points
    integer :: nc(3) ! geometric cell number
    integer :: bc(2, 3) ! (two sides, three directions)
    real(wp) :: ubc(2, 3) ! (two sides, three directions)
    real(wp) :: h(3) ! uniform dx
    real(wp) :: h1r(3) ! uniform (dx)^(-1)
    real(wp) :: h2r(3) ! uniform (dx)^(-2)

    integer(4) :: ux_xst(3) ! x-pencil, for variable ux and gx, starting index
    integer(4) :: ux_xen(3) ! x-pencil, for variable ux and gx, ending index
    integer(4) :: ux_xsz(3) ! x-pencil, for variable ux and gx, local size
    integer(4) :: uy_xst(3) ! x-pencil, for variable uy and gy, starting index
    integer(4) :: uy_xen(3) ! x-pencil, for variable uy and gy, ending index
    integer(4) :: uy_xsz(3) ! x-pencil, for variable uy and gy, local size
    integer(4) :: uz_xst(3) ! x-pencil, for variable uz and gz, starting index
    integer(4) :: uz_xen(3) ! x-pencil, for variable uz and gz, ending index
    integer(4) :: uz_xsz(3) ! x-pencil, for variable uz and gz, local size
    integer(4) :: ps_xst(3) ! x-pencil, for variable p and p, starting index
    integer(4) :: ps_xen(3) ! x-pencil, for variable p and p, ending index
    integer(4) :: ps_xsz(3) ! x-pencil, for variable p and p, local size

    integer(4) :: ux_yst(3) ! y-pencil, for variable ux and gx, starting index
    integer(4) :: ux_yen(3) ! y-pencil, for variable ux and gx, ending index
    integer(4) :: ux_ysz(3) ! y-pencil, for variable ux and gx, local size
    integer(4) :: uy_yst(3) ! y-pencil, for variable uy and gy, starting index
    integer(4) :: uy_yen(3) ! y-pencil, for variable uy and gy, ending index
    integer(4) :: uy_ysz(3) ! y-pencil, for variable uy and gy, local size
    integer(4) :: uz_yst(3) ! y-pencil, for variable uz and gz, starting index
    integer(4) :: uz_yen(3) ! y-pencil, for variable uz and gz, ending index
    integer(4) :: uz_ysz(3) ! y-pencil, for variable uz and gz, local size
    integer(4) :: ps_yst(3) ! y-pencil, for variable p and p, starting index
    integer(4) :: ps_yen(3) ! y-pencil, for variable p and p, ending index
    integer(4) :: ps_ysz(3) ! y-pencil, for variable p and p, local size

    integer(4) :: ux_zst(3) ! z-pencil, for variable ux and gx, starting index
    integer(4) :: ux_zen(3) ! z-pencil, for variable ux and gx, ending index
    integer(4) :: ux_zsz(3) ! z-pencil, for variable ux and gx, local size
    integer(4) :: uy_zst(3) ! z-pencil, for variable uy and gy, starting index
    integer(4) :: uy_zen(3) ! z-pencil, for variable uy and gy, ending index
    integer(4) :: uy_zsz(3) ! z-pencil, for variable uy and gy, local size
    integer(4) :: uz_zst(3) ! z-pencil, for variable uz and gz, starting index
    integer(4) :: uz_zen(3) ! z-pencil, for variable uz and gz, ending index
    integer(4) :: uz_zsz(3) ! z-pencil, for variable uz and gz, local size
    integer(4) :: ps_zst(3) ! z-pencil, for variable p and p, starting index
    integer(4) :: ps_zen(3) ! z-pencil, for variable p and p, ending index
    integer(4) :: ps_zsz(3) ! z-pencil, for variable p and p, local size

    type(DECOMP_INFO) :: dpcc ! eg, ux
    type(DECOMP_INFO) :: dcpc ! eg, uy
    type(DECOMP_INFO) :: dccp ! eg, uz
    type(DECOMP_INFO) :: dccc ! eg, p
    type(DECOMP_INFO) :: dppc ! eg, <ux>^y, <uy>^x
    type(DECOMP_INFO) :: dpcp ! eg, <ux>^z, <uz>^x
    type(DECOMP_INFO) :: dcpp ! eg, <uy>^z, <uz>^y


    integer(4), allocatable :: iNeighb(:, :) ! global index
    integer(4), allocatable :: jNeighb(:, :)
    integer(4), allocatable :: kNeighb(:, :)
    ! node location, mapping 
    real(wp), allocatable :: yMappingpt(:, :) ! j = 1, first coefficient in first deriviation. 1/h'
                                              ! j = 2, first coefficient in second deriviation 1/h'^2
                                              ! j = 3, second coefficient in second deriviation -h"/h'^3
    ! cell center location, mapping
    real(wp), allocatable :: yMappingcc(:, :) ! first coefficient in first deriviation. 1/h'
                                              ! first coefficient in second deriviation 1/h'^2
                                              ! second coefficient in second deriviation -h"/h'^3
    real(wp), allocatable :: yp(:)
    real(wp), allocatable :: yc(:)
  end type t_domain

  type t_flow
    real(WP) :: time
    real(WP) :: rre
    real(WP) :: fgravity                      ! 1 / Re
    real(WP), allocatable :: qx(:, :, :)  !
    real(WP), allocatable :: qy(:, :, :)
    real(WP), allocatable :: qz(:, :, :)
    real(WP), allocatable :: gx(:, :, :)
    real(WP), allocatable :: gy(:, :, :)
    real(WP), allocatable :: gz(:, :, :)

    real(WP), allocatable :: pres(:, :, :)
    real(WP), allocatable :: pcor(:, :, :)

    real(WP), allocatable :: dDens(:, :, :)
    real(WP), allocatable :: mVisc(:, :, :)
    real(WP), allocatable :: dDensm1(:, :, :)
    real(WP), allocatable :: dDensm2(:, :, :)

    real(WP), allocatable :: mx_rhs(:, :, :) ! current step rhs in x
    real(WP), allocatable :: my_rhs(:, :, :) ! current step rhs in y
    real(WP), allocatable :: mz_rhs(:, :, :) ! current step rhs in z

    real(WP), allocatable :: mx_rhs0(:, :, :)! last step rhs in x
    real(WP), allocatable :: my_rhs0(:, :, :)! last step rhs in y
    real(WP), allocatable :: mz_rhs0(:, :, :)! last step rhs in z

  end type t_flow

  type t_thermo
    real(WP) :: time
    real(WP) :: rPrRen
    real(WP), allocatable :: dh(:, :, :)
    real(WP), allocatable :: hEnth(:, :, :)
    real(WP), allocatable :: kCond(:, :, :)
    real(WP), allocatable :: tTemp(:, :, :)
  end type t_thermo

end module

module type_vars_mod
  use udf_type_mod
  implicit none

  type(t_domain), save :: domain
  type(t_flow),   save :: flow
  type(t_thermo), save :: thermo
end module

!##############################################################################
module math_mod
  use precision_mod
  use parameters_constant_mod, only : ONE, ZERO, MINP
  implicit none

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
  elemental function abs_sp ( r ) result(d)
  real(kind = S6P), intent(in) :: r
  real(kind = S6P) :: d
    d = abs ( r )
  end function

  elemental function abs_dp ( r ) result (d)
  real(kind = D15P), intent(in) :: r
  real(kind = D15P) :: d
    d = dabs ( r ) 
  end function

  ! sqrt
  pure function sqrt_sp ( r ) result(d)
    real(kind = S6P), intent(in) :: r
    real(kind = S6P) :: d
    d = sqrt ( r )
  end function

  pure function sqrt_dp ( r ) result (d)
    real(kind = D15P), intent(in) :: r
    real(kind = D15P) :: d
    d = dsqrt ( r ) 
  end function

  ! sin
  pure function sin_sp ( r ) result(d)
    real(kind = S6P), intent(in) :: r
    real(kind = S6P) :: d
    d = sin ( r )
  end function

  pure function sin_dp ( r ) result (d)
    real(kind = D15P), intent(in) :: r
    real(kind = D15P) :: d
    d = dsin ( r ) 
  end function

  ! cos
  pure function cos_sp ( r ) result(d)
    real(kind = S6P), intent(in) :: r
    real(kind = S6P) :: d
    d = cos ( r )
  end function

  pure function cos_dp ( r ) result (d)
    real(kind = D15P), intent(in) :: r
    real(kind = D15P) :: d
    d = dcos ( r ) 
  end function

  ! tanh
  pure function tanh_sp ( r ) result(d)
    real(kind = S6P), intent(in) :: r
    real(kind = S6P) :: d
    d = tanh ( r )
  end function

  pure function tanh_dp ( r ) result (d)
    real(kind = D15P), intent(in) :: r
    real(kind = D15P) :: d
    d = dtanh ( r ) 
  end function

  ! tan
  pure function tan_sp ( r ) result(d)
    real(kind = S6P), intent(in) :: r
    real(kind = S6P) :: d
    d = tan ( r )
  end function

  pure function tan_dp ( r ) result (d)
    real(kind = D15P), intent(in) :: r
    real(kind = D15P) :: d
    d = tan ( r ) 
  end function

  ! atan
  pure function atan_sp ( r ) result(d)
    real(kind = S6P), intent(in) :: r
    real(kind = S6P) :: d
    d = atan ( r )
  end function

  pure function atan_dp ( r ) result (d)
    real(kind = D15P), intent(in) :: r
    real(kind = D15P) :: d
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

module typeconvert_mod
contains
  character(len=20) function int2str(k)
    implicit none
    integer, intent(in) :: k
    write (int2str, *) k
    int2str = adjustl(int2str)
  end function int2str
  character(len=20) function real2str(r)
    use precision_mod
    implicit none
    real(wp), intent(in) :: r
    write (real2str, '(F10.4)') r
    real2str = adjustl(real2str)
  end function real2str
end module typeconvert_mod


