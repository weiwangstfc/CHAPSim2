!=============================================================================================================================================
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
!=============================================================================================================================================
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
  real(WP), parameter :: twopfive    = 2.5_WP
  real(WP), parameter :: THREE       = 3.0_WP
  real(WP), parameter :: threepfive  = 3.5_WP
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

  real(WP), parameter :: THIRTYTWO   = 32.0_WP
  real(WP), parameter :: THIRTYFIVE  = 35.0_WP
  real(WP), parameter :: THIRTYSIX   = 36.0_WP
  real(WP), parameter :: THIRTYSEVEN = 37.0_WP

  real(WP), parameter :: FOURTYFIVE  = 45.0_WP

  real(WP), parameter :: FIFTY       = 50.0_WP

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

  complex(mytype),parameter :: cx_one_one=cmplx(one, one, kind=mytype)

  real(WP), parameter, dimension(3, 3) :: KRONECKER_DELTA = &
                                            reshape( (/ &
                                            ONE, ZERO, ZERO, &
                                            ZERO, ONE, ZERO, &
                                            ZERO, ZERO, ONE  /), &
                                            (/3, 3/) )

  real(WP), parameter :: GRAVITY     = 9.80665_WP

  integer, parameter :: nvd = 3

  integer, parameter :: ICASE_CHANNEL = 1, &
                        ICASE_PIPE    = 2, &
                        ICASE_ANNUAL  = 3, &
                        ICASE_TGV3D   = 4, &
                        ICASE_TGV2D   = 5, &
                        ICASE_BURGERS = 6

  integer, parameter :: ICARTESIAN   = 1, &
                        ICYLINDRICAL = 2
                        
  integer, parameter :: ISTRET_NO     = 0, &
                        ISTRET_CENTRE = 1, &
                        ISTRET_2SIDES = 2, &
                        ISTRET_BOTTOM = 3, &
                        ISTRET_TOP    = 4
                        

  integer, parameter :: ITIME_RK3    = 3, &
                        ITIME_RK3_CN = 2, &
                        ITIME_AB2    = 1

  ! warning : Don't change below order for BC types.
  integer, parameter :: IBC_INTERIOR    = 0, &
                        IBC_PERIODIC    = 1, &
                        IBC_SYMMETRIC   = 2, &
                        IBC_ASYMMETRIC  = 3, &
                        IBC_DIRICHLET   = 4, &
                        IBC_NEUMANN     = 5, &
                        IBC_INTRPL      = 6, &
                        IBC_CONVECTIVE  = 7, &
                        IBC_TURBGEN     = 8, &
                        IBC_DATABASE    = 9
!                        IBC_INLET_MEAN  = 4, &
!                        IBC_INLET_TG    = 5, &
!                        IBC_INLET_MAP   = 6, &
!                        IBC_INLET_DB    = 7, &
!                        IBC_OUTLET_EXPO = 8, &
!                        IBC_OUTLET_CONV = 9, &
!                        IBC_INTERIOR    = 0, &
                        
  integer, parameter :: IACCU_CD2 = 2, &
                        IACCU_CD4 = 3, &
                        IACCU_CP4 = 4, &
                        IACCU_CP6 = 6

  integer, parameter :: NDIM = 3

  integer, parameter :: INITIAL_RANDOM  = 0, &
                        INITIAL_RESTART = 1, &
                        INITIAL_INTERPL = 2

  integer, parameter :: IVIS_EXPLICIT   = 1, &
                        IVIS_SEMIMPLT   = 2

  integer, parameter :: IDRVF_NO        = 0, &
                        IDRVF_MASSFLUX  = 1, &
                        IDRVF_SKINFRIC  = 2, &
                        IDRVF_PRESLOSS  = 3

  integer, parameter :: THERMAL_BC_CONST_T  = 0, &
                        THERMAL_BC_CONST_H  = 1

  integer, parameter :: ISCP_WATER      = 1, &
                        ISCP_CO2        = 2, &
                        ILIQUID_SODIUM  = 3, &
                        ILIQUID_LEAD    = 4, &
                        ILIQUID_BISMUTH = 5, &
                        ILIQUID_LBE     = 6, &
                        ILIQUID_WATER   = 7 ! to be updated 

  integer, parameter :: IPROPERTY_TABLE = 1, &
                        IPROPERTY_FUNCS = 2

  character(len = 64), parameter :: INPUT_SCP_WATER = 'NIST_WATER_23.5MP.DAT'
  character(len = 64), parameter :: INPUT_SCP_CO2   = 'NIST_CO2_8MP.DAT'

  real(WP), parameter :: TM0_Na  = 371.0  ! unit: K, melting temperature at 1 atm for Na
  real(WP), parameter :: TM0_Pb  = 600.6  ! unit: K, melting temperature at 1 atm for Lead
  real(WP), parameter :: TM0_BI  = 544.6  ! unit: K, melting temperature at 1 atm for Bismuth
  real(WP), parameter :: TM0_LBE = 398.0  ! unit: K, melting temperature at 1 atm for LBE
  real(WP), parameter :: TM0_H2O = 273.15 ! unit: K, melting temperature at 1 atm for water

  real(WP), parameter :: TB0_Na  = 1155.0 ! unit: K, boling temperature at 1 atm for Na
  real(WP), parameter :: TB0_Pb  = 2021.0 ! unit: K, boling temperature at 1 atm for Lead
  real(WP), parameter :: TB0_BI  = 1831.0 ! unit: K, boling temperature at 1 atm for Bismuth
  real(WP), parameter :: TB0_LBE = 1927.0 ! unit: K, boling temperature at 1 atm for LBE
  real(WP), parameter :: TB0_H2O = 373.15 ! unit: K, boling temperature at 1 atm for water

  real(WP), parameter :: HM0_Na  = 113.0e3 ! unit: J / Kg, latent melting heat, enthalpy for Na
  real(WP), parameter :: HM0_Pb  = 23.07e3 ! unit: J / Kg, latent melting heat, enthalpy for Lead
  real(WP), parameter :: HM0_BI  = 53.3e3  ! unit: J / Kg, latent melting heat, enthalpy for Bismuth
  real(WP), parameter :: HM0_LBE = 38.6e3  ! unit: J / Kg, latent melting heat, enthalpy for LBE
  real(WP), parameter :: HM0_H2O = 334.0e3 ! unit: J / Kg, latent melting heat, enthalpy for water

  ! D = CoD(0) + CoD(1) * T
  real(WP), parameter :: CoD_Na(0:1) = (/1014.0, -0.235/)
  real(WP), parameter :: CoD_Pb(0:1) = (/11441.0, -1.2795/)
  real(WP), parameter :: CoD_Bi(0:1) = (/10725.0, -1.22 /)
  real(WP), parameter :: CoD_LBE(0:1) = (/11065.0, 1.293 /)

  ! K = CoK(0) + CoK(1) * T + CoK(2) * T^2
  real(WP), parameter :: CoK_Na(0:2) = (/104.0, -0.047, 0.0/)
  real(WP), parameter :: CoK_Pb(0:2) = (/9.2, 0.011, 0.0/)
  real(WP), parameter :: CoK_Bi(0:2) = (/7.34, 9.5E-3, 0.0/)
  real(WP), parameter :: CoK_LBE(0:2) = (/ 3.284, 1.617E-2, -2.305E-6/)

  ! B = 1 / (CoB - T)
  real(WP), parameter :: CoB_Na = 4316.0
  real(WP), parameter :: CoB_Pb = 8942.0
  real(WP), parameter :: CoB_BI = 8791.0
  real(WP), parameter :: CoB_LBE = 8558.0

  ! Cp = CoCp(-2) * T^(-2) + CoCp(-1) * T^(-1) + CoCp(0) + CoCp(1) * T + CoCp(2) * T^2
  real(WP), parameter :: CoCp_Na(-2:2) = (/- 3.001e6, 0.0, 1658.0, -0.8479, 4.454E-4/)
  real(WP), parameter :: CoCp_Pb(-2:2) = (/- 1.524e6, 0.0, 176.2, -4.923E-2, 1.544E-5/)
  real(WP), parameter :: CoCp_Bi(-2:2) = (/7.183e6, 0.0, 118.2, 5.934E-3, 0.0/)
  real(WP), parameter :: CoCp_LBE(-2:2) = (/-4.56e5, 0.0, 164.8, - 3.94E-2, 1.25E-5/)

  ! H = HM0 + CoH(-1) * (1 / T - 1 / Tm0) + CoH(0) + CoH(1) * (T - Tm0) +  CoH(2) * (T^2 - Tm0^2) +  CoH(3) * (T^3- Tm0^3)
  real(WP), parameter :: CoH_Na(-1:3) = (/4.56e5, 0.0, 164.8, -1.97E-2, 4.167E-4/)
  real(WP), parameter :: CoH_Pb(-1:3) = (/1.524e6, 0.0, 176.2, -2.4615E-2, 5.147E-6/)
  real(WP), parameter :: CoH_Bi(-1:3) = (/-7.183e6, 0.0, 118.2, 2.967E-3, 0.0/)
  real(WP), parameter :: CoH_LBE(-1:3) = (/4.56e5, 0.0, 164.8, -1.97E-2, 4.167E-4/)! check, WRong from literature.

  ! M = vARies
  real(WP), parameter :: CoM_Na(-1:1) = (/556.835, -6.4406, -0.3958/) ! M = exp ( CoM(-1) / T + CoM(0) + CoM(1) * ln(T) )
  real(WP), parameter :: CoM_Pb(-1:1) = (/1069.0, 4.55E-4, 0.0/) ! M = CoM(0) * exp (CoM(-1) / T)
  real(WP), parameter :: CoM_Bi(-1:1) = (/780.0, 4.456E-4, 0.0/) ! M = CoM(0) * exp (CoM(-1) / T)
  real(WP), parameter :: CoM_LBE(-1:1) = (/754.1, 4.94E-4, 0.0/) ! M = CoM(0) * exp (CoM(-1) / T)
end module parameters_constant_mod
!=============================================================================================================================================
module wtformat_mod
  use iso_fortran_env
  implicit none

  character(len = 17) :: wrtfmt1i   = '(2X, A48, 1I20.1)'
  character(len = 17) :: wrtfmt2i   = '(2X, A48, 2I10.1)'
  character(len = 17) :: wrtfmt1r   = '(2X, A48, 1F20.4)'
  character(len = 17) :: wrtfmt2r   = '(2X, A48, 2F10.2)'
  character(len = 17) :: wrtfmt1e   = '(2X, A48, 1E20.4)'
  character(len = 25) :: wrtfmt1i1r = '(2X, A48, 1I10.1, 1F10.2)'
  character(len = 25) :: wrtfmt2i2r = '(2X, A48, 2I10.1, 2F10.2)'
  character(len = 14) :: wrtfmt1s   = '(2X, A48, A20)'

end module wtformat_mod
!=============================================================================================================================================
module udf_type_mod
  use precision_mod
  use mpi_mod
  implicit none
!---------------------------------------------------------------------------------------------------------------------------------------------
!  domain info
!--------------------------------------------------------------------------------------------------------------------------------------------- 
  type t_domain
    logical :: is_periodic(3)
    logical :: is_stretching(3)
    integer :: idom
    integer :: icase
    integer :: icoordinate
    integer :: ithermo
    integer :: icht
    integer :: iTimeScheme
    integer :: iviscous
    integer :: iAccuracy
    integer :: nfreqckpt
    integer :: nvisu
    integer :: nIterStatsStart
    integer :: nfreqStats
    integer :: nsubitr
    integer :: istret
    integer :: nc(3) ! geometric cell number
    integer :: np_geo(3) ! geometric points
    integer :: np(3) ! calculated points
    integer  :: ibcx(2, 5) ! bc type, (5 variables, 2 sides), u, v, w, p, T
    integer  :: ibcy(2, 5) ! bc type, (5 variables, 2 sides)
    integer  :: ibcz(2, 5) ! bc type, (5 variables, 2 sides)
    real(wp) :: fbcx(2, 5) ! bc values, (5 variables, 2 sides)
    real(wp) :: fbcy(2, 5) ! bc values, (5 variables, 2 sides)
    real(wp) :: fbcz(2, 5) ! bc values, (5 variables, 2 sides)
    real(WP) :: fbc_vism(2, 3) ! bc values for mu, in 3 direction, 2 sides.
    real(WP) :: fbc_dend(2, 3) ! bc values for density, in 3 direction, 2 sides.
    real(wp) :: lxx
    real(wp) :: lyt
    real(wp) :: lyb
    real(wp) :: lzz
    real(WP) :: rstret
    real(wp) :: dt

    real(wp) :: h(3) ! uniform dx
    real(wp) :: h1r(3) ! uniform (dx)^(-1)
    real(wp) :: h2r(3) ! uniform (dx)^(-2)
    real(wp) :: tGamma(0:3)
    real(wp) :: tZeta (0:3)
    real(wp) :: tAlpha(0:3)
    real(wp) :: sigma1p, sigma2p

    type(DECOMP_INFO) :: dccc ! eg, p
    type(DECOMP_INFO) :: dpcc ! eg, ux
    type(DECOMP_INFO) :: dcpc ! eg, uy
    type(DECOMP_INFO) :: dccp ! eg, uz
    type(DECOMP_INFO) :: dppc ! eg, <ux>^y, <uy>^x
    type(DECOMP_INFO) :: dpcp ! eg, <ux>^z, <uz>^x
    type(DECOMP_INFO) :: dcpp ! eg, <uy>^z, <uz>^y
    type(DECOMP_INFO) :: dppp

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
!---------------------------------------------------------------------------------------------------------------------------------------------
!  flow info
!--------------------------------------------------------------------------------------------------------------------------------------------- 
  type t_flow
    integer  :: idriven
    integer  :: igravity
    integer  :: irestart
    integer  :: nrsttckpt
    integer  :: nIterIniRen
    integer  :: nIterFlowStart
    integer  :: nIterFlowEnd
    integer  :: iteration

    real(WP) :: time
    real(WP) :: ren
    real(WP) :: rre
    real(WP) :: drvfc
    real(WP) :: fgravity(3)
    real(wp) :: renIni
    real(wp) :: initNoise
  
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

  
  type t_fluidThermoProperty
    real(WP) :: t  ! temperature
    real(WP) :: d  ! density
    real(WP) :: m  ! dynviscosity
    real(WP) :: k  ! thermconductivity
    real(WP) :: h  ! enthalpy
    real(WP) :: dh ! mass enthalpy
    real(WP) :: cp ! specific heat capacity 
    real(WP) :: b  ! thermal expansion
  end type t_fluidThermoProperty

  type t_thermo
    integer :: ifluid
    integer  :: irestart
    integer  :: nrsttckpt
    integer  :: iteration
    integer  :: nIterThermoStart
    integer  :: nIterThermoEnd
    real(WP) :: lenRef
    real(WP) :: t0ref ! '0' means dimensional 
    real(WP) :: t0ini
    real(WP) :: time
    real(WP) :: rPrRen
    
    type(t_fluidThermoProperty) :: ftpbcx(2)  ! undim, xbc state
    type(t_fluidThermoProperty) :: ftpbcy(2)  ! undim, ybc state
    type(t_fluidThermoProperty) :: ftpbcz(2)  ! undim, zbc state

    real(WP), allocatable :: dh(:, :, :)
    real(WP), allocatable :: hEnth(:, :, :)
    real(WP), allocatable :: kCond(:, :, :)
    real(WP), allocatable :: tTemp(:, :, :)
    real(WP), allocatable :: ene_rhs(:, :, :)  ! current step rhs
    real(WP), allocatable :: ene_rhs0(:, :, :) ! last step rhs
  end type t_thermo

  type t_fluid_parameter
    character(len = 64) :: inputProperty
    integer :: ifluid
    integer :: ipropertyState
    integer :: nlist
    real(WP) :: TM0
    real(WP) :: TB0
    real(WP) :: HM0
    real(WP) :: CoD(0:1)
    real(WP) :: CoK(0:2)
    real(WP) :: CoB
    real(WP) :: CoCp(-2:2)
    real(WP) :: CoH(-1:3)
    real(WP) :: CoM(-1:1)
    type(t_fluidThermoProperty) :: ftp0ref    ! dim, reference state
    type(t_fluidThermoProperty) :: ftpini     ! undim, initial state
  end type t_fluid_parameter



end module
!=============================================================================================================================================
!=============================================================================================================================================
module vars_df_mod
  use udf_type_mod
  implicit none

  type(t_domain), allocatable, save :: domain(:)
  type(t_flow),   allocatable, save :: flow(:)
  type(t_thermo), allocatable, save :: thermo(:)
end module
!=============================================================================================================================================
!=============================================================================================================================================
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

  interface abs_prec
    module procedure abs_sp
    module procedure abs_dp
    module procedure abs_csp
    module procedure abs_cdp
  end interface abs_prec

  interface sin_wp
    module procedure sin_sp
    module procedure sin_dp
  end interface sin_wp

  interface sin_prec
    module procedure sin_sp
    module procedure sin_dp
  end interface sin_prec

  interface cos_wp
    module procedure cos_sp
    module procedure cos_dp
  end interface cos_wp

  interface cos_prec
    module procedure cos_sp
    module procedure cos_dp
  end interface cos_prec

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

  elemental function abs_csp ( r ) result(d)
  COMPLEX(kind = S6P), intent(in) :: r
  real(kind = S6P) :: d
    d = cabs ( r )
  end function

  elemental function abs_cdp ( r ) result (d)
  COMPLEX(kind = D15P), intent(in) :: r
  real(kind = D15P) :: d
    d = cdabs ( r ) 
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
!=============================================================================================================================================
!=============================================================================================================================================
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


