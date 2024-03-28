!==========================================================================================================
module precision_mod
  use mpi_mod
  implicit none

  public
  integer, parameter :: I4 = selected_int_kind( 4 )
  integer, parameter :: I8 = selected_int_kind( 8 )
  integer, parameter :: I15 = selected_int_kind( 15 )
  integer, parameter :: S6P = selected_real_kind( p = 6, r = 37 )
  integer, parameter :: D15P = selected_real_kind( p = 15, r = 307 )
  integer, parameter :: Q33P = selected_real_kind( p = 33, r = 4931 )
#ifdef DOUBLE_PREC
  integer, parameter :: WP = D15P
  integer, parameter :: MPI_REAL_WP = MPI_DOUBLE_PRECISION
  integer, parameter :: MPI_CPLX_WP = MPI_DOUBLE_COMPLEX
#else
  integer, parameter :: WP = S6P !D15P
  integer, parameter :: MPI_REAL_WP = MPI_REAL
  integer, parameter :: MPI_CPLX_WP = MPI_COMPLEX
#endif

end module precision_mod
!==========================================================================================================
module parameters_constant_mod
  use precision_mod
  implicit none
!----------------------------------------------------------------------------------------------------------
! constants
!----------------------------------------------------------------------------------------------------------
  real(WP), parameter :: ZPONE       = 0.1_WP
  real(WP), parameter :: EIGHTH      = 0.125_WP
  real(WP), parameter :: ZPTWO       = 0.2_WP
  real(WP), parameter :: QUARTER     = 0.25_WP
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
  real(WP), parameter :: TWOPFIVE    = 2.5_WP
  real(WP), parameter :: THREE       = 3.0_WP
  real(WP), parameter :: threepfive  = 3.5_WP
  real(WP), parameter :: FOUR        = 4.0_WP
  real(WP), parameter :: FIVE        = 5.0_WP
  real(WP), parameter :: SIX         = 6.0_WP
  real(WP), parameter :: SEVEN       = 7.0_WP
  real(WP), parameter :: EIGHT       = 8.0_WP
  real(WP), parameter :: NINE        = 9.0_WP
  real(WP), parameter :: ONE_THIRD   = 0.33333333333333333333_WP
  real(WP), parameter :: TWO_THIRD   = 0.66666666666666666667_WP

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
  real(WP), parameter :: MAXVELO     = 1.0E3_WP

  real(WP), parameter :: MINN        = -1.0E20_WP
  real(WP), parameter :: MAXN        = -1.0E-20_WP


  real(WP), parameter :: TRUNCERR    = 1.0E-16_WP

  

  real(WP), parameter :: PI          = 2.0_WP*(DASIN(1.0_WP)) !3.14159265358979323846_WP !dacos( -ONE ) 
  real(WP), parameter :: TWOPI       = TWO * PI !6.28318530717958647692_WP!TWO * dacos( -ONE )

  complex(mytype),parameter :: cx_one_one=cmplx(one, one, kind=mytype)

  real(WP), parameter, dimension(3, 3) :: KRONECKER_DELTA = &
                                            reshape( (/ &
                                            ONE, ZERO, ZERO, &
                                            ZERO, ONE, ZERO, &
                                            ZERO, ZERO, ONE  /), &
                                            (/3, 3/) )

  real(WP), parameter :: GRAVITY     = 9.80665_WP
!----------------------------------------------------------------------------------------------------------
! case id
!----------------------------------------------------------------------------------------------------------
  integer, parameter :: ICASE_OTHERS = 0, &
                        ICASE_CHANNEL = 1, &
                        ICASE_PIPE    = 2, &
                        ICASE_ANNUAL  = 3, &
                        ICASE_TGV3D   = 4, &
                        ICASE_TGV2D   = 5, &
                        ICASE_BURGERS = 6, &
                        ICASE_ALGTEST = 7
  integer, parameter :: NDIM = 3
!----------------------------------------------------------------------------------------------------------
! flow initilisation
!----------------------------------------------------------------------------------------------------------     
  integer, parameter :: INIT_RESTART = 0, &
                        INIT_INTERPL = 1, &
                        INIT_RANDOM  = 2, &
                        INIT_INLET   = 3, &
                        INIT_GVCONST = 4, &
                        INIT_POISEUILLE = 5, &
                        INIT_FUNCTION = 6
!----------------------------------------------------------------------------------------------------------
! coordinates
!----------------------------------------------------------------------------------------------------------
  integer, parameter :: ICARTESIAN   = 1, &
                        ICYLINDRICAL = 2
!----------------------------------------------------------------------------------------------------------
! grid stretching
!----------------------------------------------------------------------------------------------------------               
  integer, parameter :: ISTRET_NO     = 0, &
                        ISTRET_CENTRE = 1, &
                        ISTRET_2SIDES = 2, &
                        ISTRET_BOTTOM = 3, &
                        ISTRET_TOP    = 4               
!----------------------------------------------------------------------------------------------------------
! time scheme
!----------------------------------------------------------------------------------------------------------
  integer, parameter :: ITIME_RK3    = 3, &
                        ITIME_RK3_CN = 2, &
                        ITIME_AB2    = 1
!----------------------------------------------------------------------------------------------------------
! BC
!----------------------------------------------------------------------------------------------------------
  ! warning : Don't change below order for BC types.
  integer, parameter :: IBC_INTERIOR    = 0, & ! basic and nominal, used in operations, bulk, 2 ghost layers
                        IBC_PERIODIC    = 1, & ! basic and nominal, used in operations 
                        IBC_SYMMETRIC   = 2, & ! basic and nominal, used in operations
                        IBC_ASYMMETRIC  = 3, & ! basic and nominal, used in operations
                        IBC_DIRICHLET   = 4, & ! basic and nominal, used in operations
                        IBC_NEUMANN     = 5, & ! basic and nominal, used in operations
                        IBC_INTRPL      = 6, & ! basic only, for all others, used in operations
                        IBC_CONVECTIVE  = 7, & ! nominal only, = IBC_INTPRL
                        IBC_TURBGEN     = 8, & ! nominal only, = IBC_PERIODIC, bulk, 2 ghost layers
                        IBC_PROFILE1D   = 9, & ! nominal only, = IBC_DIRICHLET
                        IBC_DATABASE    = 10, &! nominal only, = IBC_PERIODIC, bulk, 2 ghost layers 
                        IBC_OTHERS      = 11   ! exclusive
  integer, parameter :: NBC = 5! u, v, w, p, T
!----------------------------------------------------------------------------------------------------------
! numerical accuracy
!----------------------------------------------------------------------------------------------------------             
  integer, parameter :: IACCU_CD2 = 2, &
                        IACCU_CD4 = 3, &
                        IACCU_CP4 = 4, &
                        IACCU_CP6 = 6
!----------------------------------------------------------------------------------------------------------
! numerical scheme for viscous term
!---------------------------------------------------------------------------------------------------------- 
  integer, parameter :: IVIS_EXPLICIT   = 1, &
                        IVIS_SEMIMPLT   = 2
!----------------------------------------------------------------------------------------------------------
! driven force in periodic flow
!---------------------------------------------------------------------------------------------------------- 
  integer, parameter :: IDRVF_NO         = 0, &
                        IDRVF_X_MASSFLUX = 1, &
                        IDRVF_X_Cf       = 2, &
                        IDRVF_Z_MASSFLUX = 3, &
                        IDRVF_Z_Cf       = 4 
!----------------------------------------------------------------------------------------------------------
! BC for thermal
!---------------------------------------------------------------------------------------------------------- 
  integer, parameter :: THERMAL_BC_CONST_T  = 0, &
                        THERMAL_BC_CONST_H  = 1
!----------------------------------------------------------------------------------------------------------
! working fluid media
!---------------------------------------------------------------------------------------------------------- 
  integer, parameter :: ISCP_WATER      = 1, &
                        ISCP_CO2        = 2, &
                        ILIQUID_SODIUM  = 3, &
                        ILIQUID_LEAD    = 4, &
                        ILIQUID_BISMUTH = 5, &
                        ILIQUID_LBE     = 6, &
                        ILIQUID_WATER   = 7 ! to be updated 
!----------------------------------------------------------------------------------------------------------
! physical property
!---------------------------------------------------------------------------------------------------------- 
  integer, parameter :: IPROPERTY_TABLE = 1, &
                        IPROPERTY_FUNCS = 2
!----------------------------------------------------------------------------------------------------------
! database for physical property
!----------------------------------------------------------------------------------------------------------
  character(len = 64), parameter :: INPUT_SCP_WATER = 'NIST_WATER_23.5MP.DAT'
  character(len = 64), parameter :: INPUT_SCP_CO2   = 'NIST_CO2_8MP.DAT'

  real(WP), parameter :: TM0_Na  = 371.0_WP  ! unit: K, melting temperature at 1 atm for Na
  real(WP), parameter :: TM0_Pb  = 600.6_WP  ! unit: K, melting temperature at 1 atm for Lead
  real(WP), parameter :: TM0_BI  = 544.6_WP  ! unit: K, melting temperature at 1 atm for Bismuth
  real(WP), parameter :: TM0_LBE = 398.0_WP  ! unit: K, melting temperature at 1 atm for LBE
  real(WP), parameter :: TM0_H2O = 273.15_WP ! unit: K, melting temperature at 1 atm for water

  real(WP), parameter :: TB0_Na  = 1155.0_WP ! unit: K, boling temperature at 1 atm for Na
  real(WP), parameter :: TB0_Pb  = 2021.0_WP ! unit: K, boling temperature at 1 atm for Lead
  real(WP), parameter :: TB0_BI  = 1831.0_WP ! unit: K, boling temperature at 1 atm for Bismuth
  real(WP), parameter :: TB0_LBE = 1927.0_WP ! unit: K, boling temperature at 1 atm for LBE
  real(WP), parameter :: TB0_H2O = 373.15_WP ! unit: K, boling temperature at 1 atm for water

  real(WP), parameter :: HM0_Na  = 113.0e3_WP ! unit: J / Kg, latent melting heat, enthalpy for Na
  real(WP), parameter :: HM0_Pb  = 23.07e3_WP ! unit: J / Kg, latent melting heat, enthalpy for Lead
  real(WP), parameter :: HM0_BI  =  53.3e3_WP ! unit: J / Kg, latent melting heat, enthalpy for Bismuth
  real(WP), parameter :: HM0_LBE =  38.6e3_WP ! unit: J / Kg, latent melting heat, enthalpy for LBE
  real(WP), parameter :: HM0_H2O = 334.0e3_WP ! unit: J / Kg, latent melting heat, enthalpy for water

  ! D = CoD(0) + CoD(1) * T
  real(WP), parameter :: CoD_Na(0:1)  = (/ 1014.0_WP,  -0.235_WP /)
  real(WP), parameter :: CoD_Pb(0:1)  = (/11441.0_WP, -1.2795_WP /)
  real(WP), parameter :: CoD_Bi(0:1)  = (/10725.0_WP,   -1.22_WP /)
  real(WP), parameter :: CoD_LBE(0:1) = (/11065.0_WP,   1.293_WP /)

  ! K = CoK(0) + CoK(1) * T + CoK(2) * T^2
  real(WP), parameter :: CoK_Na(0:2)  = (/104.0_WP,   -0.047_WP,       0.0_WP/)
  real(WP), parameter :: CoK_Pb(0:2)  = (/  9.2_WP,    0.011_WP,       0.0_WP/)
  real(WP), parameter :: CoK_Bi(0:2)  = (/ 7.34_WP,   9.5E-3_WP,       0.0_WP/)
  real(WP), parameter :: CoK_LBE(0:2) = (/3.284_WP, 1.617E-2_WP, -2.305E-6_WP/)

  ! B = 1 / (CoB - T)
  real(WP), parameter :: CoB_Na = 4316.0_WP
  real(WP), parameter :: CoB_Pb = 8942.0_WP
  real(WP), parameter :: CoB_BI = 8791.0_WP
  real(WP), parameter :: CoB_LBE= 8558.0_WP

  ! Cp = CoCp(-2) * T^(-2) + CoCp(-1) * T^(-1) + CoCp(0) + CoCp(1) * T + CoCp(2) * T^2
  real(WP), parameter :: CoCp_Na(-2:2) = (/-3.001e6_WP, 0.0_WP, 1658.0_WP,   -0.8479_WP, 4.454E-4_WP/)
  real(WP), parameter :: CoCp_Pb(-2:2) = (/-1.524e6_WP, 0.0_WP,  176.2_WP, -4.923E-2_WP, 1.544E-5_WP/)
  real(WP), parameter :: CoCp_Bi(-2:2) = (/ 7.183e6_WP, 0.0_WP,  118.2_WP,  5.934E-3_WP,      0.0_WP/)
  real(WP), parameter :: CoCp_LBE(-2:2)= (/-4.56e5_WP, 0.0_WP,  164.8_WP, - 3.94E-2_WP,  1.25E-5_WP/)

  ! H = HM0 + CoH(-1) * (1 / T - 1 / Tm0) + CoH(0) + CoH(1) * (T - Tm0) +  CoH(2) * (T^2 - Tm0^2) +  CoH(3) * (T^3- Tm0^3)
  real(WP), parameter :: CoH_Na(-1:3)  = (/  4.56e5_WP, 0.0_WP, 164.8_WP,   -1.97E-2_WP, 4.167E-4_WP/)
  real(WP), parameter :: CoH_Pb(-1:3)  = (/ 1.524e6_WP, 0.0_WP, 176.2_WP, -2.4615E-2_WP, 5.147E-6_WP/)
  real(WP), parameter :: CoH_Bi(-1:3)  = (/-7.183e6_WP, 0.0_WP, 118.2_WP,   2.967E-3_WP,      0.0_WP/)
  real(WP), parameter :: CoH_LBE(-1:3) = (/  4.56e5_WP, 0.0_WP, 164.8_WP,   -1.97E-2_WP, 4.167E-4_WP/)! check, WRong from literature.

  ! M = vARies
  real(WP), parameter :: CoM_Na(-1:1) = (/556.835_WP,  -6.4406_WP, -0.3958_WP/) ! M = exp ( CoM(-1) / T + CoM(0) + CoM(1) * ln(T) )
  real(WP), parameter :: CoM_Pb(-1:1) = (/ 1069.0_WP,  4.55E-4_WP,     0.0_WP/) ! M = CoM(0) * exp (CoM(-1) / T)
  real(WP), parameter :: CoM_Bi(-1:1) = (/  780.0_WP, 4.456E-4_WP,     0.0_WP/) ! M = CoM(0) * exp (CoM(-1) / T)
  real(WP), parameter :: CoM_LBE(-1:1)= (/  754.1_WP,  4.94E-4_WP,     0.0_WP/) ! M = CoM(0) * exp (CoM(-1) / T)
end module parameters_constant_mod
!==========================================================================================================
module wtformat_mod
  !use iso_fortran_env
  implicit none

  character(len = 17) :: wrtfmt1i   = '(2X, A48, 1I20.1)'
  character(len = 17) :: wrtfmt2i   = '(2X, A48, 2I10.1)'
  character(len = 17) :: wrtfmt3i   = '(2X, A48, 3I10.1)'
  character(len = 17) :: wrtfmt4i   = '(2X, A48, 4I10.1)'
  character(len = 17) :: wrtfmt1r   = '(2X, A48, 1F20.6)'
  character(len = 17) :: wrtfmt2r   = '(2X, A48, 2F10.6)'
  character(len = 18) :: wrtfmt3r   = '(2X, A48, 3F23.15)'
  character(len = 19) :: wrtfmt1e   = '(2X, A48, 1ES23.15)'
  character(len = 34) :: wrtfmt2e   = '(2X, A24, 1ES23.15, A24, 1ES23.15)'
  character(len = 25) :: wrtfmt1i1r = '(2X, A48, 1I10.1, 1F10.6)'
  character(len = 25) :: wrtfmt2i2r = '(2X, A48, 2I10.1, 2F10.6)'
  character(len = 14) :: wrtfmt3l   = '(2X, A48, 3L3)'
  character(len = 3)  :: wrtfmt1s   = '(A)'

end module wtformat_mod
!==========================================================================================================
module udf_type_mod
  use parameters_constant_mod, only: NDIM, NBC, WP
  use mpi_mod
  implicit none
!----------------------------------------------------------------------------------------------------------
!  fluid thermal property info
!---------------------------------------------------------------------------------------------------------- 
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
!----------------------------------------------------------------------------------------------------------
!  parameters to calculate the fluid thermal property 
!---------------------------------------------------------------------------------------------------------- 
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
    real(WP) :: dhmax
    real(WP) :: dhmin
    type(t_fluidThermoProperty) :: ftp0ref    ! dim, reference state
    type(t_fluidThermoProperty) :: ftpini     ! undim, initial state
  end type t_fluid_parameter
!----------------------------------------------------------------------------------------------------------
!  domain info
!---------------------------------------------------------------------------------------------------------- 
  type t_domain
    logical :: is_periodic(NDIM)        ! is this direction periodic bc?
    logical :: is_stretching(NDIM)      ! is this direction of stretching grids?
    logical :: is_compact_scheme     ! is compact scheme applied?
    logical :: is_thermo             ! is thermal field considered? 
    integer :: idom                  ! domain id
    integer :: icase                 ! case id
    integer :: icoordinate           ! coordinate type
    
    integer :: icht
    integer :: iTimeScheme
    integer :: iviscous
    integer :: iAccuracy
    integer :: ckpt_nfre
    integer :: visu_nfre
    integer :: visu_idim
    integer :: visu_nskip(NDIM)
    integer :: stat_istart
    integer :: stat_nskip(NDIM)
    integer :: nsubitr
    integer :: istret
    integer :: nc(NDIM) ! geometric cell number
    integer :: np_geo(NDIM) ! geometric points
    integer :: np(NDIM) ! calculated points
    integer :: proben   ! global number of probed points
    integer  :: ibcx(2, NBC) ! real bc type, (5 variables, 2 sides), u, v, w, p, T
    integer  :: ibcy(2, NBC) ! real bc type, (5 variables, 2 sides)
    integer  :: ibcz(2, NBC) ! real bc type, (5 variables, 2 sides)
    integer  :: ibcx_nominal(2, NBC) ! nominal (given) bc type, (5 variables, 2 sides), u, v, w, p, T
    integer  :: ibcy_nominal(2, NBC) ! nominal (given) bc type, (5 variables, 2 sides)
    integer  :: ibcz_nominal(2, NBC) ! nominal (given) bc type, (5 variables, 2 sides)
    real(wp) :: fbcx_const(2, NBC) ! bc values, (5 variables, 2 sides)
    real(wp) :: fbcy_const(2, NBC) ! bc values, (5 variables, 2 sides)
    real(wp) :: fbcz_const(2, NBC) ! bc values, (5 variables, 2 sides)

    real(wp) :: lxx
    real(wp) :: lyt
    real(wp) :: lyb
    real(wp) :: lzz
    real(WP) :: rstret
    real(wp) :: dt

    real(wp) :: h(NDIM) ! uniform dx
    real(wp) :: h1r(NDIM) ! uniform (dx)^(-1)
    real(wp) :: h2r(NDIM) ! uniform (dx)^(-2)
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
    ! cell centre location, mapping
    real(wp), allocatable :: yMappingcc(:, :) ! first coefficient in first deriviation. 1/h'
                                              ! first coefficient in second deriviation 1/h'^2
                                              ! second coefficient in second deriviation -h"/h'^3
    real(wp), allocatable :: yp(:)
    real(wp), allocatable :: yc(:)
    real(wp), allocatable :: fbcx_var(:, :, :, :) ! variable bc
    real(wp), allocatable :: fbcy_var(:, :, :, :) ! variable bc
    real(wp), allocatable :: fbcz_var(:, :, :, :) ! variable bc
    type(t_fluidThermoProperty), allocatable :: ftpbcx_var(:, :, :)  ! undim, xbc state
    type(t_fluidThermoProperty), allocatable :: ftpbcy_var(:, :, :)  ! undim, ybc state
    type(t_fluidThermoProperty), allocatable :: ftpbcz_var(:, :, :)  ! undim, zbc state
    real(WP), allocatable :: probexyz(:, :) ! (1:3, xyz coord)
    logical,  allocatable :: probe_is_in(:)
    integer,  allocatable :: probexid(:, :) ! (1:3, local index)
  end type t_domain
!----------------------------------------------------------------------------------------------------------
!  flow info
!---------------------------------------------------------------------------------------------------------- 
  type t_flow
    integer  :: idriven
    integer  :: igravity
    integer  :: inittype
    integer  :: iterfrom
    integer  :: initReTo
    integer  :: nIterFlowStart
    integer  :: nIterFlowEnd
    integer  :: iteration

    real(WP) :: time
    real(WP) :: ren
    real(WP) :: rre
    real(WP) :: init_velo3d(NDIM)
    real(wp) :: reninit
    real(WP) :: drvfc
    real(WP) :: fgravity(NDIM)

    real(wp) :: noiselevel
  
    real(WP), allocatable :: qx(:, :, :)  !
    real(WP), allocatable :: qy(:, :, :)
    real(WP), allocatable :: qz(:, :, :)
    real(WP), allocatable :: gx(:, :, :)
    real(WP), allocatable :: gy(:, :, :)
    real(WP), allocatable :: gz(:, :, :)

    real(WP), allocatable :: pres(:, :, :)
    real(WP), allocatable :: pcor(:, :, :)
    real(WP), allocatable :: pcor_zpencil_ggg(:, :, :)

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

    real(WP), allocatable :: u_vector_mean(:, :, :, :) ! u, v, w
    real(WP), allocatable :: pr_mean(:, :, :)
    real(WP), allocatable :: uu_tensor6_mean(:, :, :, :) ! uu, vv, ww, uv, uw, vw

  end type t_flow


  type t_thermo
    integer :: ifluid
    integer  :: inittype
    integer  :: iterfrom
    integer  :: iteration
    integer  :: nIterThermoStart
    integer  :: nIterThermoEnd
    real(WP) :: ref_l0
    real(WP) :: ref_T0 ! '0' means dimensional 
    real(WP) :: init_T0
    real(WP) :: time
    real(WP) :: rPrRen
    
    real(WP), allocatable :: dh(:, :, :)
    real(WP), allocatable :: hEnth(:, :, :)
    real(WP), allocatable :: kCond(:, :, :)
    real(WP), allocatable :: tTemp(:, :, :)
    real(WP), allocatable :: ene_rhs(:, :, :)  ! current step rhs
    real(WP), allocatable :: ene_rhs0(:, :, :) ! last step rhs

    real(WP), allocatable :: t_mean(:, :, :)
    real(WP), allocatable :: tt_mean(:, :, :)

  end type t_thermo


end module
!==========================================================================================================
!==========================================================================================================
module vars_df_mod
  use udf_type_mod
  implicit none

  type(t_domain), allocatable, save :: domain(:)
  type(t_flow),   allocatable, save :: flow(:)
  type(t_thermo), allocatable, save :: thermo(:)
end module
!==========================================================================================================
module files_io_mod
  implicit none
  character(9) :: dir_data='1_data'
  character(6) :: dir_visu='2_visu'
  character(9) :: dir_moni='3_monitor'
  character(9) :: dir_chkp='4_check'
  public :: create_directory
contains
  subroutine create_directory
    implicit none
    call system('mkdir -p '//dir_data)
    call system('mkdir -p '//dir_visu)
    call system('mkdir -p '//dir_moni)
    call system('mkdir -p '//dir_chkp)
    return
  end subroutine
end module
!==========================================================================================================
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
    d = abs ( r )
  end function

  elemental function abs_cdp ( r ) result (d)
  COMPLEX(kind = D15P), intent(in) :: r
  real(kind = D15P) :: d
    d = abs ( r ) 
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
!==========================================================================================================
!==========================================================================================================
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

module EvenOdd_mod
  implicit none
contains
  logical function is_even(number)  
    implicit none
    integer, intent(in) :: number  
    ! Check if the number is even or odd
    if (mod(number, 2) == 0) then
        is_even = .true.
    else
        is_even = .false.
    end if
  end function
end module

!==========================================================================================================
module flatten_index_mod
 implicit none 
 
 interface flatten_index
   module procedure flatten_3d_to_1d
   module procedure flatten_2d_to_1d
 end interface
 
contains

 function flatten_3d_to_1d(i, j, k, Nx, Ny, Nz) result(n)
   integer, intent(in) :: i, j, k, Nx, Ny, Nz
   integer :: n
   n = i + Nx * (j - 1)  + Nx * Ny * (k - 1)
 end function
 
 function flatten_2d_to_1d(i, j, Nx, Ny) result(n)
   integer, intent(in) :: i, j, Nx, Ny
   integer :: n
   n = i + Nx * (j - 1)
 end function
 
end module

