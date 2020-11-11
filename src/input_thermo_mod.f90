module input_thermo_mod
  use precision_mod
  use input_general_mod, only: ifluid, t0ref, tiref
  implicit none

  integer, parameter :: ISCP_WATER      = 1, &
                        ISCP_CO2        = 2, &
                        ILIQUID_SODIUM  = 3, &
                        ILIQUID_LEAD    = 4, &
                        ILIQUID_BISMUTH = 5, &
                        ILIQUID_LBE     = 6

  integer, parameter :: IPROPERTY_TABLE = 1, &
                        IPROPERTY_FUNCS = 2

  character(len = 64), parameter :: INPUT_SCP_WATER = 'NIST_WATER_23.5MP.DAT'
  character(len = 64), parameter :: INPUT_SCP_CO2   = 'NIST_CO2_8MP.DAT'
  character(len = 64) :: inputProperty

  integer :: ipropertyState

  real(WP), parameter :: TM0_Na = 371.0 ! unit: K, melting temperature at 1 atm for Na
  real(WP), parameter :: TM0_Pb = 600.6 ! unit: K, melting temperature at 1 atm for Lead
  real(WP), parameter :: TM0_BI = 544.6 ! unit: K, melting temperature at 1 atm for Bismuth
  real(WP), parameter :: TM0_LBE = 398.0 ! unit: K, melting temperature at 1 atm for LBE
  real(WP) :: TM0

  real(WP), parameter :: TB0_Na = 1155.0 ! unit: K, boling temperature at 1 atm for Na
  real(WP), parameter :: TB0_Pb = 2021.0 ! unit: K, boling temperature at 1 atm for Lead
  real(WP), parameter :: TB0_BI = 1831.0 ! unit: K, boling temperature at 1 atm for Bismuth
  real(WP), parameter :: TB0_LBE = 1927.0 ! unit: K, boling temperature at 1 atm for LBE
  real(WP) :: TB0

  real(WP), parameter :: HM0_Na = 113.0e3 ! unit: J / Kg, latent melting heat, enthalpy
  real(WP), parameter :: HM0_Pb = 23.07e3 ! unit: J / Kg, latent melting heat, enthalpy
  real(WP), parameter :: HM0_BI = 53.3e3 ! unit: J / Kg, latent melting heat, enthalpy
  real(WP), parameter :: HM0_LBE = 38.6e3 ! unit: J / Kg, latent melting heat, enthalpy
  real(WP) :: HM0
  ! D = CoD(0) + CoD(1) * T
  real(WP), parameter :: CoD_Na(0:1) = (/1014.0, -0.235/)
  real(WP), parameter :: CoD_Pb(0:1) = (/11441.0, -1.2795/)
  real(WP), parameter :: CoD_Bi(0:1) = (/10725.0, -1.22 /)
  real(WP), parameter :: CoD_LBE(0:1) = (/11065.0, 1.293 /)
  real(WP) :: CoD(0:1)
  ! K = CoK(0) + CoK(1) * T + CoK(2) * T^2
  real(WP), parameter :: CoK_Na(0:2) = (/104.0, -0.047, 0.0/)
  real(WP), parameter :: CoK_Pb(0:2) = (/9.2, 0.011, 0.0/)
  real(WP), parameter :: CoK_Bi(0:2) = (/7.34, 9.5E-3, 0.0/)
  real(WP), parameter :: CoK_LBE(0:2) = (/ 3.284, 1.617E-2, -2.305E-6/)
  real(WP) :: CoK(0:2)
  ! B = 1 / (CoB - T)
  real(WP), parameter :: CoB_Na = 4316.0
  real(WP), parameter :: CoB_Pb = 8942.0
  real(WP), parameter :: CoB_BI = 8791.0
  real(WP), parameter :: CoB_LBE = 8558.0
  real(WP) :: CoB
  ! Cp = CoCp(-2) * T^(-2) + CoCp(-1) * T^(-1) + CoCp(0) + CoCp(1) * T + CoCp(2) * T^2
  real(WP), parameter :: CoCp_Na(-2:2) = (/- 3.001e6, 0.0, 1658.0, -0.8479, 4.454E-4/)
  real(WP), parameter :: CoCp_Pb(-2:2) = (/- 1.524e6, 0.0, 176.2, -4.923E-2, 1.544E-5/)
  real(WP), parameter :: CoCp_Bi(-2:2) = (/7.183e6, 0.0, 118.2, 5.934E-3, 0.0/)
  real(WP), parameter :: CoCp_LBE(-2:2) = (/-4.56e5, 0.0, 164.8, - 3.94E-2, 1.25E-5/)
  real(WP) :: CoCp(-2:2)
  ! H = HM0 + CoH(-1) * (1 / T - 1 / Tm0) + CoH(0) + CoH(1) * (T - Tm0) +  CoH(2) * (T^2 - Tm0^2) +  CoH(3) * (T^3- Tm0^3)
  real(WP), parameter :: CoH_Na(-1:3) = (/4.56e5, 0.0, 164.8, -1.97E-2, 4.167E-4/)
  real(WP), parameter :: CoH_Pb(-1:3) = (/1.524e6, 0.0, 176.2, -2.4615E-2, 5.147E-6/)
  real(WP), parameter :: CoH_Bi(-1:3) = (/-7.183e6, 0.0, 118.2, 2.967E-3, 0.0/)
  real(WP), parameter :: CoH_LBE(-1:3) = (/4.56e5, 0.0, 164.8, -1.97E-2, 4.167E-4/)! check, WRong from literature.
  real(WP) :: CoH(-1:3)
  ! M = vARies
  real(WP), parameter :: CoM_Na(-1:1) = (/556.835, -6.4406, -0.3958/) ! M = exp ( CoM(-1) / T + CoM(0) + CoM(1) * ln(T) )
  real(WP), parameter :: CoM_Pb(-1:1) = (/1069.0, 4.55E-4, 0.0/) ! M = CoM(0) * exp (CoM(-1) / T)
  real(WP), parameter :: CoM_Bi(-1:1) = (/780.0, 4.456E-4, 0.0/) ! M = CoM(0) * exp (CoM(-1) / T)
  real(WP), parameter :: CoM_LBE(-1:1) = (/754.1, 4.94E-4, 0.0/) ! M = CoM(0) * exp (CoM(-1) / T)
  real(WP) :: CoM(-1:1)

  type thermoProperty_t
    real(WP) :: t !temperature
    real(WP) :: d !density
    real(WP) :: m !dynviscosity
    real(WP) :: k !thermconductivity
    real(WP) :: h !enthalpy
    real(WP) :: dh ! mass enthalpy
    real(WP) :: cp ! specific heat capacity 
    real(WP) :: b ! thermal expansion

    contains

    private

    procedure, public :: is_T_in_scope
    procedure, public :: Get_initilized_tp
    procedure, public :: GetTP_from_T
    procedure, public :: GetTP_from_H
    procedure, public :: GetTP_from_DH

    procedure :: Print_debug
    generic :: Print => Print_debug
    generic :: write(formatted) => Print_debug
  end type thermoProperty_t

  integer :: nlist

  type(thermoProperty_t), save, allocatable, dimension(:) :: listTP
  type(thermoProperty_t) :: tpRef0 ! dim
  type(thermoProperty_t) :: tpIni0 ! dim


  private :: Sort_listTP_Tsmall2big
  private :: GetTP_from_list_T
  private :: GetTP_from_function_T
  private :: Check_monotonicity_DH_of_HT_list
  private :: Initialize_thermo_parameters
  private :: Buildup_property_relations_from_table
  private :: Buildup_property_relations_from_function

  public  :: Initialize_thermo_input

contains
  !!--------------------------------
  subroutine Get_initilized_tp ( tp )
    use parameters_constant_mod, only: ZERO, ONE
    class(thermoProperty_t) :: tp

    tp%t = ONE
    tp%d = ONE
    tp%m = ONE
    tp%k = ONE
    tp%cp = ONE
    tp%b = ONE
    tp%h = ZERO
    tp%dh = ZERO

  end subroutine Get_initilized_tp


  subroutine is_T_in_scope ( tp )
    class(thermoProperty_t) :: tp

    if(ipropertyState == IPROPERTY_TABLE) then
      if ( ( tp%t < listTP(1)%t     )  .OR. &
           ( tp%t > listTP(nlist)%t ) ) then
        print*, tp%t, listTP(1)%t, listTP(nlist)%t
        stop 'temperature exceeds specified range.'
      end if
    end if

    if(ipropertyState == IPROPERTY_FUNCS) then
      if ( ( tp%t < ( TM0 / tpRef0%t ) ) .OR. &
           ( tp%t > ( TB0 / tpRef0%t ) ) ) then 
        print*, tp%t, TM0 / tpRef0%t, TB0 / tpRef0%t
        stop 'temperature exceeds specified range.'
      end if
    end if

  end subroutine is_T_in_scope

  !!--------------------------------
  subroutine Sort_listTP_Tsmall2big(list)
    type(thermoProperty_t),intent(inout) :: list(:)
    integer :: i, n, k
    real(WP) :: buf

    n = size( list )

    do i = 1, n
      k = minloc( list(i:n)%t, dim = 1) + i - 1

      buf = list(i)%t
      list(i)%t = list(k)%t
      list(k)%t = buf

      buf = list(i)%d
      list(i)%d = list(k)%d
      list(k)%d = buf

      buf = list(i)%m
      list(i)%m = list(k)%m
      list(k)%m = buf

      buf = list(i)%k
      list(i)%k = list(k)%k
      list(k)%k = buf

      buf = list(i)%b
      list(i)%b = list(k)%b
      list(k)%b = buf

      buf = list(i)%cp
      list(i)%cp = list(k)%cp
      list(k)%cp = buf

      buf = list(i)%h
      list(i)%h = list(k)%h
      list(k)%h = buf

      buf = list(i)%dh
      list(i)%dh = list(k)%dh
      list(k)%dh = buf

    end do
  end subroutine Sort_listTP_Tsmall2big
  !!--------------------------------
  subroutine GetTP_from_H(a)
    use parameters_constant_mod, only : MINP, ONE
    class(thermoProperty_t), intent(inout) :: a

    integer :: i1, i2, im
    real(WP) :: d1, dm
    real(WP) :: w1, w2

    i1 = 1
    i2 = nlist

    do while ( (i2 - i1) > 1)
      im = i1 + (i2 - i1) / 2
      d1 = listTP(i1)%h - a%h
      dm = listTP(im)%h - a%h
      if ( (d1 * dm) > MINP ) then
        i1 = im
      else
        i2 = im
      end if
    end do

    w1 = (listTP(i2)%h - a%h) / (listTP(i2)%h - listTP(i1)%h) 
    w2 = ONE - w1

    a%d = w1 * listTP(i1)%d + w2 * listTP(i2)%d
    a%m = w1 * listTP(i1)%m + w2 * listTP(i2)%m
    a%k = w1 * listTP(i1)%k + w2 * listTP(i2)%k
    a%t = w1 * listTP(i1)%t + w2 * listTP(i2)%t
    a%b = w1 * listTP(i1)%b + w2 * listTP(i2)%b
    a%cp = w1 * listTP(i1)%cp + w2 * listTP(i2)%cp
    a%dh = a%d * a%h

  end subroutine GetTP_from_H

  !!--------------------------------
  subroutine GetTP_from_DH(a)
    use parameters_constant_mod, only : MINP, ONE
    class(thermoProperty_t), intent(inout) :: a

    integer :: i1, i2, im
    real(WP) :: d1, dm
    real(WP) :: w1, w2

    i1 = 1
    i2 = nlist

    do while ( (i2 - i1) > 1)
      im = i1 + (i2 - i1) / 2
      d1 = listTP(i1)%dh - a%dh
      dm = listTP(im)%dh - a%dh
      if ( (d1 * dm) > MINP ) then
        i1 = im
      else
        i2 = im
      end if
    end do

    w1 = (listTP(i2)%dh - a%dh) / (listTP(i2)%dh - listTP(i1)%dh) 
    w2 = ONE - w1
    
    if(ipropertyState == IPROPERTY_TABLE) then 
      a%h = w1 * listTP(i1)%h + w2 * listTP(i2)%h
      call a%GetTP_from_H
    else if (ipropertyState == IPROPERTY_FUNCS) then 
      a%t = w1 * listTP(i1)%t + w2 * listTP(i2)%t
      call a%GetTP_from_T
    else  
      STOP 'Error. No such option of ipropertyState.'
    end if
  end subroutine GetTP_from_DH
  !!--------------------------------
  subroutine GetTP_from_list_T(a)
    use parameters_constant_mod, only : MINP, ONE
    type(thermoProperty_t), intent(inout) :: a

    integer :: i1, i2, im
    real(WP) :: d1, dm
    real(WP) :: w1, w2

    i1 = 1
    i2 = nlist

    do while ( (i2 - i1) > 1)
      im = i1 + (i2 - i1) / 2
      d1 = listTP(i1)%t - a%t
      dm = listTP(im)%t - a%t
      if ( (d1 * dm) > MINP ) then
        i1 = im
      else
        i2 = im
      end if
    end do

    w1 = (listTP(i2)%t - a%t) / (listTP(i2)%t - listTP(i1)%t) 
    w2 = ONE - w1

    a%d = w1 * listTP(i1)%d + w2 * listTP(i2)%d
    a%m = w1 * listTP(i1)%m + w2 * listTP(i2)%m
    a%k = w1 * listTP(i1)%k + w2 * listTP(i2)%k
    a%h = w1 * listTP(i1)%h + w2 * listTP(i2)%h
    a%b = w1 * listTP(i1)%b + w2 * listTP(i2)%b
    a%cp = w1 * listTP(i1)%cp + w2 * listTP(i2)%cp
    a%dh = a%d * a%h

  end subroutine GetTP_from_list_T

  !!--------------------------------
  subroutine GetTP_from_function_T(a, dim)
    use parameters_constant_mod, only: ONE
    type(thermoProperty_t) :: a
    integer, intent(in), optional :: dim ! without = undim, with = dim.
    

    real(WP) :: t1, dummy

    if (present(dim)) then 
      t1 = a%t
    else 
      ! convert undim to dim 
      t1 = a%t * tpRef0%t
    end if

    ! D = density = f(T)
    dummy = CoD(0) + CoD(1) * t1
    if (present(dim)) then 
      a%d = dummy
    else 
      a%d = dummy / tpRef0%d
    end if

    ! K = thermal conductivity = f(T)
    dummy = CoK(0) + CoK(1) * t1 + CoK(2) * t1**2
    if (present(dim)) then 
      a%k = dummy
    else 
      a%k = dummy / tpRef0%k
    end if

    ! Cp = f(T)
    dummy = CoCp(-2) * t1**(-2) + CoCp(-1) * t1**(-1) + CoCp(0) + CoCp(1) * t1 + CoCp(2) * t1**2
    if (present(dim)) then 
      a%cp = dummy
    else 
      a%cp = dummy / tpRef0%cp
    end if

    ! H = entropy = f(T)
    dummy = Hm0 + &
      CoH(-1) * (ONE / t1 - ONE / Tm0) + &
      CoH(0) + &
      CoH(1) * (t1 - Tm0) + &
      CoH(2) * (t1**2 - Tm0**2) + &
      CoH(3) * (t1**3 - Tm0**3)
    if (present(dim)) then 
      a%h = dummy
    else 
      a%h = (dummy - tpRef0%h) / (tpRef0%cp * tpRef0%t)
    end if

    ! B = f(T)
    dummy = ONE / (CoB - t1)
    if (present(dim)) then 
      a%b = dummy
    else 
      a%b = dummy / tpRef0%b
    end if

    ! dynamic viscosity = f(T)
    select case (ifluid)
      ! unit: T(Kelvin), M(Pa S)
      case (ILIQUID_SODIUM)
        dummy = EXP( CoM_Na(-1) / t1 + CoM_Na(0) + CoM_Na(1) * LOG(t1) )
      case (ILIQUID_LEAD)
        dummy = CoM_Pb(0) * EXP (CoM_Pb(-1) / t1)
      case (ILIQUID_BISMUTH)
        dummy = CoM_Bi(0) * EXP (CoM_Bi(-1) / t1)
      case (ILIQUID_LBE)
        dummy = CoM_LBE(0) * EXP (CoM_LBE(-1) / t1)
      case default
        dummy = EXP( CoM_Na(-1) / t1 + CoM_Na(0) + CoM_Na(1) * LOG(t1) )
    end select
    if (present(dim)) then 
      a%m = dummy
    else 
      a%m = dummy / tpRef0%m
    end if


    a%dh = a%d * a%h

  end subroutine GetTP_from_function_T

  !!--------------------------------
  subroutine GetTP_from_T ( tp )
    class(thermoProperty_t) :: tp

    if(ipropertyState == IPROPERTY_TABLE) then 
      call GetTP_from_list_T(tp)
    end if 

    if(ipropertyState == IPROPERTY_FUNCS) then 
      call GetTP_from_function_T(tp)
    end if
  end subroutine GetTP_from_T
  
  !!--------------------------------
  subroutine Initialize_thermo_parameters

    select case (ifluid)
    case (ISCP_WATER)
      ipropertyState = IPROPERTY_TABLE
      inputProperty = TRIM(INPUT_SCP_WATER)

    case (ISCP_CO2)
      ipropertyState = IPROPERTY_TABLE
      inputProperty = TRIM(INPUT_SCP_CO2)

    case (ILIQUID_SODIUM)
      ipropertyState = IPROPERTY_FUNCS
      TM0 = TM0_Na
      TB0 = TB0_Na
      HM0 = HM0_Na
      CoD(0:1) = CoD_Na(0:1)
      CoK(0:2) = CoK_Na(0:2)
      CoB = CoB_Na
      CoCp(-2:2) = CoCp_Na(-2:2)
      CoH(-1:3) = CoH_Na(-1:3)
      CoM(-1:1) = CoM_Na(-1:1)

    case (ILIQUID_LEAD)
      ipropertyState = IPROPERTY_FUNCS
      TM0 = TM0_Pb
      TB0 = TB0_Pb
      HM0 = HM0_Pb
      CoD(0:1) = CoD_Pb(0:1)
      CoK(0:2) = CoK_Pb(0:2)
      CoB = CoB_Pb
      CoCp(-2:2) = CoCp_Pb(-2:2)
      CoH(-1:3) = CoH_Pb(-1:3)
      CoM(-1:1) = CoM_Pb(-1:1)

    case (ILIQUID_BISMUTH)
      ipropertyState = IPROPERTY_FUNCS
      TM0 = TM0_BI
      TB0 = TB0_BI
      HM0 = HM0_BI
      CoD(0:1) = CoD_BI(0:1)
      CoK(0:2) = CoK_BI(0:2)
      CoB = CoB_BI
      CoCp(-2:2) = CoCp_BI(-2:2)
      CoH(-1:3) = CoH_BI(-1:3)
      CoM(-1:1) = CoM_BI(-1:1)

    case (ILIQUID_LBE)
      ipropertyState = IPROPERTY_FUNCS
      TM0 = TM0_LBE
      TB0 = TB0_LBE
      HM0 = HM0_LBE
      CoD(0:1) = CoD_LBE(0:1)
      CoK(0:2) = CoK_LBE(0:2)
      CoB = CoB_LBE
      CoCp(-2:2) = CoCp_LBE(-2:2)
      CoH(-1:3) = CoH_LBE(-1:3)
      CoM(-1:1) = CoM_LBE(-1:1)
    case default
      ipropertyState = IPROPERTY_FUNCS
      TM0 = TM0_Na
      TB0 = TB0_Na
      HM0 = HM0_Na
      CoD(0:1) = CoD_Na(0:1)
      CoK(0:2) = CoK_Na(0:2)
      CoB = CoB_Na
      CoCp(-2:2) = CoCp_Na(-2:2)
      CoH(-1:3) = CoH_Na(-1:3)
      CoM(-1:1) = CoM_Na(-1:1)

    end select

    call tpRef0%Get_initilized_tp()
    call tpIni0%Get_initilized_tp()
    tpRef0%t = t0Ref
    tpIni0%t = tiRef

  end subroutine Initialize_thermo_parameters

  !!--------------------------------
  subroutine Buildup_property_relations_from_table ( )
    use mpi_mod
    use iso_fortran_env, only : ERROR_UNIT, IOSTAT_END
    use parameters_constant_mod, only : ZERO

    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit
    character(len = 80) :: str
    real(WP) :: rtmp
    integer :: i
 
    open ( newunit = inputUnit, file = inputProperty, status = 'old', action  = 'read', &
          iostat = ioerr, iomsg = iotxt)
    if(ioerr /= 0) then
      write (ERROR_UNIT, *) 'Problem openning : ', inputProperty, ' for reading.'
      write (ERROR_UNIT, *) 'Message: ', trim (iotxt)
      stop 4
    end if

    nlist = 0
    read(inputUnit, *, iostat = ioerr) str
    do
      read(inputUnit, *, iostat = ioerr) rtmp, rtmp, rtmp, rtmp, rtmp, rtmp, rtmp, rtmp
      if(ioerr /= 0) exit
      nlist = nlist + 1
    end do
    rewind(inputUnit)

    allocate ( listTP (nlist) )

    read(inputUnit, *, iostat = ioerr) str
    block_tablereading: do i = 1, nlist
     call listTP(i)%Get_initilized_tp()
      read(inputUnit, *, iostat = ioerr) rtmp, listTP(i)%h, listTP(i)%t, listTP(i)%d, &
      listTP(i)%m, listTP(i)%k, listTP(i)%cp, listTP(i)%b
      listTP(i)%dh = listTP(i)%d * listTP(i)%h
    end do block_tablereading

    call Sort_listTP_Tsmall2big (listTP(:))

    call GetTP_from_list_T(tpRef0)
    call GetTP_from_list_T(tpIni0)

    listTP(:)%t = listTP(:)%t / tpRef0%t
    listTP(:)%d = listTP(:)%d / tpRef0%d
    listTP(:)%m = listTP(:)%m / tpRef0%m
    listTP(:)%k = listTP(:)%k / tpRef0%k
    listTP(:)%b = listTP(:)%b / tpRef0%b
    listTP(:)%cp = listTP(:)%cp / tpRef0%cp
    listTP(:)%h = (listTP(:)%h - tpRef0%h) / tpRef0%t / tpRef0%cp
    listTP(:)%dh = listTP(:)%d * listTP(:)%h

  end subroutine Buildup_property_relations_from_table

  !!--------------------------------
  subroutine Buildup_property_relations_from_function ( )
    integer :: i

    call GetTP_from_function_T(tpRef0, 1)
    call GetTP_from_function_T(tpIni0, 1)
    
    nlist = 1024
    allocate ( listTP (nlist) )
    
    do i = 1, nlist
      call listTP(i)%Get_initilized_tp()
      listTP(i)%t = ( Tm0 + (Tb0 - Tm0) * real(i, WP) / real(nlist, WP) ) / tpRef0%t
      call GetTP_from_function_T( listTP(i) )
    end do

  end subroutine Buildup_property_relations_from_function

  !!--------------------------------
  subroutine Print_debug(this, unit, iotype, v_list, iostat, iomsg)
    use iso_fortran_env, only : error_unit
    class(thermoProperty_t), intent(in) :: this
    integer, intent(in) :: unit
    character(len = *), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(len = *), intent(inout) :: iomsg
    integer :: i_pass

    iostat = 0
    iomsg = ""
    
    this_block: do i_pass = 1, 1
      !write(unit, *, iostat = iostat, iomsg = iomsg) 'thermalProperty'
      !if(iostat /= 0) exit this_block
      if(iotype(1:2) == 'DT' .and. len(iotype) > 2) &
        write(unit, *, iostat = iostat, iomsg = iomsg) iotype(3:)
      if(iostat /= 0) exit this_block

      write(unit, *, iostat = iostat, iomsg = iomsg) &
      this%h, this%t, this%d, this%m, this%k, this%cp, this%b, this%dh
      
      if(iostat /= 0) exit this_block
    end do this_block

    if(iostat /= 0) then
      write (error_unit, "(A)") "print error : " // trim(iomsg)
      write (error_unit, "(A, I0)") "  iostat : ", iostat
    end if
  end subroutine Print_debug

  subroutine Check_monotonicity_DH_of_HT_list
    use parameters_constant_mod, only : MINP
    integer :: i
    real(WP) :: ddh1, dt1, dh1
    real(WP) :: ddh2, dt2, dh2
    real(WP) :: ddt, ddh

    do i = 2, nlist - 1
        ddh1 = listTP(i)%dh - listTP(i - 1)%dh
        dt1 = listTP(i)%t - listTP(i - 1)%t
        dh1 = listTP(i)%h - listTP(i - 1)%h

        ddh2 = listTP(i + 1)%dh - listTP(i)%dh
        dt2 = listTP(i + 1)%t - listTP(i)%t
        dh2 = listTP(i + 1)%h - listTP(i)%h

        ddt = ddh1 / dt1 * ddh2 / dt2 
        ddh = ddh1 / dh1 * ddh2 / dh2

        if (ddt < MINP) STOP 'Error. The relation (rho * h) = FUNCTION (T) is not monotonicity.'
        if (ddh < MINP) STOP 'Error. The relation (rho * h) = FUNCTION (H) is not monotonicity.'

    end do

  end subroutine


  subroutine Write_thermo_property
    use mpi_mod ! for test
    use parameters_constant_mod, only : ZERO, TRUNCERR
    type(thermoProperty_t) :: tp
    integer :: n, i
    real(WP) :: dhmax1, dhmin1
    real(WP) :: dhmax, dhmin
    integer :: tp_unit

    if (nrank /= 0) return

    n = 128
    call tp%Get_initilized_tp

    open (newunit = tp_unit, file = 'check_tp_from_dh.dat')
    

    dhmin = ZERO
    dhmax = ZERO
    if(ipropertyState == IPROPERTY_TABLE) then 
      dhmin = listTP(1)%dh + TRUNCERR
      dhmax = listTP(nlist)%dh - TRUNCERR
    end if

    if(ipropertyState == IPROPERTY_FUNCS) then 
      tp%t  = TB0 / tpRef0%t
      call tp%GetTP_from_T
      dhmin1 = tp%dh

      tp%t  = TM0 / tpRef0%t
      call tp%GetTP_from_T
      dhmax1 = tp%dh
      
      dhmin = dmin1( dhmin1, dhmax1) + TRUNCERR
      dhmax = dmax1( dhmin1, dhmax1) - TRUNCERR
    end if

    do i = 1, n
      tp%dh = dhmin + (dhmax - dhmin) * real(i - 1, WP) / real(n - 1, WP)
      call tp%GetTP_from_DH
      call tp%is_T_in_scope
      write(tp_unit, '(dt)') tp
    end do
    close (tp_unit)

  end subroutine 

  !!--------------------------------
  subroutine Initialize_thermo_input
    use input_general_mod, only : ithermo
    
    if (ithermo /= 1) return 
    call Initialize_thermo_parameters
    if (ipropertyState == IPROPERTY_TABLE) call Buildup_property_relations_from_table
    if (ipropertyState == IPROPERTY_FUNCS) call Buildup_property_relations_from_function
    call Check_monotonicity_DH_of_HT_list
    call Write_thermo_property ! for test

  end subroutine Initialize_thermo_input



end module input_thermo_mod



